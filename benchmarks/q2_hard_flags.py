"""Round 6, Q2 (development): flags for hard-region errors from the sample's own reads.

`docs/design/q2-hard-flags.md` fixes the features, the flag rule, the tuning
and the gate. Per platform, depth and chromosome a subprocess (with that
platform's environment) reads the stored fragments and phaser VCFs and writes
per-variant depth and allele balance, and per phased pair its support, its
error against the truth and its window MEC, to `data/r6_q2/rows/`. The k-mer
stage (`--kmers`) counts flanking reference 31-mers of every HG002 variant in
the HG002 HiFi reads and in the reference with jellyfish. The main run tunes
flags on odd chromosomes, scores them on even ones, and compares Learn then
Test with and without the flagged pairs.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/q2_hard_flags.py [--kmers]`
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from itertools import combinations
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import genome_common  # noqa: E402

ROOT = genome_common.ROOT  # the repository, wherever this file is run from
OUT = ROOT / "data" / "r6_q2"
CHROMS = tuple(f"chr{i}" for i in range(1, 23) if i != 20)
PLATFORMS = {"hprc_HG002": ("30", "20", "10"), "w4_grch38": ("30",),
             "cal_v5": ("30", "20", "10"), "hg001_v5": ("30", "20", "10")}
GATE, SOURCE, TRANSFER = "hprc_HG002", "cal_v5", "hg001_v5"
PHASERS = ("whatshap", "longphase")
HARD = ("segdup", "lowmap", "tandem", "difficult", "indel")
FEATURES = ("depth", "balance", "mec", "kmer")
WINDOW = 2
ABSTAIN_TUNE = 0.09
ABSTAIN_LIMIT = 0.10
RECALL = 0.5
MIN_CONFIDENT = 10
MIN_READS = 10
K = 31
KMER_PLATFORMS = ("hprc_HG002", "w4_grch38")
JELLYFISH = Path.home() / ".local" / "opt" / "jellyfish" / "jellyfish-2.3.0"
READS = ROOT / "data" / "w4" / "reads.fq.gz"
GENOME = ROOT / "data" / "reference" / "GRCh38_no_alt_analysis_set.fna"
WORKERS = 14  # with the k-mer stage's 8 jellyfish threads and pigz, within 24 cores


# Pure parts, checked in tests/test_q2_hard_flags.py.

def window_columns(a: int, b: int, n: int, w: int = WINDOW) -> list[int]:
    """The pair's columns `a < b` and `w` columns on each side outside `[a, b]`, within `[0, n)`."""
    return list(range(max(0, a - w), a)) + [a, b] + list(range(b + 1, min(n, b + 1 + w)))


@lru_cache(maxsize=None)
def _haplotypes(k: int) -> np.ndarray:
    """Every haplotype over `k` columns with the first allele 0, as `exact` enumerates them."""
    codes = np.arange(2 ** (k - 1))
    return np.array([[0] + [(c >> (k - 2 - i)) & 1 for i in range(k - 1)] for c in codes], float)


def window_mec(alleles: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    """`(optimum wMEC, total weight)` over reads with at least two calls in the window.

    `alleles` is reads by columns with 0, 1, or -1 for no call; `weights` the
    same shape. The second haplotype is the complement of the first, and each
    read takes its cheaper side, so this is `exact.solve_bruteforce`'s cost on
    the window's matrix (single-call reads cost nothing there either).
    """
    called = alleles >= 0
    keep = called.sum(axis=1) >= 2
    if not keep.any():
        return 0.0, 0.0
    m = np.where(called, weights, 0.0)[keep]
    x = (alleles[keep] == 1).astype(float)
    h = _haplotypes(alleles.shape[1])
    mx = m * x
    mismatch = mx.sum(axis=1)[:, None] + m @ h.T - 2 * (mx @ h.T)
    total = m.sum(axis=1)[:, None]
    cost = np.minimum(mismatch, total - mismatch).sum(axis=0)
    return float(cost.min()), float(total.sum())


def balance_score(n0: np.ndarray, n1: np.ndarray) -> np.ndarray:
    """`-log10` of the two-sided binomial p-value of the smaller allele count under p = 1/2."""
    from scipy.stats import binom

    n = n0 + n1
    p = np.minimum(1.0, 2 * binom.cdf(np.minimum(n0, n1), n, 0.5))
    return np.where(n > 0, -np.log10(np.maximum(p, 1e-300)), 0.0)


def tail_rank(reference_sorted: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Share of the reference values at least as large as each value (1 for -inf)."""
    n = len(reference_sorted)
    return (n - np.searchsorted(reference_sorted, values, side="left")) / max(1, n)


def choose_tau(m: np.ndarray, share: float) -> float:
    """The largest value `tau` of `m` with at most `share` of `m` at or below it; -inf if none."""
    s = np.sort(m)
    u = np.unique(s)
    count = np.searchsorted(s, u, side="right")
    ok = u[count <= share * len(s)]
    return float(ok.max()) if len(ok) else -math.inf


def kmer_spans(pos1: int, ref: str, k: int = K) -> tuple[tuple[int, int], tuple[int, int]]:
    """0-based half-open spans of the reference k-mer ending before the variant and the one
    starting after its reference allele."""
    start0 = pos1 - 1
    end0 = start0 + len(ref)
    return (start0 - k, start0), (end0, end0 + k)


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def canonical(kmer: str) -> str:
    return min(kmer, kmer.translate(COMPLEMENT)[::-1])


def flat_calls(matrix) -> tuple:
    """`(offsets, column, read, allele, weight)` of every called allele, sorted by column;
    column `c`'s calls are `offsets[c]:offsets[c + 1]`."""
    col, read, allele, weight = [], [], [], []
    for r, fragment in enumerate(matrix):
        for j, a, w in fragment.observed():
            col.append(j)
            read.append(r)
            allele.append(a)
            weight.append(w)
    order = np.argsort(np.array(col, np.int64), kind="stable")
    col = np.array(col, np.int64)[order]
    offsets = np.searchsorted(col, np.arange(matrix.n_variants + 1))
    return (offsets, col, np.array(read, np.int64)[order], np.array(allele, np.int8)[order],
            np.array(weight, float)[order])


def window_arrays(calls: tuple, cols: list[int]) -> tuple[np.ndarray, np.ndarray]:
    """Reads by `cols` alleles (-1 for no call) and weights, for reads calling any of `cols`."""
    offsets, _, read, allele, weight = calls
    parts = [np.arange(offsets[c], offsets[c + 1]) for c in cols]
    where = np.concatenate([np.full(len(p), i, np.int64) for i, p in enumerate(parts)])
    idx = np.concatenate(parts).astype(np.int64)
    reads, row = np.unique(read[idx], return_inverse=True)
    al = np.full((len(reads), len(cols)), -1, np.int8)
    wt = np.zeros((len(reads), len(cols)))
    al[row, where] = allele[idx]
    wt[row, where] = weight[idx]
    return al, wt


# Per chromosome, in the platform's environment.

def build(chrom: str, depth: str) -> dict:
    import genome_common as common
    from phase_confidence import phasing_from_vcf
    from phasebudget.confidence import pair_support

    extractor = common._module("extract_fragments.py")
    sc = common._module("score_hg002.py")
    variants = sc.variant_table(common.hetsnv_path(chrom))
    matrix, positions, _, _ = extractor.load_fragments(
        str(common.fragments_path(chrom, depth, "realign")))
    truth = common.truth_at(chrom, variants, positions)
    n = matrix.n_variants
    cov = np.array(matrix.coverage(), float)
    median = max(1.0, float(np.median(cov)))
    logdepth = np.abs(np.log2(np.maximum(cov, 0.5) / median))
    calls = flat_calls(matrix)
    col, allele = calls[1], calls[3]
    n1 = np.bincount(col, weights=(allele == 1), minlength=n)
    n0 = np.bincount(col, minlength=n) - n1
    balance = balance_score(n0, n1)
    mec: dict = {}

    def window(a: int, b: int) -> float:
        if (a, b) not in mec:
            cost, total = window_mec(*window_arrays(calls, window_columns(a, b, n)))
            mec[(a, b)] = cost / total if total > 0 else -math.inf
        return mec[(a, b)]

    pairs = {}
    for phaser in PHASERS:
        vcf = common.chrom_dir(chrom) / "external" / f"{phaser}_{depth}.vcf"
        vcf = vcf if vcf.exists() else vcf.with_name(vcf.name + ".gz")
        phasing = phasing_from_vcf(str(vcf), common.VCF_SAMPLE, chrom, positions)
        rows = []
        for s in pair_support(matrix, phasing):
            a, b = s.first, s.second
            if truth[a] is None or truth[b] is None:
                continue
            wrong = (phasing.alleles[a] == phasing.alleles[b]) != (truth[a] == truth[b])
            rows.append([a, b, s.agree, s.disagree, bool(wrong), window(a, b)])
        pairs[phaser] = rows
    return {"pos1": [p + 1 for p in positions], "logdepth": logdepth.tolist(),
            "balance": balance.tolist(), "pairs": pairs}


def chrom_path(platform: str, depth: str, chrom: str) -> Path:
    """The cached output of `build`, computed in a subprocess with the platform's environment."""
    path = OUT / "rows" / f"{platform}_{depth}_{chrom}.json"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        env = dict(os.environ, PHASEBUDGET_PLATFORM=platform,
                   PYTHONPATH=f"{ROOT / 'src'}:{ROOT / 'benchmarks'}:{HERE}")
        subprocess.run([sys.executable, __file__, "--build", platform, depth, chrom, str(path)],
                       env=env, check=True)
    return path


# K-mer stage.

def kmers() -> None:
    """Write `data/r6_q2/kmer/excess.json`: per chromosome and variant position, the largest
    log2 copy excess of its two flanking reference k-mers in the HG002 HiFi reads."""
    import gzip

    import pysam

    d = OUT / "kmer"
    d.mkdir(parents=True, exist_ok=True)
    site_kmers: dict = {}
    for chrom in CHROMS:
        fasta = pysam.FastaFile(str(ROOT / "data" / "reference" / f"{chrom}.fa"))
        length = fasta.get_reference_length(chrom)
        sites = {}
        for platform in KMER_PLATFORMS:
            with gzip.open(ROOT / "data" / f"ext_{platform}" / chrom / "hetsnv.vcf.gz", "rt") as h:
                for line in h:
                    if not line.startswith("#"):
                        f = line.split("\t", 5)
                        sites[int(f[1])] = f[3]
        out = {}
        for pos1, ref in sites.items():
            mers = []
            for s, e in kmer_spans(pos1, ref):
                if s >= 0 and e <= length:
                    seq = fasta.fetch(chrom, s, e).upper()
                    if len(seq) == K and set(seq) <= set("ACGT"):
                        mers.append(canonical(seq))
            out[pos1] = mers
        site_kmers[chrom] = out
        print(f"  {chrom}: {len(out)} variants", flush=True)
    every = sorted({m for c in site_kmers.values() for ms in c.values() for m in ms})
    with open(d / "sites.fa", "w") as handle:
        for i, m in enumerate(every):
            handle.write(f">{i}\n{m}\n")
    print(f"  {len(every)} distinct canonical {K}-mers", flush=True)
    size = str(2 * len(every))
    (d / "reads.gen").write_text(f"pigz -dc -p 4 {READS}\n")
    for name, inputs in (("reads", ["-g", str(d / "reads.gen"), "-G", "1"]),
                         ("genome", [str(GENOME)])):
        if not (d / f"{name}.jf").exists():
            subprocess.run([str(JELLYFISH), "count", "-m", str(K), "-C", "-s", size, "-t", "8",
                            "--if", str(d / "sites.fa"), "-o", str(d / f"{name}.part.jf"),
                            *inputs], check=True)
            (d / f"{name}.part.jf").rename(d / f"{name}.jf")
        print(f"  counted {name}", flush=True)
    counts = {}
    for name in ("reads", "genome"):
        dump = subprocess.run([str(JELLYFISH), "dump", "-c", "-t", str(d / f"{name}.jf")],
                              check=True, capture_output=True, text=True).stdout
        counts[name] = {m: int(c) for m, c in (line.split("\t") for line in dump.splitlines())}
    median = float(np.median([counts["reads"].get(m, 0) for m in every]))
    excess = {}
    for chrom, sites in site_kmers.items():
        excess[chrom] = {
            p: max(math.log2(max(counts["reads"].get(m, 0), 0.5)
                             / (median * max(1, counts["genome"].get(m, 0)))) for m in ms)
            for p, ms in sites.items() if ms}
    (d / "excess.json").write_text(json.dumps({"median": median, "excess": excess}))
    print(f"  median read count {median}", flush=True)


# Tuning, flags and scoring.

def cell_arrays(platform: str, depth: str, phaser: str, stratum, indel: set, kmer) -> dict:
    """Per pair of the cell: chromosome, positions, group, key, wrong, confident, features."""
    from phasebudget.confidence import wilson_lower

    out = {k: [] for k in ("chrom", "p1", "group", "key", "wrong", "confident", "f")}
    for chrom in CHROMS:
        t = json.loads(chrom_path(platform, depth, chrom).read_text())
        pos1, logd, bal = t["pos1"], t["logdepth"], t["balance"]
        km = (kmer or {}).get(chrom, {})
        for a, b, agree, disagree, wrong, mec in t["pairs"][phaser]:
            p, q = pos1[a], pos1[b]
            out["chrom"].append(chrom)
            out["p1"].append(p)
            out["group"].append("indel" if (chrom, p) in indel or (chrom, q) in indel
                                else str(stratum(chrom, p, q)))
            out["key"].append((wilson_lower(agree, agree + disagree), agree + disagree))
            out["wrong"].append(wrong)
            out["confident"].append(disagree == 0 and agree >= MIN_READS)
            ka = km.get(str(p), -math.inf) if kmer is not None else -math.inf
            kb = km.get(str(q), -math.inf) if kmer is not None else -math.inf
            out["f"].append((max(logd[a], logd[b]), max(bal[a], bal[b]), mec, max(ka, kb)))
    arr = {"chrom": out["chrom"], "p1": out["p1"], "key": out["key"],
           "group": np.array(out["group"]), "wrong": np.array(out["wrong"], bool),
           "confident": np.array(out["confident"], bool), "f": np.array(out["f"], float)}
    arr["odd"] = np.array([int(c[3:]) % 2 == 1 for c in arr["chrom"]])
    return arr


def subsets(features: tuple) -> list[tuple]:
    return [s for r in range(1, len(features) + 1) for s in combinations(features, r)]


def flag(f: np.ndarray, reference: np.ndarray, subset: tuple, tau: float) -> np.ndarray:
    """Pairs whose smallest tail rank over `subset`, against `reference` values, is <= tau."""
    ranks = [tail_rank(np.sort(reference[:, FEATURES.index(n)]), f[:, FEATURES.index(n)])
             for n in subset]
    return np.min(ranks, axis=0) <= tau


def tune(odd_f: np.ndarray, odd_wrong: np.ndarray, available: tuple) -> dict:
    """Per subset of `available`: `(tau, odd recall, odd flagged share)`; and the chosen one."""
    out = {}
    for s in subsets(available):
        ranks = [tail_rank(np.sort(odd_f[:, FEATURES.index(n)]), odd_f[:, FEATURES.index(n)])
                 for n in s]
        m = np.min(ranks, axis=0)
        tau = choose_tau(m, ABSTAIN_TUNE)
        flagged = m <= tau
        out[s] = (tau, flagged[odd_wrong].sum() / max(1, odd_wrong.sum()), flagged.mean())
    best = max(out, key=lambda s: (out[s][1], -len(s), [-available.index(n) for n in s]))
    return {"subsets": out, "chosen": best}


def evaluate(platform: str, depth: str, phaser: str, arr: dict, available: tuple,
             settings: dict | None = None) -> dict:
    """Tune (or take `settings`) per hard group, flag, and score baseline and flagged LTT."""
    import confidence_score as cs
    import stratified_score as ss

    from phasebudget.confidence import Units, coverage_grid, learn_then_test_by_group

    odd, f, group = arr["odd"], arr["f"], arr["group"]
    flagged = np.zeros(len(group), bool)
    groups, chosen = {}, {}
    for g in HARD:
        go, ge = (group == g) & odd, (group == g) & ~odd
        if not go.any():
            continue
        if settings is None:
            t = tune(f[go], arr["wrong"][go], available)
            subset = t["chosen"]
            tau = t["subsets"][subset][0]
        else:
            t = None
            subset, tau = settings[g]
        chosen[g] = (subset, tau)
        flagged[go] = flag(f[go], f[go], subset, tau)
        flagged[ge] = flag(f[ge], f[go], subset, tau)
        singles = {}
        for n in available:
            tau_n = choose_tau(np.min([tail_rank(np.sort(f[go][:, FEATURES.index(n)]),
                                                 f[go][:, FEATURES.index(n)])], axis=0),
                               ABSTAIN_TUNE)
            fe = flag(f[ge], f[go], (n,), tau_n)
            w = arr["wrong"][ge]
            singles[n] = (fe[w].sum() / max(1, w.sum()), fe.mean())
        keys = [arr["key"][i] for i in np.flatnonzero(go)]
        top = coverage_grid(keys, start=0.5, step=0.0025)[0]
        ekeys = [arr["key"][i] for i in np.flatnonzero(ge)]
        high = np.array([k >= top for k in ekeys], bool)
        w, c, fl = arr["wrong"][ge], arr["confident"][ge], flagged[ge]
        groups[g] = {
            "subset": subset, "tau": tau,
            "odd_recall": None if t is None else t["subsets"][subset][1],
            "odd_flagged": flagged[go].mean(),
            "pairs": int(ge.sum()), "wrong": int(w.sum()), "conf": int((w & c).sum()),
            "high": int((w & high).sum()),
            "flagged": int(fl.sum()), "fw": int((fl & w).sum()), "fc": int((fl & w & c).sum()),
            "fh": int((fl & w & high).sum()),
            "rate_flagged": (fl & w).sum() / max(1, fl.sum()),
            "rate_unflagged": (~fl & w).sum() / max(1, (~fl).sum()),
            "singles": singles,
        }
    rows = [(str(g), (c, p // cs.UNIT), k, bool(w))
            for g, c, p, k, w in zip(group, arr["chrom"], arr["p1"], arr["key"], arr["wrong"])]
    cal_all = [r for r, o in zip(rows, odd) if o]
    test_all = [r for r, o in zip(rows, odd) if not o]
    cal_kept = [r for r, o, x in zip(rows, odd, flagged) if o and not x]
    test_kept = [r for r, o, x in zip(rows, odd, flagged) if not o and not x]
    for g in groups:
        groups[g]["odd_units"] = len(Units.from_pairs([r[1:] for r in cal_all if r[0] == g]).best)
        groups[g]["odd_units_kept"] = len(
            Units.from_pairs([r[1:] for r in cal_kept if r[0] == g]).best)
    ltt = {}
    for beta in cs.BETAS:
        for name, cal, test in (("baseline", cal_all, test_all), ("flagged", cal_kept, test_kept)):
            t = learn_then_test_by_group(cal, beta, cs.DELTA)
            ltt[(name, beta)] = ss.cells(t, test, beta)
    return {"groups": groups, "ltt": ltt, "chosen": chosen}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", nargs=4, help=argparse.SUPPRESS)
    parser.add_argument("--kmers", action="store_true", help="run the k-mer stage only")
    parser.add_argument("--platforms", nargs="+", default=list(PLATFORMS))
    parser.add_argument("--rows-only", action="store_true", help="build the cached rows and stop")
    args = parser.parse_args()
    if args.build:
        platform, depth, chrom, path = args.build
        Path(path + ".part").write_text(json.dumps(build(chrom, depth)))
        Path(path + ".part").rename(path)
        return
    if args.kmers:
        kmers()
        return
    import confidence_score as cs
    import stratified_score as ss
    from stratified_gate import strata

    excess_path = OUT / "kmer" / "excess.json"
    kmer = json.loads(excess_path.read_text())["excess"] if excess_path.exists() else None
    print(f"Q2 (development, docs/design/q2-hard-flags.md): flags from the sample's own reads; "
          f"platforms {', '.join(f'{p} {PLATFORMS[p]}' for p in args.platforms)}; phasers "
          f"{PHASERS}; chromosomes {' '.join(CHROMS)}, odd tune and calibrate, even scored; hard "
          f"groups {HARD}, easy never flagged; features {FEATURES} (depth |log2| against the "
          f"chromosome median, allele balance -log10 binomial p, window MEC with {WINDOW} columns "
          f"each side, k-mer {K}-mer log2 copy excess "
          f"{'from data/r6_q2/kmer' if kmer else 'NOT AVAILABLE'} on {KMER_PLATFORMS} only); "
          f"tuning abstention {ABSTAIN_TUNE}; gate on {GATE}: recall >= {RECALL} (confident "
          f"wrong too where >= {MIN_CONFIDENT}; confident = >= {MIN_READS} reads, none "
          f"disagreeing), abstention <= {ABSTAIN_LIMIT} per cell; LTT units {cs.UNIT // 1000} kb, "
          f"delta {cs.DELTA}, betas {cs.BETAS}, grid from 0.5 step 0.0025; transfer "
          f"{SOURCE} settings to {TRANSFER}", flush=True)
    jobs = [(p, d, c) for p in args.platforms for d in PLATFORMS[p] for c in CHROMS]
    with ThreadPoolExecutor(WORKERS) as pool:
        list(pool.map(lambda job: chrom_path(*job), jobs))
    if args.rows_only:
        return
    stratum = strata(list(CHROMS))
    results = {}
    for platform in [p for p in (GATE, "w4_grch38", SOURCE, TRANSFER) if p in args.platforms]:
        indel = ss.indels(platform, CHROMS)
        use_kmer = platform in KMER_PLATFORMS and kmer is not None
        available = FEATURES if use_kmer else FEATURES[:3]
        for depth in PLATFORMS[platform]:
            for phaser in PHASERS:
                arr = cell_arrays(platform, depth, phaser, stratum, indel,
                                  kmer if use_kmer else None)
                settings = None
                if platform == TRANSFER and (SOURCE, depth, phaser) in results:
                    settings = results[(SOURCE, depth, phaser)]["chosen"]
                results[(platform, depth, phaser)] = evaluate(platform, depth, phaser, arr,
                                                              available, settings)
                report(platform, depth, phaser, results[(platform, depth, phaser)],
                       transferred=settings is not None)
    if GATE in args.platforms:
        gate(results)
    if TRANSFER in args.platforms:
        summary(results, TRANSFER)
    if SOURCE in args.platforms:
        summary(results, SOURCE)


def report(platform: str, depth: str, phaser: str, r: dict, transferred: bool) -> None:
    print(f"\n{platform} {depth}x {phaser}{' (settings from ' + SOURCE + ')' if transferred else ''}"
          ": per hard group, what the tuning optimised (subset, tau, odd recall at odd flagged "
          "share) | even: pairs, wrong, confident wrong, wrong in the top half; flagged share; "
          "recall of wrong, confident wrong, top-half wrong; error flagged / unflagged; odd units "
          "all / unflagged", flush=True)
    for g, x in r["groups"].items():
        odd = "transferred" if x["odd_recall"] is None else f"{x['odd_recall']:.3f}"
        print(f"  {g:<10} {'+'.join(x['subset']):<24} tau {x['tau']:.4f} odd recall {odd} at "
              f"{x['odd_flagged']:.3f} | {x['pairs']:>8} {x['wrong']:>6} {x['conf']:>5} "
              f"{x['high']:>5}; flagged {x['flagged'] / max(1, x['pairs']):.3f}; recall "
              f"{x['fw'] / max(1, x['wrong']):.3f} {x['fc'] / max(1, x['conf']):.3f} "
              f"{x['fh'] / max(1, x['high']):.3f}; error {x['rate_flagged']:.2e} / "
              f"{x['rate_unflagged']:.2e}; units {x['odd_units']} / {x['odd_units_kept']}")
        print("             single features (even recall, flagged share): " + ", ".join(
            f"{n} {a:.3f} {b:.3f}" for n, (a, b) in x["singles"].items()))
    print("  LTT per beta and group: baseline | flagged: certified, trusted share of all even "
          "group pairs, wrong units/units, p_exceed")
    for (name, beta), cells in sorted(r["ltt"].items(), key=lambda kv: (kv[0][1], kv[0][0])):
        if name != "baseline":
            continue
        flagged = r["ltt"][("flagged", beta)]
        for g in HARD + ("easy",):
            if g not in cells:
                continue
            b, f = cells[g], flagged.get(g)
            total = b[3]
            line = (f"  beta {beta:<7g} {g:<10} {'yes' if b[7] else 'no ':<3} "
                    f"{b[4] / max(1, total):.4f} {b[1]}/{b[0]} {b[2]:.2e} | ")
            line += ("-" if f is None else
                     f"{'yes' if f[7] else 'no ':<3} {f[4] / max(1, total):.4f} {f[1]}/{f[0]} "
                     f"{f[2]:.2e}")
            print(line)


def gate(results: dict) -> None:
    cells = [(d, ph) for d in PLATFORMS[GATE] for ph in PHASERS]
    level = 0.05 / (len(HARD) * len(cells) * 3)
    print(f"\nGate on {GATE}, even chromosomes, {len(cells)} cells; violation level {level:.2e}")
    violating = []
    for d, ph in cells:
        for (name, beta), c in results[(GATE, d, ph)]["ltt"].items():
            for g in HARD:
                if name == "flagged" and g in c and c[g][7] and c[g][0] and c[g][2] < level:
                    violating.append(f"{d}x {ph} beta {beta:g} {g} p {c[g][2]:.2e}")
    base_bad = sum(1 for d, ph in cells for (name, beta), c in results[(GATE, d, ph)]["ltt"].items()
                   for g in HARD if name == "baseline" and g in c and c[g][7] and c[g][0]
                   and c[g][2] < level)
    b_ok = not violating
    print(f"  (b) flagged violating cells {len(violating)} of {len(HARD) * len(cells) * 3} "
          f"(baseline {base_bad})" + "".join(f"\n      {v}" for v in violating))
    passing = []
    for g in HARD:
        xs = [results[(GATE, d, ph)]["groups"][g] for d, ph in cells
              if g in results[(GATE, d, ph)]["groups"]]
        wrong, fw = sum(x["wrong"] for x in xs), sum(x["fw"] for x in xs)
        conf, fc = sum(x["conf"] for x in xs), sum(x["fc"] for x in xs)
        recalls = [x["fw"] / max(1, x["wrong"]) for x in xs]
        shares = [x["flagged"] / max(1, x["pairs"]) for x in xs]
        recall, crecall = fw / max(1, wrong), fc / max(1, conf)
        a_ok = (recall >= RECALL and (conf < MIN_CONFIDENT or crecall >= RECALL)
                and max(shares) <= ABSTAIN_LIMIT)
        before = [(d, ph) for d, ph in cells
                  if not results[(GATE, d, ph)]["ltt"][("baseline", 1e-3)][g][7]]
        after = [(d, ph) for d, ph in before
                 if results[(GATE, d, ph)]["ltt"][("flagged", 1e-3)].get(g, (0,) * 8)[7]]
        c_ok = bool(before) and len(after) >= len(before) / 2
        print(f"  {g:<10} (a) recall {recall:.3f} (cells {min(recalls):.3f} to {max(recalls):.3f}), "
              f"confident {crecall:.3f} of {conf}{' (not judged)' if conf < MIN_CONFIDENT else ''}, "
              f"abstention mean {np.mean(shares):.3f} worst {max(shares):.3f} -> "
              f"{'ok' if a_ok else 'FAIL'}; (c) baseline certified nothing at 1e-3 in "
              f"{len(before)} cells, flagged certifies {len(after)} of them -> "
              f"{'ok' if c_ok else 'FAIL'}")
        if a_ok and c_ok:
            passing.append(g)
    print(f"  (b) {'ok' if b_ok else 'FAIL'}; groups meeting (a) and (c): "
          f"{', '.join(passing) or 'none'}")
    print(f"Q2 gate: {'passed' if b_ok and passing else 'not passed'}")


def summary(results: dict, platform: str) -> None:
    """Pooled recall and abstention per hard group, and certification at 1e-3, for a platform."""
    cells = [(d, ph) for d in PLATFORMS[platform] for ph in PHASERS]
    level = 0.05 / (len(HARD) * len(cells) * 3)
    print(f"\n{platform} summary (reported, not gated), level {level:.2e}")
    for g in HARD:
        xs = [results[(platform, d, ph)]["groups"].get(g) for d, ph in cells]
        xs = [x for x in xs if x]
        wrong, fw = sum(x["wrong"] for x in xs), sum(x["fw"] for x in xs)
        conf, fc = sum(x["conf"] for x in xs), sum(x["fc"] for x in xs)
        shares = [x["flagged"] / max(1, x["pairs"]) for x in xs]
        bad = sum(1 for d, ph in cells for (name, beta), c in
                  results[(platform, d, ph)]["ltt"].items()
                  if name == "flagged" and g in c and c[g][7] and c[g][0] and c[g][2] < level)
        before = [(d, ph) for d, ph in cells
                  if not results[(platform, d, ph)]["ltt"][("baseline", 1e-3)][g][7]]
        after = [k for k in before
                 if results[(platform, *k)]["ltt"][("flagged", 1e-3)].get(g, (0,) * 8)[7]]
        print(f"  {g:<10} recall {fw / max(1, wrong):.3f} of {wrong}, confident "
              f"{fc / max(1, conf):.3f} of {conf}, abstention mean {np.mean(shares):.3f} worst "
              f"{max(shares):.3f}; flagged violating cells {bad}; uncertified at 1e-3 before "
              f"{len(before)}, certified after {len(after)}")


if __name__ == "__main__":
    main()
