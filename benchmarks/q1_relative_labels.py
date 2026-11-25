"""Round 6, Q1 (development): siblings without parents as calibration labels.

`docs/design/q1-relative-labels.md` states the claim and the gate. For a
child A and a sibling B, a three-state HMM (IBD0, IBD1, IBD2) over the two
genotypes alone finds segments where they share one haplotype. Inside a
confident IBD1 run, a site where A is heterozygous and B homozygous puts B's
allele on the shared haplotype of A; two such sites in one run give A's
relation between them. Children are the five Platinum Pedigree ONT runs
(`ped_*`); sibling genotypes come from the Pedigree truthset with phase
removed, exactly (`exact`) or with random genotype errors (`e1e-3`, `e1e-2`).
Odd chromosomes calibrate on eligible pairs against the labels; even
chromosomes are scored against the Pedigree truth on eligible and all pairs.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/q1_relative_labels.py`
"""

from __future__ import annotations

import bisect
import math
import random
import sys
from itertools import combinations
from pathlib import Path
from typing import Optional, Sequence

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

CHILDREN = ("NA12879", "NA12881", "NA12882", "NA12885", "NA12886")
GROUPS = ("easy", "difficult", "lowmap", "segdup", "tandem", "indel")
LAMBDA = 2e-8  # flips per bp of each parent's IBD indicator: two meioses at about 1 cM/Mb
TAU = 0.99  # posterior P(IBD1) every site of a run must reach
MARGIN = 100_000  # bp trimmed from each end of a run
EM_ITERATIONS = 4
VERSIONS = (("exact", 0.0), ("e1e-3", 1e-3), ("e1e-2", 1e-2))
SEED = 17
GATE_LABEL_ERROR = 1e-4
WINDOW = 250_000  # amendment 1: bp each side of a labelled site
MIN_IBS1 = 0.25  # amendment 1: least share of IBS1 sites in the window
LABELLER = "both"  # amendment 1; "hmm" reproduces the first run


def category(ga: int, gb: int) -> Optional[int]:
    """Observation from two alt-allele counts: 0 opposite homozygotes (IBS0), 1 one heterozygous
    and one homozygous (IBS1), 2 both heterozygous, 3 both homozygous alt; None for both 0/0."""
    if ga < 0 or gb < 0 or ga == gb == 0:
        return None
    if ga == gb:
        return 2 if ga == 1 else 3
    return 1 if 1 in (ga, gb) else 0


def transition(d: float, lam: float = LAMBDA) -> tuple:
    """Rows from IBD0, IBD1, IBD2 over `d` bp: two independent indicators that flip at rate `lam`."""
    s = 0.5 * (1 + math.exp(-2 * lam * d))
    f = 1 - s
    return ((s * s, 2 * s * f, f * f), (s * f, s * s + f * f, s * f), (f * f, 2 * s * f, s * s))


INIT = (0.25, 0.5, 0.25)
EMIT0 = ((0.05, 0.45, 0.25, 0.25), (1e-3, 0.5, 0.25, 0.249), (1e-3, 1e-3, 0.5, 0.498))


def forward_backward(pos: Sequence[int], obs: Sequence[int], emit, lam: float = LAMBDA
                     ) -> tuple[list, list, list]:
    """Posterior state probabilities per site, scaled forward and backward passes."""
    n = len(pos)
    e0, e1, e2 = emit
    f0, f1, f2 = [0.0] * n, [0.0] * n, [0.0] * n
    x0, x1, x2 = INIT[0] * e0[obs[0]], INIT[1] * e1[obs[0]], INIT[2] * e2[obs[0]]
    z = x0 + x1 + x2
    f0[0], f1[0], f2[0] = x0 / z, x1 / z, x2 / z
    exp = math.exp
    for t in range(1, n):
        s = 0.5 * (1 + exp(-2 * lam * (pos[t] - pos[t - 1])))
        f = 1 - s
        ss, ff, sf = s * s, f * f, s * f
        x0, x1, x2 = f0[t - 1], f1[t - 1], f2[t - 1]
        o = obs[t]
        y0 = (x0 * ss + x1 * sf + x2 * ff) * e0[o]
        y1 = ((x0 + x2) * 2 * sf + x1 * (ss + ff)) * e1[o]
        y2 = (x0 * ff + x1 * sf + x2 * ss) * e2[o]
        z = y0 + y1 + y2
        f0[t], f1[t], f2[t] = y0 / z, y1 / z, y2 / z
    b0 = b1 = b2 = 1.0
    for t in range(n - 1, -1, -1):
        g0, g1, g2 = f0[t] * b0, f1[t] * b1, f2[t] * b2
        z = g0 + g1 + g2
        f0[t], f1[t], f2[t] = g0 / z, g1 / z, g2 / z  # f now holds the posterior
        if t == 0:
            break
        s = 0.5 * (1 + exp(-2 * lam * (pos[t] - pos[t - 1])))
        f = 1 - s
        ss, ff, sf = s * s, f * f, s * f
        o = obs[t]
        w0, w1, w2 = e0[o] * b0, e1[o] * b1, e2[o] * b2
        b0, b1, b2 = ss * w0 + 2 * sf * w1 + ff * w2, sf * w0 + (ss + ff) * w1 + sf * w2, \
            ff * w0 + 2 * sf * w1 + ss * w2
        z = b0 + b1 + b2
        b0, b1, b2 = b0 / z, b1 / z, b2 / z
    return f0, f1, f2


def ibd_posterior(chunks: Sequence[tuple[Sequence[int], Sequence[int]]],
                  iterations: int = EM_ITERATIONS) -> tuple[list, tuple]:
    """P(IBD1) per site of each `(positions, observations)` chunk (a chromosome), with one set
    of emissions fitted by EM over all chunks (transitions fixed); and the emissions. Fitting
    per chromosome fails when a chromosome lies in one state."""
    emit = EMIT0
    for it in range(iterations + 1):
        posts = [forward_backward(pos, obs, emit) for pos, obs in chunks]
        if it == iterations:
            break
        counts = [[1e-9] * 4 for _ in range(3)]
        for post, (_, obs) in zip(posts, chunks):
            for k in range(3):
                row, p = counts[k], post[k]
                for t, o in enumerate(obs):
                    row[o] += p[t]
        emit = tuple(tuple(c / sum(row) for c in row) for row in counts)
    return [post[1] for post in posts], emit


def site_labels(pos: Sequence[int], ga: Sequence[int], gb: Sequence[int], p1: Sequence[float],
                tau: float = TAU, margin: int = MARGIN) -> dict:
    """`{pos: (allele on the shared haplotype, run)}` at sites where A is heterozygous, B is
    homozygous, inside a run of consecutive sites with P(IBD1) >= `tau`, at least `margin` bp
    from both ends of the run."""
    out: dict = {}
    n, t, run = len(pos), 0, 0
    while t < n:
        if p1[t] < tau:
            t += 1
            continue
        end = t
        while end + 1 < n and p1[end + 1] >= tau:
            end += 1
        lo, hi = pos[t] + margin, pos[end] - margin
        for k in range(t, end + 1):
            if ga[k] == 1 and gb[k] in (0, 2) and lo <= pos[k] <= hi:
                out[pos[k]] = (gb[k] // 2, run)
        run += 1
        t = end + 1
    return out


def sibling_relation(labels: dict, a: int, b: int) -> Optional[bool]:
    """A's relation (same first allele) between sites `a` and `b` from one sibling's labels."""
    x, y = labels.get(a), labels.get(b)
    if x is None or y is None or x[1] != y[1]:
        return None
    return x[0] == y[0]


def combined(relations: Sequence[Optional[bool]]) -> Optional[bool]:
    """The relation when every sibling that gives one agrees; None otherwise."""
    given = {r for r in relations if r is not None}
    return given.pop() if len(given) == 1 else None


def window_labels(pos: Sequence[int], ga: Sequence[int], gb: Sequence[int], w: int = WINDOW,
                  min_ibs1: float = MIN_IBS1) -> dict:
    """Amendment 1: `{pos: (allele on the shared haplotype, segment)}` at sites where A is
    heterozygous and B homozygous, whose sites within `w` bp on both sides hold no IBS0 site and
    at least `min_ibs1` IBS1 sites; the segment counts the IBS0 sites before the site, so two
    labels relate only with no IBS0 site between them. Sites must have a category."""
    cats = [category(x, y) for x, y in zip(ga, gb)]
    ibs0, ibs1 = [0], [0]
    for c in cats:
        ibs0.append(ibs0[-1] + (c == 0))
        ibs1.append(ibs1[-1] + (c == 1))
    out: dict = {}
    for k, p in enumerate(pos):
        if ga[k] != 1 or gb[k] not in (0, 2):
            continue
        lo, hi = bisect.bisect_left(pos, p - w), bisect.bisect_right(pos, p + w)
        if ibs0[hi] - ibs0[lo] == 0 and ibs1[hi] - ibs1[lo] >= min_ibs1 * (hi - lo):
            out[p] = (gb[k] // 2, ibs0[k])
    return out


def pair_labels(chroms: Sequence[tuple[Sequence[int], Sequence[int], Sequence[int]]],
                tau: float = TAU, margin: int = MARGIN, labeller: Optional[str] = None
                ) -> tuple[list, list, Optional[tuple]]:
    """Per chromosome `(positions, A's alt counts, B's)`: labels of A from B and of B from A,
    and the fitted emissions (None for the window rule). Sites where both are 0/0 or either is
    missing are skipped."""
    kept = []
    for pos, ga, gb in chroms:
        keep = [k for k in range(len(pos)) if category(ga[k], gb[k]) is not None]
        kept.append(([pos[k] for k in keep], [ga[k] for k in keep], [gb[k] for k in keep]))
    kind = labeller or LABELLER
    if kind == "window":
        return ([window_labels(p, a, b) for p, a, b in kept],
                [window_labels(p, b, a) for p, a, b in kept], None)
    p1s, emit = ibd_posterior([(p, [category(x, y) for x, y in zip(a, b)])
                               for p, a, b in kept if p])
    la, lb, it = [], [], iter(p1s)
    for p, a, b in kept:
        p1 = next(it) if p else []
        la.append(site_labels(p, a, b, p1, tau, margin))
        lb.append(site_labels(p, b, a, p1, tau, margin))
    if kind == "both":  # amendment 1: both rules must label the site, in one run and segment
        la = [both(h, window_labels(p, a, b)) for h, (p, a, b) in zip(la, kept)]
        lb = [both(h, window_labels(p, b, a)) for h, (p, a, b) in zip(lb, kept)]
    return la, lb, emit


def both(hmm: dict, window: dict) -> dict:
    """Sites both labellers label, with the pair (HMM run, window segment) as the segment."""
    return {p: (x[0], (x[1], window[p][1])) for p, x in hmm.items() if p in window}


def perturb(g: int, rate: float, rng: random.Random) -> int:
    """With probability `rate`, one of the other two genotypes at random."""
    if g < 0 or rng.random() >= rate:
        return g
    return rng.choice([x for x in (0, 1, 2) if x != g])


def simulate_family(n_sites: int, length: int, children: int, seed: int, rate: float = 1e-8
                    ) -> tuple[list, list, list]:
    """Two unrelated parents, `children` children with Poisson crossovers at `rate` per bp per
    meiosis. Returns positions, per child alt counts, per child (paternal, maternal) alleles."""
    rng = random.Random(seed)
    pos = sorted(rng.sample(range(length), n_sites))
    freq = [min(0.95, 10 ** rng.uniform(-2.3, 0)) for _ in pos]
    parents = [[[int(rng.random() < f) for f in freq] for _ in range(2)] for _ in range(2)]

    def gamete(parent):
        cuts = sorted(rng.uniform(0, length) for _ in range(_poisson(rng, rate * length)))
        src, out, c = rng.randrange(2), [], 0
        for k, p in enumerate(pos):
            while c < len(cuts) and cuts[c] <= p:
                src, c = 1 - src, c + 1
            out.append(parent[src][k])
        return out

    genotypes, haplotypes = [], []
    for _ in range(children):
        pat, mat = gamete(parents[0]), gamete(parents[1])
        genotypes.append([x + y for x, y in zip(pat, mat)])
        haplotypes.append((pat, mat))
    return pos, genotypes, haplotypes


def _poisson(rng: random.Random, mean: float) -> int:
    k, t = 0, rng.expovariate(1.0)
    while t < mean:
        k, t = k + 1, t + rng.expovariate(1.0)
    return k


def simulate_check(sites: int, length: int, seed: int, tau: float = TAU, margin: int = MARGIN,
                   rate: float = 0.0, labeller: Optional[str] = None) -> tuple[int, int, int]:
    """`(A's consecutive het pairs, labelled, wrong)` for child 0 labelled by child 1."""
    pos, g, h = simulate_family(sites, length, 2, seed)
    rng = random.Random(seed + 1)
    obs = [[perturb(x, rate, rng) for x in gi] for gi in g]
    (labels,), _, _ = pair_labels([(pos, obs[0], obs[1])], tau, margin, labeller)
    het = [k for k in range(len(pos)) if g[0][k] == 1]
    labelled = wrong = 0
    for i, j in zip(het, het[1:]):
        r = sibling_relation(labels, pos[i], pos[j])
        if r is not None:
            labelled += 1
            wrong += r != (h[0][0][i] == h[0][0][j])
    return len(het) - 1, labelled, wrong


# The run on the server.

PARAMS: dict = {}
TRUTH: dict = {}  # chrom -> (positions, {sample: alt counts}, {child: first allele or -1})
LABELS: dict = {}  # (version, child, sibling, chrom) -> site labels


def read_truthset(chroms) -> list:
    """Biallelic truthset records inside the high-quality regions, one per position."""
    import gzip

    import genome_common as common
    import pysam

    beds: dict = {c: ([], []) for c in chroms}
    with gzip.open(common.PEDIGREE / "hq_regions_final.bed.gz", "rt") as handle:
        for line in handle:
            f = line.split("\t")
            if f[0] in beds:
                beds[f[0]][0].append(int(f[1]))
                beds[f[0]][1].append(int(f[2]))
    with pysam.TabixFile(str(common.PEDIGREE / "CEPH1463.GRCh38.family-truthset.ov.vcf.gz")) as tb:
        samples = tb.header[-1].split("\t")[9:]
        for chrom in chroms:
            starts, ends = beds[chrom]
            seen: dict = {}
            for line in tb.fetch(chrom):
                f = line.rstrip("\n").split("\t")
                pos = int(f[1])
                k = bisect.bisect_right(starts, pos - 1) - 1
                if "," in f[4] or k < 0 or pos - 1 >= ends[k]:
                    continue
                gts = [x.split(":")[0] for x in f[9:]]
                seen[pos] = None if pos in seen else gts
            pos_list = sorted(p for p, g in seen.items() if g is not None)
            counts = {s: [] for s in samples}
            first = {c: [] for c in CHILDREN}
            for p in pos_list:
                for s, gt in zip(samples, seen[p]):
                    alleles = gt.replace("|", "/").split("/")
                    ok = len(alleles) == 2 and set(alleles) <= {"0", "1"}
                    counts[s].append(sum(map(int, alleles)) if ok else -1)
                    if s in first:
                        first[s].append(int(gt[0]) if ok and gt in ("0|1", "1|0") else -1)
            TRUTH[chrom] = (pos_list, counts, first)
    return samples


def _labels_task(task):
    version, rate, a, b, chroms = task

    def genotypes(sample, chrom):
        return _perturbed(version, rate, sample, chrom) if rate else TRUTH[chrom][1][sample]

    la, lb, emit = pair_labels([(TRUTH[c][0], genotypes(a, c), genotypes(b, c)) for c in chroms])
    return task, la, lb, emit


def _perturbed(version, rate, sample, chrom):
    """A sample's genotypes with errors, the same whichever pair asks for them."""
    rng = random.Random(f"{SEED}-{version}-{sample}-{chrom}")
    return [perturb(x, rate, rng) for x in TRUTH[chrom][1][sample]]


def _strata():
    from stratified_gate import strata

    import genome_common as common
    return strata(list(common.CONFIRMATORY))


STRATUM = None


def _cell_task(task):
    """Everything for one (child, depth, phaser): per version and sibling set, label error,
    eligibility and the Learn then Test cells."""
    import confidence_score as cs
    import genome_common as common
    import stratified_score as ss

    from phasebudget.confidence import learn_then_test_by_group

    child, depth, phaser = task
    chroms = list(common.CONFIRMATORY)
    odd = {c for c in chroms if int(c[3:]) % 2}
    platform = f"ped_{child}"
    ind = ss.indels(platform, chroms)
    sibs = [c for c in CHILDREN if c != child]
    index = {c: {p: k for k, p in enumerate(TRUTH[c][0])} for c in chroms}
    rows = []
    for c, a, key, wrong, b, _ in cs.load_rows(platform, depth, phaser, chroms):
        g = "indel" if (c, a) in ind or (c, b) in ind else STRATUM(c, a, b)
        first = TRUTH[c][2][child]
        ia, ib = index[c].get(a), index[c].get(b)
        t = None if ia is None or ib is None or first[ia] < 0 or first[ib] < 0 \
            else first[ia] == first[ib]
        rows.append((g, (c, a // cs.UNIT), key, wrong, c, a, b, t))
    out = {"pairs": len(rows), "groups": {}, "coverage": {}, "label": {}, "cells": []}
    for g in GROUPS:
        out["groups"][g] = sum(r[0] == g for r in rows)
    sets = {"all": tuple(sibs), "two": tuple(sibs[:2]), "one": tuple(sibs[:1])}
    for version, _ in VERSIONS:
        rel = [tuple(sibling_relation(LABELS[(version, child, s, r[4])], r[5], r[6])
                     for s in sibs) for r in rows]
        if version == "exact":
            for size in (1, 2, len(sibs)):
                shares = []
                for chosen in combinations(range(len(sibs)), size):
                    n = sum(r[7] is not None and combined([x[i] for i in chosen]) is not None
                            for r, x in zip(rows, rel))
                    shares.append(n / max(1, len(rows)))
                out["coverage"][size] = (sum(shares) / len(shares), min(shares), max(shares))
        for name, chosen in sets.items():
            if version != "exact" and name != "all":
                continue
            idx = [sibs.index(s) for s in chosen]
            lab = [None if r[7] is None else combined([x[i] for i in idx])
                   for r, x in zip(rows, rel)]
            conflicts = sum(r[7] is not None and len({x[i] for i in idx} - {None}) > 1
                            for r, x in zip(rows, rel))
            per = {}
            for g in GROUPS:
                mine = [(r, lb) for r, lb in zip(rows, lab) if r[0] == g and lb is not None]
                per[g] = (len(mine), sum(lb != r[7] for r, lb in mine))
            out["label"][(version, name)] = (per, conflicts)
            cal = [(r[0], r[1], r[2], r[3] != (lb != r[7]))
                   for r, lb in zip(rows, lab) if lb is not None and r[4] in odd]
            test_eligible = [r[:4] for r, lb in zip(rows, lab) if lb is not None and r[4] not in odd]
            test_all = [r[:4] for r in rows if r[4] not in odd]
            for beta in cs.BETAS:
                t = learn_then_test_by_group(cal, beta, cs.DELTA)
                for scope, test in (("eligible", test_eligible), ("all", test_all)):
                    for g, cell in ss.cells(t, test, beta).items():
                        out["cells"].append((version, name, scope, child, depth, phaser, beta, g,
                                             cell))
    return task, out


def main() -> None:
    import argparse
    import time
    from multiprocessing import Pool

    import confidence_score as cs
    import genome_common as common

    global STRATUM, LABELLER
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--labeller", choices=("both", "window", "hmm"), default=LABELLER)
    args = parser.parse_args()
    LABELLER = args.labeller
    chroms = list(common.CONFIRMATORY)
    print(f"Q1 (development, docs/design/q1-relative-labels.md): children {CHILDREN} (ONT ped_*), "
          "each labelled by the other four as siblings from the Platinum Pedigree truthset "
          "(CEPH1463.GRCh38.family-truthset.ov.vcf.gz, phase removed; biallelic records inside "
          "hq_regions_final.bed, one per position); HMM states IBD0/1/2, two IBD indicators "
          f"flipping at {LAMBDA:g} per bp, initial {INIT}, emissions over (IBS0, IBS1, both het, "
          f"both hom alt) fitted per pair over all chromosomes jointly by EM ({EM_ITERATIONS} "
          "iterations, transitions fixed, start "
          f"{EMIT0}); runs with P(IBD1) >= {TAU}, {MARGIN} bp trimmed from each run end; sibling "
          "relations combined when all agree; genotype versions "
          f"{VERSIONS} (seed {SEED}, each genotype replaced by another with that probability); "
          "sibling sets all, two (first two of the others in "
          f"{CHILDREN} order), one (first); calibrate on eligible odd-chromosome pairs against "
          "the labels, score even-chromosome pairs against the Pedigree truth, eligible only and "
          f"all; chromosomes {len(chroms)} confirmatory; {cs.UNIT // 1000} kb units, delta "
          f"{cs.DELTA}; betas {cs.BETAS}; depths {cs.DEPTHS}; phasers {cs.PHASERS}; groups {GROUPS}; "
          f"gate label error <= {GATE_LABEL_ERROR:g}; workers {args.workers}; labeller "
          f"{LABELLER}" + (f" (amendment 1: the HMM run rule above and, at the same site, no "
                           f"IBS0 site and at least {MIN_IBS1:g} IBS1 sites within {WINDOW} bp "
                           f"each side)" if LABELLER == "both" else ""), flush=True)
    start = time.time()
    read_truthset(chroms)
    print(f"truthset read: {sum(len(v[0]) for v in TRUTH.values())} records "
          f"({time.time() - start:.0f} s)", flush=True)
    tasks = [(v, r, a, b, tuple(chroms)) for v, r in VERSIONS
             for a, b in combinations(CHILDREN, 2)]
    print("fitted emissions per version and pair "
          "(rows IBD0, IBD1, IBD2; columns IBS0, IBS1, both het, both hom alt):", flush=True)
    with Pool(args.workers) as pool:
        for (v, _, a, b, _), la, lb, emit in pool.imap_unordered(_labels_task, tasks):
            for c, x, y in zip(chroms, la, lb):
                LABELS[(v, a, b, c)], LABELS[(v, b, a, c)] = x, y
            if emit is not None:
                print(f"  {v} {a}-{b}: " + "; ".join(" ".join(f"{x:.2e}" for x in row)
                                                    for row in emit), flush=True)
    print(f"labels done ({time.time() - start:.0f} s)", flush=True)
    STRATUM = _strata()
    cell_tasks = [(c, d, p) for c in CHILDREN for d in cs.DEPTHS for p in cs.PHASERS]
    with Pool(min(args.workers, len(cell_tasks))) as pool:
        results = dict(pool.imap_unordered(_cell_task, cell_tasks))
    print(f"cells done ({time.time() - start:.0f} s)", flush=True)
    report(results, cs)


def report(results: dict, cs) -> None:
    print("\nlabelled-pair coverage (share of scored pairs with a label, 30x whatshap), mean "
          "[min, max] over sibling choices, exact genotypes")
    for child in CHILDREN:
        cov = results[(child, "30", "whatshap")]["coverage"]
        print(f"  {child}: " + "; ".join(f"{k} sibling{'s' if k > 1 else ''} {m:.3f} "
                                         f"[{lo:.3f}, {hi:.3f}]" for k, (m, lo, hi) in cov.items()))
    keys = [(v, "all") for v, _ in VERSIONS] + [("exact", "two"), ("exact", "one")]
    print("\nlabel error per pair (wrong labels / labelled pairs) and eligible share, per "
          "version and sibling set, per depth and phaser pooled over children; conflicts are "
          "pairs where siblings disagree")
    worst = {}
    for v, name in keys:
        for depth in cs.DEPTHS:
            for phaser in cs.PHASERS:
                n = k = pairs = conf = 0
                for child in CHILDREN:
                    r = results[(child, depth, phaser)]
                    per, c = r["label"][(v, name)]
                    n += sum(x[0] for x in per.values())
                    k += sum(x[1] for x in per.values())
                    pairs += r["pairs"]
                    conf += c
                worst.setdefault((v, name), []).append(k / max(1, n))
                print(f"  {v:<6} {name:<4} {depth}x {phaser:<10} {k}/{n} = {k / max(1, n):.2e}; "
                      f"eligible {n / max(1, pairs):.3f}; conflicts {conf}")
    print("\nper child, exact genotypes, all siblings, 30x whatshap: per group eligible share, "
          "label errors/labelled")
    for child in CHILDREN:
        r = results[(child, "30", "whatshap")]
        per, _ = r["label"][("exact", "all")]
        print(f"  {child}: " + "; ".join(f"{g} {per[g][0] / max(1, r['groups'][g]):.3f} "
                                         f"{per[g][1]}/{per[g][0]}" for g in GROUPS))
    cells = [c for r in results.values() for c in r["cells"]]
    for v, name in keys:
        for scope in ("eligible", "all"):
            mine_cells = [c for c in cells if c[:3] == (v, name, scope)]
            level = 0.05 / max(1, len(mine_cells))
            print(f"\n{v} genotypes, {name} siblings, scored on {scope} even-chromosome pairs; "
                  f"violation level {level:.2e} over {len(mine_cells)} cells; per beta and group, "
                  "pooled over children, depths and phasers: cells certifying, violating, min p, "
                  "trusted share of the scored pairs, wrong units/units")
            for beta in cs.BETAS:
                for g in GROUPS:
                    mine = [c[8] for c in mine_cells if c[6] == beta and c[7] == g]
                    if not mine:
                        continue
                    cert = [c for c in mine if c[7]]
                    bad = [c for c in cert if c[0] and c[2] < level]
                    n, k = sum(c[0] for c in mine), sum(c[1] for c in mine)
                    share = sum(c[4] for c in mine) / max(1, sum(c[3] for c in mine))
                    minp = min((c[2] for c in cert), default=1.0)
                    print(f"  {beta:<7g} {g:<10} {len(cert):>2} of {len(mine):<2} {len(bad):>3}  "
                          f"{minp:.2e}  {share:.4f}  {k}/{n} = {k / max(1, n):.2e}")
            if scope == "eligible":
                for c in mine_cells:
                    if c[8][7] and c[8][0] and c[8][2] < level:
                        print(f"  violating: {c[3]} {c[4]}x {c[5]} beta {c[6]:g} {c[7]} "
                              f"p {c[8][2]:.2e} ({c[8][1]}/{c[8][0]})")
        easy = {(c[2], c[3], c[4], c[5]): c[8] for c in cells
                if c[:2] == (v, name) and c[6] == 1e-3 and c[7] == "easy"}
        print(f"\n{v} {name}: easy at beta 1e-3 per child, depth, phaser: certified; trusted "
              "share of eligible easy pairs; trusted eligible over all easy pairs")
        for child in CHILDREN:
            parts = []
            for depth in cs.DEPTHS:
                for phaser in cs.PHASERS:
                    e, a = easy.get(("eligible", child, depth, phaser)), \
                        easy.get(("all", child, depth, phaser))
                    if e is None or a is None:
                        parts.append(f"{depth}x {phaser[:2]} -")
                        continue
                    parts.append(f"{depth}x {phaser[:2]} {'yes' if e[7] else 'no'} "
                                 f"{e[4] / max(1, e[3]):.3f} {e[4] / max(1, a[3]):.3f}")
            print(f"  {child}: " + "; ".join(parts))
    # Gate: exact genotypes, all siblings.
    label_max = max(worst[("exact", "all")])
    gate_cells = [c for c in cells if c[:3] == ("exact", "all", "eligible")]
    level = 0.05 / max(1, len(gate_cells))
    bad = [c for c in gate_cells if c[8][7] and c[8][0] and c[8][2] < level]
    easy = [c for c in cells if c[:2] == ("exact", "all") and c[6] == 1e-3 and c[7] == "easy"]
    cert = {(c[3], c[4], c[5]) for c in easy if c[2] == "eligible" and c[8][7]}
    trusted = sum(c[8][4] for c in easy if c[2] == "eligible")
    total = sum(c[8][3] for c in easy if c[2] == "all")
    eligible_share = sum(sum(x[0] for x in r["label"][("exact", "all")][0].values()) / r["pairs"]
                         for r in results.values()) / len(results)
    a_ok = label_max <= GATE_LABEL_ERROR
    b_ok = not bad
    need = len(CHILDREN) * len(cs.DEPTHS) * len(cs.PHASERS)
    c_ok = len(cert) == need
    print(f"\nQ1 gate (exact genotypes, all siblings): (a) worst label error over depth and "
          f"phaser {label_max:.2e} {'ok' if a_ok else 'FAIL'}; (b) {len(bad)} violating of "
          f"{len(gate_cells)} cells at level {level:.2e} {'ok' if b_ok else 'FAIL'}; (c) easy at "
          f"1e-3 certifies in {len(cert)} of {need} cells {'ok' if c_ok else 'FAIL'}, trusted "
          f"share of all even-chromosome easy pairs {trusted / max(1, total):.4f}; mean eligible "
          f"share of pairs {eligible_share:.3f}")
    print(f"Q1 gate: {'passed' if a_ok and b_ok and c_ok else 'not passed'}")


if __name__ == "__main__":
    main()
