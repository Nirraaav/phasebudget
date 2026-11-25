"""L2 (`docs/design/l2-end-variant-calls.md`): within-gene cis or trans calls trusted when the
person's assembly (primary arm) or the reads spanning both ends (descriptive arm) give the
phaser's relation at the two end variants; per-call error calibrated on K4b's six-genome pool.

Per (platform, chromosome) a subprocess in the platform's environment writes, for truth-phased
input variants inside a listed gene, `[pos1, hard, truth, whatshap, longphase, assembly]`
(phaser and assembly values from K1b's tables) and the read support `[i, j, same, diff]` of
pairs at most 100 kb apart, to `data/l2/`. Then each (genome, chromosome) is reduced to
per-group tallies and per-(genome, gene) CRC picks, in a process pool.

`--dry-run` reads no test genome: HG002's odd chromosomes in two pseudo-genomes plus the five
Pedigree genomes form the pool, HG002's even chromosomes the test.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/l2_end_calls.py [--dry-run]`
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import json
import os
import subprocess
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

ROOT = HERE.parent
DEPTH = "30"
CHROMS = tuple(f"chr{i}" for i in range(1, 23) if i != 20)
PHASERS = ("whatshap", "longphase")
ARMS = ("assembly", "reads")
TARGETS = (1e-3, 1e-2)
DELTA = 0.05
WINDOW = 100_000
DISTANCES = (1_000, 10_000, WINDOW)
CLASSES = ("easy", "hard")
GROUPS = tuple((d, c) for d in DISTANCES for c in CLASSES)
POOL = {"HG002": "hprc_HG002", **{s: f"phifi_{s}" for s in
                                  ("NA12879", "NA12881", "NA12882", "NA12885", "NA12886")}}
TESTS = {"NA12877": "phifi_NA12877", "HG005": "k1_hg005"}
CELLS = len(TESTS) * len(PHASERS) * len(TARGETS)  # 8, as K4b
LEVEL = 0.05 / CELLS
GOAL_R, GOAL_SHARE, GOAL_BIN = 1e-2, 0.70, WINDOW
K4B_LONG = {("NA12877", "whatshap"): 0.5529, ("NA12877", "longphase"): 0.5994,
            ("HG005", "whatshap"): 0.8621, ("HG005", "longphase"): 0.8704}
NAMED = ("PMS2", "SMN1", "SMN2", "GBA1", "CYP21A2", "HBA1", "HBA2", "STRC", "TNXB", "NEB",
         "CYP11B1", "CYP11B2", "NCF1", "PKD1", "CYP2D6", "SDHA", "CHEK2", "OTOA", "IKBKG",
         "HYDIN")
CACHE = ROOT / "data" / "l2"


# ---- pure rules (`tests/test_l2_end_calls.py`) ----

def relation(a, b):
    """1 cis, 0 trans, None unrelated; values are `[allele, block]` or None."""
    if a is None or b is None or a[1] is None or a[1] != b[1]:
        return None
    return int(a[0] == b[0])


def read_agrees(rel, support) -> bool:
    """At least 2 reads call both ends and a strict majority give relation `rel`."""
    if rel is None or support is None:
        return False
    same, diff = support
    return same + diff >= 2 and (same > diff if rel == 1 else diff > same)


def certified(n: int, k: int, r: float, delta: float = DELTA) -> bool:
    """`P(X <= k | n, r) <= delta`: the per-call rate is below `r` at level `delta`."""
    from phasebudget.confidence import binomial_cdf

    return n > 0 and binomial_cdf(k, n, r) <= delta


def p_exceed(k: int, n: int, r: float) -> float:
    from phasebudget.confidence import binomial_cdf

    return 1.0 if k == 0 else 1.0 - binomial_cdf(k - 1, n, r)


def pairs(variants, genes):
    """`(i, j, names)` for `i < j` inside one listed gene, at most `WINDOW` apart; `variants`
    sorted by position (`v[0]`), `genes` `[(start, end, name)]` 1-based inclusive."""
    positions = [v[0] for v in variants]
    holding: list = [[] for _ in variants]
    for start, end, name in genes:
        for k in range(bisect.bisect_left(positions, start), bisect.bisect_right(positions, end)):
            holding[k].append((end, name))
    for i, pos_i in enumerate(positions):
        if not holding[i]:
            continue
        limit = min(pos_i + WINDOW, max(e for e, _ in holding[i]))
        for j in range(i + 1, len(positions)):
            if positions[j] > limit:
                break
            yield i, j, [g for e, g in holding[i] if e >= positions[j]]


# ---- extraction, in the platform's environment ----

def extract(chrom: str) -> dict:
    import chain_relations as cr
    import genome_common as common
    import k1b_score as kb

    stratify = common._module("stratify.py")
    extractor = common._module("extract_fragments.py")
    table = json.loads(kb.table_path(common.PLATFORM, chrom).read_text())
    variants = [tuple(v) for v in table["variants"]]
    truth = common.primary_truth(chrom, variants)
    genes = cr.genes([chrom], cr.clinical_symbols())[chrom]
    positions = [v[0] for v in variants]
    ingene = [False] * len(variants)
    for start, end, _ in genes:
        for k in range(bisect.bisect_left(positions, start), bisect.bisect_right(positions, end)):
            ingene[k] = True
    keep = [k for k in range(len(variants)) if ingene[k] and truth[k] is not None]
    beds = {name: stratify.load_bed(path, [chrom]) for name, path in stratify.fetch().items()}
    zero = [positions[k] - 1 for k in keep]
    named = [stratify.inside(beds[n][chrom], zero) for n in ("segdups", "lowmap", "TR+homopol")]
    easy = stratify.inside(beds["easy"][chrom], zero)
    out_vars = []
    for x, k in enumerate(keep):
        _, ref, alt = variants[k]
        hard = len(ref) != 1 or len(alt) != 1 or any(b[x] for b in named) or not easy[x]
        out_vars.append([positions[k], hard, truth[k], table["whatshap"][k],
                         table["longphase"][k], table["assembly"][k]])
    local = {positions[k] - 1: x for x, k in enumerate(keep)}  # 0-based pos to local index
    matrix, frag_pos, _, _ = extractor.load_fragments(
        str(common.fragments_path(chrom, DEPTH, "realign")))
    col = [local.get(p) for p in frag_pos]
    counts: dict = {}
    for f in matrix:
        calls = [(col[f.start + o], a) for o, a in enumerate(f.alleles)
                 if a is not None and col[f.start + o] is not None]
        for u in range(len(calls)):
            i, ai = calls[u]
            for v in range(u + 1, len(calls)):
                j, aj = calls[v]
                if out_vars[j][0] - out_vars[i][0] > WINDOW:
                    break
                c = counts.setdefault((i, j), [0, 0])
                c[0 if ai == aj else 1] += 1
    return {"variants": out_vars, "reads": [[i, j, s, d] for (i, j), (s, d) in counts.items()]}


def cache_path(platform: str, chrom: str) -> Path:
    return CACHE / f"{platform}_{DEPTH}_{chrom}.json.gz"


def extract_job(job) -> str:
    platform, chrom = job
    path = cache_path(platform, chrom)
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        env = dict(os.environ, PHASEBUDGET_PLATFORM=platform,
                   PYTHONPATH=f"{ROOT / 'src'}:{HERE}")
        env.pop("PHASEBUDGET_TRUTH", None)
        part = path.with_name(path.name + f".{os.getpid()}.part")
        subprocess.run([sys.executable, __file__, "--extract", chrom, str(part)], env=env,
                       check=True)
        part.replace(path)
    return str(path)


# ---- one (genome, chromosome): tallies and picks ----

def tally(job) -> dict:
    """Per (phaser, arm, group): agreeing pairs, wrong agreeing, and per gene the agreeing pair
    with the smallest CRC32 `(crc, wrong)`; per distance bin all pairs; named-gene counts."""
    label, platform, chrom = job
    import chain_relations as cr

    with gzip.open(cache_path(platform, chrom), "rt") as handle:
        data = json.load(handle)
    variants = data["variants"]
    reads = {(i, j): (s, d) for i, j, s, d in data["reads"]}
    genes = cr.genes([chrom], cr.clinical_symbols())[chrom]
    keys = [(ph, arm, g) for ph in PHASERS for arm in ARMS for g in GROUPS]
    agree = {k: 0 for k in keys}
    wrong = {k: 0 for k in keys}
    picks: dict = {}  # (phaser, arm, group, gene): (crc, wrong)
    named: dict = {}  # (gene, None | (phaser, arm, group)): pairs
    bins = {d: 0 for d in DISTANCES}
    for i, j, names in pairs(variants, genes):
        vi, vj = variants[i], variants[j]
        d = DISTANCES[bisect.bisect_left(DISTANCES, vj[0] - vi[0])]
        group = (d, "hard" if vi[1] or vj[1] else "easy")
        bins[d] += 1
        truth = int(vi[2] == vj[2])
        asm = relation(vi[5], vj[5])
        support = reads.get((i, j))
        crc = None
        for col, ph in ((3, "whatshap"), (4, "longphase")):
            rel = relation(vi[col], vj[col])
            if rel is None:
                continue
            bad = rel != truth
            for arm, ok in (("assembly", asm is not None and asm == rel),
                            ("reads", read_agrees(rel, support))):
                if not ok:
                    continue
                k = (ph, arm, group)
                agree[k] += 1
                wrong[k] += bad
                if crc is None:
                    crc = [zlib.crc32(f"{n}:{vi[0]}:{vj[0]}".encode()) for n in names]
                for n, c in zip(names, crc):
                    old = picks.get((ph, arm, group, n))
                    if old is None or c < old[0]:
                        picks[(ph, arm, group, n)] = (c, bad)
                    if n in NAMED:
                        named[(n, k)] = named.get((n, k), 0) + 1
        for n in names:
            if n in NAMED:
                named[(n, None)] = named.get((n, None), 0) + 1
    return {"label": label, "agree": agree, "wrong": wrong, "picks": picks, "named": named,
            "bins": bins}


def merge(parts: list) -> dict:
    out = {"agree": {}, "wrong": {}, "picks": {}, "named": {}, "bins": {d: 0 for d in DISTANCES}}
    for p in parts:
        for key in ("agree", "wrong", "named"):
            for k, v in p[key].items():
                out[key][k] = out[key].get(k, 0) + v
        for d, v in p["bins"].items():
            out["bins"][d] += v
        for (ph, arm, g, n), pick in p["picks"].items():
            k = (ph, arm, g, (p["label"], n))  # a pool pick is per (genome, gene)
            if k not in out["picks"] or pick[0] < out["picks"][k][0]:
                out["picks"][k] = pick
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extract", nargs=2, metavar=("CHROM", "OUT"))
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()
    if args.extract:
        with gzip.open(args.extract[1], "wt") as handle:
            json.dump(extract(args.extract[0]), handle)
        return

    started = time.time()
    odd = [c for c in CHROMS if int(c[3:]) % 2]
    even = [c for c in CHROMS if not int(c[3:]) % 2]
    if args.dry_run:
        pool = [("HG002a", "hprc_HG002", odd[0::2]), ("HG002b", "hprc_HG002", odd[1::2])]
        pool += [(s, p, list(CHROMS)) for s, p in POOL.items() if s != "HG002"]
        tests = [("HG002even", "hprc_HG002", even)]
    else:
        pool = [(s, p, list(CHROMS)) for s, p in POOL.items()]
        tests = [(s, p, list(CHROMS)) for s, p in TESTS.items()]
    mode = "DRY RUN, no test genome read" if args.dry_run else "development, rules stated first"
    print(f"L2 ({mode}, docs/design/l2-end-variant-calls.md): {DEPTH}x realigned "
          f"fragments, phasers {PHASERS}; pool "
          f"{', '.join(f'{s} ({p}, {len(c)} chromosomes)' for s, p, c in pool)}; tests "
          f"{', '.join(f'{s} ({p}, {len(c)} chromosomes)' for s, p, c in tests)}; pairs: "
          f"truth-phased variants in one clinical gene (chain_relations.clinical_symbols), at "
          f"most {WINDOW} bp apart; call = phaser relation in one block; arms: assembly (K1b "
          f"tables, Q6 contig rule, same relation at both ends), reads (>= 2 fragments call both "
          f"ends, strict majority gives the call; descriptive); groups (distance bin "
          f"{DISTANCES}, easy/hard), hard = an indel end or stratum not easy; calibration: one "
          f"CRC32 pick per (pool genome, gene, group) among agreeing pairs, group certified at r "
          f"when P(X <= k | n, r) <= {DELTA}; targets {TARGETS}; {CELLS} cells, level "
          f"{LEVEL:.3e}; L2-2: share of NA12877 pairs 10-100 kb certified at r {GOAL_R:g} >= "
          f"{GOAL_SHARE} for one phaser (K4b {K4B_LONG})", flush=True)

    jobs = sorted({(p, c) for _, p, cs_ in pool + tests for c in cs_})
    with ThreadPoolExecutor(args.workers) as ex:
        list(ex.map(extract_job, jobs))
    print(f"extraction done: {len(jobs)} files ({time.time() - started:.0f} s)", flush=True)
    work = [(label, p, c) for label, p, cs_ in pool + tests for c in cs_]
    with ProcessPoolExecutor(args.workers) as ex:
        parts = list(ex.map(tally, work))
    pool_labels = {label for label, _, _ in pool}
    cal = merge([p for p in parts if p["label"] in pool_labels])

    cert: dict = {}  # (phaser, arm, r): set of groups
    print("\npool calibration: per phaser, arm, group: agreeing pairs, wrong, picks n, wrong "
          "picks k, certified at r")
    for ph in PHASERS:
        for arm in ARMS:
            for g in GROUPS:
                got = [v for (a, b, c, _), v in cal["picks"].items() if (a, b, c) == (ph, arm, g)]
                n, k = len(got), sum(w for _, w in got)
                ok = [r for r in TARGETS if certified(n, k, r)]
                for r in ok:
                    cert.setdefault((ph, arm, r), set()).add(g)
                key = (ph, arm, g)
                print(f"  {ph:<9} {arm:<8} <= {g[0] // 1000:>3} kb {g[1]:<4} agreeing "
                      f"{cal['agree'].get(key, 0):>9} wrong {cal['wrong'].get(key, 0):>5} picks "
                      f"{n:>6} wrong picks {k:>3} certified at {ok or 'none'}")

    verdict_1, goal = True, {}
    for label, _, _ in tests:
        res = merge([p for p in parts if p["label"] == label])
        total = sum(res["bins"].values())
        print(f"\n== test {label}: pairs {total} ("
              + ", ".join(f"<= {d // 1000} kb {n}" for d, n in res["bins"].items()) + ")")
        for ph in PHASERS:
            for arm in ARMS:
                for r in TARGETS:
                    groups = cert.get((ph, arm, r), set())
                    kept = sum(res["agree"].get((ph, arm, g), 0) for g in groups)
                    bad = sum(res["wrong"].get((ph, arm, g), 0) for g in groups)
                    best: dict = {}
                    for (a, b, g, (_, n)), pick in res["picks"].items():
                        if (a, b) == (ph, arm) and g in groups:
                            if n not in best or pick[0] < best[n][0]:
                                best[n] = pick
                    n, k = len(best), sum(w for _, w in best.values())
                    p = p_exceed(k, n, r)
                    violates = arm == "assembly" and n > 0 and p < LEVEL
                    verdict_1 &= not violates
                    by_bin = {d: sum(res["agree"].get((ph, arm, (d, c)), 0) for c in CLASSES
                                     if (d, c) in groups) / max(1, res["bins"][d])
                              for d in DISTANCES}
                    if arm == "assembly" and r == GOAL_R:
                        goal[(label, ph)] = by_bin[GOAL_BIN]
                    print(f"  {ph:<9} {arm:<8} r {r:<6g} certified {kept:>9} share "
                          f"{kept / max(1, total):.4f} wrong certified {bad:>5} genes picked "
                          f"{n:>5} wrong picks {k:>3} p {p:.2e}"
                          f"{'  VIOLATES' if violates else ''}{'' if arm == 'assembly' else ' (descriptive)'}")
                    print("      by distance: " + ", ".join(
                        f"<= {d // 1000} kb {s:.4f}" for d, s in by_bin.items())
                        + (f" (K4b 10-100 kb {K4B_LONG[(label, ph)]:.4f})"
                           if (label, ph) in K4B_LONG else ""))
                    if r == GOAL_R and arm == "assembly":
                        names = sorted({n for n, _ in res["named"]})
                        print("      named genes at r {:g}: ".format(r) + ", ".join(
                            f"{n} {sum(res['named'].get((n, (ph, arm, g)), 0) for g in groups)}"
                            f"/{res['named'].get((n, None), 0)}" for n in names))

    print()
    worst = min(goal.items(), key=lambda kv: kv[1], default=None)
    best = max(goal.items(), key=lambda kv: kv[1], default=None)
    print(f"10-100 kb certified share at r {GOAL_R:g}, assembly arm: "
          + ", ".join(f"{t} {ph} {s:.4f}" for (t, ph), s in goal.items())
          + f"; best {best}, worst {worst}")
    if args.dry_run:
        print(f"dry run: violating cells at {LEVEL:.3e}: {'none' if verdict_1 else 'some'}")
    else:
        ok2 = any(goal.get(("NA12877", ph), 0) >= GOAL_SHARE for ph in PHASERS)
        print(f"L2-1 validity: {'holds' if verdict_1 else 'fails'}")
        print(f"L2-2 NA12877 10-100 kb >= {GOAL_SHARE}: {'holds' if ok2 else 'fails'}")
        print(f"L2: {'succeeds' if verdict_1 and ok2 else 'does not succeed'}")
    print(f"runtime {time.time() - started:.0f} s")


if __name__ == "__main__":
    main()
