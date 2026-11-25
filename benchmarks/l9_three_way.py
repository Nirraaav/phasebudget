"""L9 (`docs/design/l9-three-way-calls.md`): L2's end-variant calls with a three-way rule
(whatshap, LongPhase and the assembly all relate the pair and agree), beside L2's two assembly
arms, with per-gene shares. Reads L2's cached tables in `data/l2/` (extracted by
`l2_end_calls.py`, which this script calls for any missing file).

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/l9_three_way.py`
"""

from __future__ import annotations

import bisect
import gzip
import json
import statistics
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import l2_end_calls as l2  # noqa: E402

ARMS = ("whatshap+asm", "longphase+asm", "three-way")
GOAL_R, GOAL_SHARE = 1e-3, 0.5
LONG = (l2.WINDOW, "easy"), (l2.WINDOW, "hard")


def call(vi, vj, arm):
    """The arm's relation for variants `vi`, `vj` (L2 cache rows), or None if it makes no call."""
    asm = l2.relation(vi[5], vj[5])
    wh, lp = l2.relation(vi[3], vj[3]), l2.relation(vi[4], vj[4])
    rel = {"whatshap+asm": wh, "longphase+asm": lp, "three-way": wh if wh == lp else None}[arm]
    return rel if rel is not None and asm == rel else None


def tally(job) -> dict:
    """Per (arm, group): called pairs, wrong, per gene the smallest-CRC call `(crc, wrong)`, and
    per (arm, group, gene) called pairs; per gene all pairs; per distance bin all pairs."""
    import chain_relations as cr

    label, platform, chrom = job
    with gzip.open(l2.cache_path(platform, chrom), "rt") as handle:
        variants = json.load(handle)["variants"]
    genes = cr.genes([chrom], cr.clinical_symbols())[chrom]
    called, wrong, picks, gene_called, gene_all = {}, {}, {}, {}, {}
    bins = {d: 0 for d in l2.DISTANCES}
    for i, j, names in l2.pairs(variants, genes):
        vi, vj = variants[i], variants[j]
        d = l2.DISTANCES[bisect.bisect_left(l2.DISTANCES, vj[0] - vi[0])]
        g = (d, "hard" if vi[1] or vj[1] else "easy")
        bins[d] += 1
        for n in names:
            gene_all[(n, g)] = gene_all.get((n, g), 0) + 1
        truth = int(vi[2] == vj[2])
        crc = None
        for arm in ARMS:
            rel = call(vi, vj, arm)
            if rel is None:
                continue
            bad = rel != truth
            called[(arm, g)] = called.get((arm, g), 0) + 1
            wrong[(arm, g)] = wrong.get((arm, g), 0) + bad
            if crc is None:
                crc = [zlib.crc32(f"{n}:{vi[0]}:{vj[0]}".encode()) for n in names]
            for n, c in zip(names, crc):
                key = (arm, g, n)
                gene_called[key] = gene_called.get(key, 0) + 1
                if key not in picks or c < picks[key][0]:
                    picks[key] = (c, bad)
    return {"label": label, "called": called, "wrong": wrong, "picks": picks,
            "gene_called": gene_called, "gene_all": gene_all, "bins": bins}


def merge(parts: list) -> dict:
    out = {"called": {}, "wrong": {}, "picks": {}, "gene_called": {}, "gene_all": {},
           "bins": {d: 0 for d in l2.DISTANCES}}
    for p in parts:
        for key in ("called", "wrong", "gene_called", "gene_all"):
            for k, v in p[key].items():
                out[key][k] = out[key].get(k, 0) + v
        for d, v in p["bins"].items():
            out["bins"][d] += v
        for (arm, g, n), pick in p["picks"].items():
            k = (arm, g, (p["label"], n))
            if k not in out["picks"] or pick[0] < out["picks"][k][0]:
                out["picks"][k] = pick
    return out


def gene_shares(res: dict, arm: str, groups: set, bins=None) -> tuple[float, float, int]:
    """(share of genes with a certified pair, median per-gene certified share, genes), over genes
    with a pair in `bins` (all bins when None)."""
    keep = lambda g: bins is None or g[0] in bins  # noqa: E731
    total: dict = {}
    for (n, g), v in res["gene_all"].items():
        if keep(g):
            total[n] = total.get(n, 0) + v
    got = {n: 0 for n in total}
    for (a, g, n), v in res["gene_called"].items():
        if a == arm and g in groups and keep(g) and n in got:
            got[n] += v
    if not total:
        return 0.0, 0.0, 0
    shares = [got[n] / total[n] for n in total]
    return sum(s > 0 for s in shares) / len(shares), statistics.median(shares), len(shares)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()
    started = time.time()
    pool = [(s, p, list(l2.CHROMS)) for s, p in l2.POOL.items()]
    tests = [(s, p, list(l2.CHROMS)) for s, p in l2.TESTS.items()]
    print(f"L9 (development, rules stated first, docs/design/l9-three-way-calls.md): L2's cached "
          f"tables ({l2.DEPTH}x), arms {ARMS} (three-way = whatshap, LongPhase and the assembly "
          f"all relate the pair, each in one block or contig, with one relation); pool "
          f"{', '.join(f'{s} ({p})' for s, p, _ in pool)}; tests "
          f"{', '.join(f'{s} ({p})' for s, p, _ in tests)}; 21 chromosomes; pairs in one clinical "
          f"gene at most {l2.WINDOW} bp apart; groups (bin {l2.DISTANCES}, easy/hard); one CRC32 "
          f"pick per (pool genome, gene, group); certified at r when P(X <= k | n, r) <= "
          f"{l2.DELTA}; targets {l2.TARGETS}; test level {l2.LEVEL:.3e}; L9-1 both 10-100 kb "
          f"groups certify three-way at {GOAL_R:g}; L9-2 NA12877 10-100 kb share >= {GOAL_SHARE} "
          f"three-way at {GOAL_R:g}; L9-3 no three-way test cell below the level", flush=True)
    jobs = sorted({(p, c) for _, p, cs_ in pool + tests for c in cs_})
    with ThreadPoolExecutor(args.workers) as ex:
        list(ex.map(l2.extract_job, jobs))
    work = [(label, p, c) for label, p, cs_ in pool + tests for c in cs_]
    with ProcessPoolExecutor(args.workers) as ex:
        parts = list(ex.map(tally, work))
    cal = merge([p for p in parts if p["label"] in l2.POOL])

    cert: dict = {}
    print("\npool: arm, group: called pairs, wrong, picks n, wrong picks k, certified at r")
    for arm in ARMS:
        for g in l2.GROUPS:
            got = [v for (a, b, _), v in cal["picks"].items() if (a, b) == (arm, g)]
            n, k = len(got), sum(w for _, w in got)
            ok = [r for r in l2.TARGETS if l2.certified(n, k, r)]
            for r in ok:
                cert.setdefault((arm, r), set()).add(g)
            print(f"  {arm:<14} <= {g[0] // 1000:>3} kb {g[1]:<4} called "
                  f"{cal['called'].get((arm, g), 0):>9} wrong {cal['wrong'].get((arm, g), 0):>5} "
                  f"picks {n:>6} wrong picks {k:>3} certified at {ok or 'none'}")
    l9_1 = all(g in cert.get(("three-way", GOAL_R), set()) for g in LONG)

    l9_3, long_share = True, {}
    for label, _, _ in tests:
        res = merge([p for p in parts if p["label"] == label])
        total = sum(res["bins"].values())
        print(f"\n== test {label}: pairs {total} ("
              + ", ".join(f"<= {d // 1000} kb {n}" for d, n in res["bins"].items()) + ")")
        for arm in ARMS:
            for r in l2.TARGETS:
                groups = cert.get((arm, r), set())
                kept = sum(res["called"].get((arm, g), 0) for g in groups)
                bad = sum(res["wrong"].get((arm, g), 0) for g in groups)
                best: dict = {}
                for (a, g, (_, n)), pick in res["picks"].items():
                    if a == arm and g in groups and (n not in best or pick[0] < best[n][0]):
                        best[n] = pick
                n, k = len(best), sum(w for _, w in best.values())
                p = l2.p_exceed(k, n, r)
                bad_cell = n > 0 and p < l2.LEVEL
                if arm == "three-way":
                    l9_3 &= not bad_cell
                by_bin = {d: sum(res["called"].get((arm, (d, c)), 0) for c in l2.CLASSES
                                 if (d, c) in groups) / max(1, res["bins"][d])
                          for d in l2.DISTANCES}
                long_share[(label, arm, r)] = by_bin[l2.WINDOW]
                gs, gm, ng = gene_shares(res, arm, groups)
                ls, lm, nl = gene_shares(res, arm, groups, {l2.WINDOW})
                print(f"  {arm:<14} r {r:<6g} certified {kept:>9} share {kept / max(1, total):.4f}"
                      f" wrong certified {bad:>5} genes picked {n:>5} wrong picks {k:>3} p "
                      f"{p:.2e}{'  VIOLATES' if bad_cell else ''}")
                print("      by distance: " + ", ".join(
                    f"<= {d // 1000} kb {s:.4f}" for d, s in by_bin.items()))
                print(f"      genes: all bins {ng} genes, with a certified pair {gs:.4f}, median "
                      f"per-gene share {gm:.4f}; 10-100 kb {nl} genes, with a certified pair "
                      f"{ls:.4f}, median per-gene share {lm:.4f}")

    l9_2 = long_share.get(("NA12877", "three-way", GOAL_R), 0.0) >= GOAL_SHARE
    print("\n10-100 kb share at r 1e-3: " + ", ".join(
        f"{t} {a} {long_share[(t, a, GOAL_R)]:.4f}" for t, _, _ in tests for a in ARMS)
          + "; worst three-way " + f"{min(long_share[(t, 'three-way', GOAL_R)] for t, _, _ in tests):.4f}")
    print(f"L9-1 pool, both 10-100 kb groups certify three-way at {GOAL_R:g}: "
          f"{'holds' if l9_1 else 'fails'}")
    print(f"L9-2 NA12877 10-100 kb three-way share at {GOAL_R:g} >= {GOAL_SHARE}: "
          f"{'holds' if l9_2 else 'fails'}")
    print(f"L9-3 validity: {'holds' if l9_3 else 'fails'}")
    print(f"L9: {'succeeds' if l9_1 and l9_2 and l9_3 else 'does not succeed'}")
    print(f"runtime {time.time() - started:.0f} s")


if __name__ == "__main__":
    main()
