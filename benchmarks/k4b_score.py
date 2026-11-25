"""K4b (registered, `docs/design/k4b-clinical-agreement-plan.md`): within-gene cis or trans calls
whose links must agree with the person's assembly, thresholds from K1b's six-genome pool.

The pool's rows are K1b's (`k1b_score.pair_rows` on `data/k1/tables/`). Test calls are Q4's
(`q4_clinical_cis_trans`), on Q4 chain files extracted for the test platforms, with every link
the test genome's assembly does not agree with set below every threshold
(`k1_assembly_agreement.rewrite_links`).

`--dry-run` reads no test genome: HG002's odd chromosomes, split in two pseudo-genomes, and the
five Pedigree genomes form the pool, and HG002's even chromosomes stand in for a test genome.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/k4b_score.py [--dry-run]`
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import k1_assembly_agreement as k1  # noqa: E402
import k1b_score as kb  # noqa: E402
import q4_clinical_cis_trans as q4  # noqa: E402

ROOT = HERE.parent
DEPTH, CHROMS, PHASERS = kb.DEPTH, kb.CHROMS, kb.PHASERS
CELLS = len(kb.TESTS) * len(PHASERS) * len(q4.TARGETS)  # 8
LEVEL = 0.05 / CELLS
USEFUL_R, USEFUL_SHARE = 1e-2, 0.5


def assembly_of(table: dict) -> dict:
    """`{pos1: [allele, block]}` of a K1b table's assembly values."""
    return {v[0]: a for v, a in zip(table["variants"], table["assembly"])}


def tally(platform: str, phaser: str, chroms, thresholds: dict, assembly: dict | None,
          stratum, genes: dict) -> dict:
    """Q4's tallies for one test platform and phaser; `assembly` (`{chrom: {pos1: value}}`)
    applies the agreement rule, None leaves the links as they are."""
    import confidence_score as cs

    nt = len(q4.TARGETS)
    out = {"pairs": 0, "chain": 0, "certified": [0] * nt, "wrong": [0] * nt,
           "bins": {d: [0, 0] for d in q4.DISTANCES}, "named": {}}
    chosen: dict = {}
    useful_x = q4.TARGETS.index(USEFUL_R)
    for chrom in chroms:
        path = q4._extract_job((ROOT, platform, DEPTH, phaser, chrom))
        with gzip.open(path, "rt") as handle:
            raw = json.load(handle)
        if assembly is not None:
            raw = k1.rewrite_links(raw, assembly[chrom])
        phased = [k for k, v in enumerate(raw) if v[1] is not None]
        nxt = dict(zip(phased, phased[1:]))
        variants = []
        for k, (pos, allele, truth, key, indel) in enumerate(raw):
            link = None
            if key is not None:
                other = raw[nxt[k]]
                group = "indel" if indel or other[4] else stratum(chrom, pos, other[0])
                t = thresholds.get(group, [None] * len(q4.BETAS))
                link = (q4.link_level(tuple(key), t), (group, pos // cs.UNIT), group)
            variants.append((pos, allele, truth, link))
        reach, holding = q4.gene_index(variants, genes[chrom])
        for i, j, status, cells, level, _, _ in q4.pairs(variants, reach, q4.WINDOW):
            pos_i, pos_j = variants[i][0], variants[j][0]
            names = [g for e, g in holding[i] if e >= pos_j]
            dist = next(d for d in q4.DISTANCES if pos_j - pos_i <= d)
            out["pairs"] += 1
            out["bins"][dist][0] += 1
            for name in names:
                if name in q4.NAMED:
                    out["named"].setdefault(name, [0, 0])[0] += 1
            if status != "chain":
                continue
            out["chain"] += 1
            wrong = (variants[i][1] ^ variants[i][2]) != (variants[j][1] ^ variants[j][2])
            for x, r in enumerate(q4.TARGETS):
                b = q4.b_max(len(cells), r)
                if b is None or level > b:
                    continue
                out["certified"][x] += 1
                out["wrong"][x] += wrong
                if x == useful_x:
                    out["bins"][dist][1] += 1
                    for name in names:
                        if name in q4.NAMED:
                            out["named"][name][1] += 1
                for name in names:
                    crc = zlib.crc32(f"{name}:{pos_i}:{pos_j}".encode())
                    old = chosen.get((name, x))
                    if old is None or crc < old[0]:
                        chosen[(name, x)] = (crc, wrong)
    out["validity"] = [(len(p), sum(p)) for p in
                       ([w for (_, y), (_, w) in chosen.items() if y == x] for x in range(nt))]
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()

    import confidence_score as cs
    from stratified_gate import strata

    started = time.time()
    odd = [c for c in CHROMS if int(c[3:]) % 2]
    even = [c for c in CHROMS if not int(c[3:]) % 2]
    if args.dry_run:
        pool = [("HG002a", "hprc_HG002", odd[0::2]), ("HG002b", "hprc_HG002", odd[1::2])]
        pool += [(s, p, list(CHROMS)) for s, p in kb.POOL.items() if s != "HG002"]
        tests = [("HG002even", "hprc_HG002", even)]
    else:
        pool = [(s, p, list(CHROMS)) for s, p in kb.POOL.items()]
        tests = [(s, p, list(CHROMS)) for s, p in kb.TESTS.items()]
    print(f"K4b ({'DRY RUN, no test genome read' if args.dry_run else 'registered'}, "
          f"docs/design/k4b-clinical-agreement-plan.md): {DEPTH}x, whatshap and LongPhase; pool "
          f"{', '.join(f'{s} ({p}, {len(c)} chromosomes)' for s, p, c in pool)}; tests "
          f"{', '.join(f'{s} ({p}, {len(c)} chromosomes)' for s, p, c in tests)}; thresholds "
          f"from pooled agreeing pairs (baseline: all pooled pairs), betas {q4.BETAS}, delta "
          f"{cs.DELTA}, unit (genome, chromosome, {cs.UNIT // 1000} kb); Q4 call on clinical "
          f"genes, window {q4.WINDOW}, targets {q4.TARGETS}; links the test genome's assembly "
          f"does not agree with are never trusted; {CELLS} registered cells, level {LEVEL:.3e}; "
          f"K4b-2/3 certified share >= {USEFUL_SHARE} at r {USEFUL_R:g} for one phaser",
          flush=True)

    stratum = strata(list(CHROMS))
    genes = q4.clinical_genes(list(CHROMS))
    jobs = sorted({(p, c) for _, p, cs_ in pool + tests for c in cs_})
    with ThreadPoolExecutor(args.workers) as ex:
        loaded = dict(zip(jobs, ex.map(kb.load_table, jobs)))
        # Chains, extracted in parallel before the single-threaded tallies.
        list(ex.map(q4._extract_job, [(ROOT, p, DEPTH, ph, c) for _, p, cs_ in tests
                                      for ph in PHASERS for c in cs_]))

    verdict_1, useful = True, {}
    for phaser in PHASERS:
        cal = []
        for label, platform, chroms in pool:
            rows = cs.load_rows(platform, DEPTH, phaser, list(CHROMS))
            cal += kb.pair_rows(label, rows, {c: loaded[(platform, c)] for c in chroms},
                                phaser, stratum)
        thresholds = {"agree": q4.thresholds_for([r[:4] for r in cal if r[4]]),
                      "all": q4.thresholds_for([r[:4] for r in cal])}
        print(f"\n== {phaser}: pool pairs {len(cal)}, agreeing {sum(r[4] for r in cal)}")
        for mode, t in thresholds.items():
            print(f"  {mode:<5} thresholds per group (beta {q4.BETAS}): " + "; ".join(
                f"{g} " + " ".join("-" if x is None else f"{x[0]:.3f}/{x[1]}" for x in ts)
                for g, ts in sorted(t.items())))
        for label, platform, chroms in tests:
            assembly = {c: assembly_of(loaded[(platform, c)]) for c in chroms}
            for mode in ("agree", "all"):
                res = tally(platform, phaser, chroms, thresholds[mode],
                            assembly if mode == "agree" else None, stratum, genes)
                print(f"\n  test {label} {phaser} {mode}: pairs {res['pairs']}, unbroken chains "
                      f"{res['chain']}")
                for x, r in enumerate(q4.TARGETS):
                    n, k = res["validity"][x]
                    p = q4.p_exceed(k, n, r)
                    share = res["certified"][x] / max(1, res["pairs"])
                    bad = mode == "agree" and n > 0 and p < LEVEL
                    verdict_1 &= not bad
                    if mode == "agree" and r == USEFUL_R:
                        useful.setdefault(label, []).append((phaser, share))
                    print(f"    r {r:<6g} certified {res['certified'][x]:>8} share {share:.4f} "
                          f"wrong certified {res['wrong'][x]:>4} genes picked {n:>5} wrong picks "
                          f"{k:>3} p {p:.2e}{'  VIOLATES' if bad else ''}")
                print(f"    certified share at r {USEFUL_R:g} by distance: " + ", ".join(
                    f"<= {d // 1000} kb {c / max(1, n):.4f} ({c}/{n})"
                    for d, (n, c) in res["bins"].items()))
                if res["named"]:
                    print("    named paralog genes at r {:g}: ".format(USEFUL_R) + ", ".join(
                        f"{g} {c}/{n}" for g, (n, c) in sorted(res["named"].items())))

    print()
    if args.dry_run:
        print(f"dry run: violating cells at {LEVEL:.3e}: {'none' if verdict_1 else 'some'}; "
              f"certified share at r {USEFUL_R:g}: {useful}")
    else:
        print(f"K4b-1 validity: {'supported' if verdict_1 else 'not supported'}")
        for k, label in enumerate(kb.TESTS, start=2):
            got = useful.get(label, [])
            ok = any(s >= USEFUL_SHARE for _, s in got)
            print(f"K4b-{k} usefulness on {label}: " + ", ".join(
                f"{ph} {s:.4f}" for ph, s in got) + f": {'supported' if ok else 'not supported'}")
    print(f"runtime {time.time() - started:.0f} s")


if __name__ == "__main__":
    main()
