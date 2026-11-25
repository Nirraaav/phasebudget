"""Decide B1 to B4 of `docs/design/w1b-pooled-plan.md` and report B5 and B6.

Per depth and phaser, per-group thresholds come from `learn_then_test_by_group`
on the pooled pairs of HG002 (v5.0q) and the first 20 consensus-truth HPRC
genomes, a unit being (genome, chromosome, 100 kb index). Test genomes: five
Platinum Pedigree members' HiFi pairs against the Pedigree truth (B1, B4, B5),
and the other 19 consensus-truth HPRC genomes (B2, testable cells only, by W1
amendment 1's gate). B3 asks whether `segdup` certifies at beta 1e-3 at
depths 30 and 20. B6 calibrates on the same pooled rows with a unit
identified by locus only, wrong when any genome is wrong there.

With `--dry-run` no genome but HG002 is read: HG002's odd chromosomes,
alternately assigned to two pseudo-genomes, form the pool, and its even
chromosomes stand in for every test genome. That run checks the code, not the
hypotheses.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/w1b_score.py [--dry-run]`
"""

from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import confidence_score as cs  # noqa: E402
import genome_common as common  # noqa: E402
from stratified_score import DIFFICULT, cells, fewer_errors, grouped  # noqa: E402
from w1_score import CAL, finished, prefetch  # noqa: E402

from phasebudget.confidence import learn_then_test_by_group  # noqa: E402

POOL = ("HG00735", "HG01243", "HG03516", "HG00438", "HG02622", "HG02145", "HG02723", "HG01358",
        "HG02148", "NA20129", "HG01071", "HG03453", "HG01978", "HG03579", "HG02818", "HG01109",
        "HG00741", "HG02572", "HG01175", "HG00621")
HELD_OUT = ("HG02055", "HG01106", "HG00673", "HG02630", "HG02080", "HG01928", "HG01361",
            "NA18906", "HG01891", "HG03540", "HG01258", "HG02717", "HG02257", "HG01952",
            "HG02886", "HG00733", "NA19240", "HG03486", "HG03098")
PEDIGREE = ("NA12879", "NA12881", "NA12882", "NA12885", "NA12886")
GROUPS = ("easy", *DIFFICULT)
B1_CELLS = len(PEDIGREE) * 3 * 2 * 3 * len(GROUPS)
B2_CELLS = len(HELD_OUT) * 3 * 2 * 3 * len(GROUPS)
B3_DEPTHS, B3_BETA = ("30", "20"), 1e-3
B4_SHARE = 0.99


def pooled(rows_by_genome: dict) -> list:
    """`(group, (genome, unit), score, wrong)`: one genome's unit never merges with another's."""
    return [(g, (genome, u), k, w) for genome, rows in rows_by_genome.items() for g, u, k, w in rows]


def by_locus(rows_by_genome: dict) -> list:
    """B6: the same rows with a unit identified by locus only."""
    return [r for rows in rows_by_genome.values() for r in rows]


def cell(item) -> list:
    """`[(kind, test, depth, phaser, beta, group, cell or threshold)]` for one depth and phaser."""
    from stratified_gate import strata

    depth, phaser, pool, tests, cal_chroms, test_chroms = item
    chroms = list(common.CONFIRMATORY)
    stratum = strata(chroms)
    rows = {name: grouped(platform, depth, phaser, part, stratum)
            for name, (platform, part) in pool.items()}
    single = grouped(CAL, depth, phaser, cal_chroms, stratum)
    test_rows = {t: grouped(t, depth, phaser, test_chroms, stratum) for t in tests}
    out = []
    for beta in cs.BETAS:
        thresholds = learn_then_test_by_group(pooled(rows), beta, cs.DELTA)
        locus = learn_then_test_by_group(by_locus(rows), beta, cs.DELTA)
        alone = learn_then_test_by_group(single, beta, cs.DELTA)
        for group in GROUPS:
            out.append(("threshold", "pooled", depth, phaser, beta, group, thresholds.get(group)))
            out.append(("threshold", "locus", depth, phaser, beta, group, locus.get(group)))
        for t, r in test_rows.items():
            for group, c in cells(thresholds, r, beta).items():
                out.append(("pooled", t, depth, phaser, beta, group, c))
            for group, c in cells(alone, r, beta).items():
                out.append(("alone", t, depth, phaser, beta, group, c))
    print(f"  {depth}x {phaser} scored", flush=True)
    return out


def main() -> None:
    import argparse

    import w1_consensus_gate
    from stratified_gate import strata

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    chroms = list(common.CONFIRMATORY)
    odd = [c for c in chroms if int(c[3:]) % 2]
    even = [c for c in chroms if int(c[3:]) % 2 == 0]
    print(f"W1b scoring{' (DRY RUN: HG002 only)' if args.dry_run else ''}; units "
          f"{cs.UNIT // 1000} kb per group and genome; betas {cs.BETAS}; delta {cs.DELTA}; depths "
          f"{cs.DEPTHS}; phasers {cs.PHASERS}", flush=True)
    testable = w1_consensus_gate.decide(chroms, strata(chroms))
    if args.dry_run:
        pool = {"HG002a": (CAL, odd[0::2]), "HG002b": (CAL, odd[1::2])}
        pedigree, held = (CAL,), (CAL,)
        cal_chroms, test_chroms = odd, even
    else:
        names = [s for s in POOL if finished(f"hprc_{s}")]
        pool = {"HG002": (CAL, chroms), **{s: (f"hprc_{s}", chroms) for s in names}}
        pedigree = tuple(f"phifi_{s}" for s in PEDIGREE if finished(f"phifi_{s}"))
        held = tuple(f"hprc_{s}" for s in HELD_OUT if finished(f"hprc_{s}"))
        cal_chroms = test_chroms = chroms
        print(f"\npool: {len(pool)} genomes ({' '.join(pool)}); Pedigree tests: {len(pedigree)}; "
              f"held-out HPRC tests: {len(held)}", flush=True)
        prefetch([CAL, *(p for p, _ in pool.values()), *pedigree, *held], chroms)
    tests = tuple(dict.fromkeys(pedigree + held))
    items = [(d, ph, pool, tests, cal_chroms, test_chroms) for d in cs.DEPTHS for ph in cs.PHASERS]
    with ProcessPoolExecutor(len(items)) as executor:
        out = [r for part in executor.map(cell, items) for r in part]
    thresholds = [r for r in out if r[0] == "threshold"]
    results = [r[1:] for r in out if r[0] == "pooled"]
    alone = [r[1:] for r in out if r[0] == "alone"]
    print("\ntest              depth phaser     beta    group      units  wrong  p_exceed  pairs"
          "     trusted  wrong  wrong_trusted")
    for t, depth, phaser, beta, group, c in results:
        print(f"{t:<17} {depth:>5} {phaser:<10} {beta:<7g} {group:<10} {c[0]:>6} {c[1]:>6}"
              f"  {c[2]:<8.2e}  {c[3]:>8} {c[4]:>8} {c[5]:>6} {c[6]:>14}")

    def decide(name, chosen, registered):
        level = 0.05 / registered
        bad = [r for r in chosen if r[5][0] and r[5][2] < level]
        certifying = sum(1 for r in chosen if r[5][7])
        print(f"{name}: {len(chosen)} cells with test pairs of {registered} registered, "
              f"{certifying} certifying, level {level:.2e}, violating {len(bad)}, smallest p "
              f"{min((r[5][2] for r in chosen), default=1.0):.3g}"
              + "".join(f"\n  {r[0]} {r[1]}x {r[2]} beta {r[3]:g} {r[4]} p {r[5][2]:.2e} "
                        f"({r[5][1]}/{r[5][0]})" for r in bad))
        return bad

    print()
    b1 = [r for r in results if r[0] in pedigree]
    b2 = [r for r in results if r[0] in held and r[3] in testable.get(r[4], [])]
    bad1 = decide("B1 (Pedigree HiFi, all groups)", b1, B1_CELLS)
    bad2 = decide("B2 (held-out HPRC, testable cells)", b2, B2_CELLS)
    b3 = [r for r in thresholds if r[1] == "pooled" and r[5] == "segdup" and r[4] == B3_BETA
          and r[2] in B3_DEPTHS]
    print(f"B3 (segdup certifies at beta {B3_BETA:g}, depths {B3_DEPTHS}): "
          + "; ".join(f"{r[2]}x {r[3]} {'certifies' if r[6] else 'abstains'}" for r in b3))
    b4 = [r for r in b1 if r[4] == "easy" and r[1] in B3_DEPTHS and r[3] == B3_BETA]
    ok4 = [r for r in b4 if r[5][3] and r[5][4] / r[5][3] >= B4_SHARE and fewer_errors(r[5])]
    print(f"B4 (Pedigree, easy, beta {B3_BETA:g}, depths {B3_DEPTHS}): {len(ok4)} of {len(b4)} "
          "cells succeed")
    print("\nB5 (descriptive), Pedigree cells at beta 1e-3, pooled over genomes, depths and "
          "phasers: trusted share per group, HG002 alone against the pool; share of pairs in "
          "groups certifying nothing")
    for label, source in (("HG002 alone", alone), ("pool", results)):
        mine = [r for r in source if r[0] in pedigree and r[3] == 1e-3]
        shares = "; ".join(
            f"{g} {sum(r[5][4] for r in mine if r[4] == g) / max(1, sum(r[5][3] for r in mine if r[4] == g)):.4f}"
            for g in GROUPS)
        none = sum(r[5][3] for r in mine if not r[5][7]) / max(1, sum(r[5][3] for r in mine))
        print(f"  {label:<12} {shares}; abstaining {none:.3f}")
    print("\nB6 (descriptive): (group, beta) certifying, per depth and phaser, pooled units "
          "against locus-merged units")
    for depth in cs.DEPTHS:
        for phaser in cs.PHASERS:
            for kind in ("pooled", "locus"):
                got = [f"{r[5]} {r[4]:g}" for r in thresholds
                       if r[1] == kind and r[2] == depth and r[3] == phaser and r[6]]
                print(f"  {depth}x {phaser:<9} {kind:<6}: {', '.join(got) or 'none'}")
    print(f"\nB1 {'supported' if b1 and not bad1 else 'not supported'}; B2 "
          f"{'supported' if b2 and not bad2 else 'not supported' if b2 else 'no testable cell'}; B3 "
          + ("supported" if b3 and all(r[6] for r in b3) else "not supported")
          + "; B4 " + ("supported" if b4 and len(ok4) == len(b4) else
                       "partly supported" if len(ok4) >= 18 else "not supported"))


if __name__ == "__main__":
    main()
