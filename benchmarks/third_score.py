"""Decide V1 to V4 and V6 of `docs/design/third-test-plan.md` and report V5 and V7.

Thresholds come from `learn_then_test_by_group` on `cal_v5` (HG002, sup v5.0.0), grouped as
in `stratified_score.py`; `hg001_v5` (HG001, sup v5.0.0) is scored against the Platinum
Pedigree truth. A cell is (depth, phaser, beta, group). V7 repeats the cells with thresholds
from `cal_r10` (HG002, the 2023 sup reads) and is descriptive.

With `--dry-run` no test genome is read: `--calibration` (default `cal_r10`) calibrates on its
odd chromosomes and is scored on its even ones. That run checks the code, not the hypotheses.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/third_score.py [--dry-run]`
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import confidence_score as cs  # noqa: E402
import genome_common as common  # noqa: E402
from stratified_score import DIFFICULT, cells, fewer_errors, grouped, s7  # noqa: E402

from phasebudget.confidence import learn_then_test_by_group  # noqa: E402

CAL, TEST, OLD_CAL = "cal_v5", "hg001_v5", "cal_r10"
# Levels are fixed by the registered cell counts, whatever number of cells has test pairs.
V1_CELLS, V2_CELLS = 3 * 2 * 3, 3 * 2 * 3 * len(DIFFICULT)
V3_DEPTHS, V3_BETA, V3_SHARE = ("30", "20"), 1e-3, 0.99


def main() -> None:
    import argparse

    from stratified_gate import strata

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--calibration", default=OLD_CAL, help="dry run only")
    args = parser.parse_args()
    chroms = list(common.CONFIRMATORY)
    odd = [c for c in chroms if int(c[3:]) % 2]
    even = [c for c in chroms if int(c[3:]) % 2 == 0]
    stratum = strata(chroms)
    if args.dry_run:
        cal, old, test, cal_chroms, test_chroms = (args.calibration, None, args.calibration,
                                                   odd, even)
    else:
        cal, old, test, cal_chroms, test_chroms = CAL, OLD_CAL, TEST, chroms, chroms
    print(f"third test scoring{' (DRY RUN: ' + cal + ' only)' if args.dry_run else ''}; "
          f"calibration {cal}, test {test}; units {cs.UNIT // 1000} kb per group; betas "
          f"{cs.BETAS}; delta {cs.DELTA}; depths {cs.DEPTHS}; phasers {cs.PHASERS}")
    results, old_results, v4 = [], [], []
    for depth in cs.DEPTHS:
        for phaser in cs.PHASERS:
            test_rows = grouped(test, depth, phaser, test_chroms, stratum)
            for name, source, out in ((cal, cal, results), (old, old, old_results)):
                if source is None:
                    continue
                cal_rows = grouped(source, depth, phaser, cal_chroms, stratum)
                for beta in cs.BETAS:
                    thresholds = learn_then_test_by_group(cal_rows, beta, cs.DELTA)
                    for group, cell in cells(thresholds, test_rows, beta).items():
                        out.append((name, depth, phaser, beta, group, cell))
                    if source == cal and depth in V3_DEPTHS and beta == V3_BETA:
                        pooled = learn_then_test_by_group([("all",) + r[1:] for r in cal_rows],
                                                          beta, cs.DELTA)
                        per = cells(thresholds, test_rows, beta)
                        one = cells(pooled, [("all",) + r[1:] for r in test_rows], beta)["all"]
                        share = sum(c[4] for c in per.values()) / sum(c[3] for c in per.values())
                        v4.append((depth, phaser, share, one[4] / max(1, one[3])))
    print("\ncalibration  depth phaser     beta    group      units  wrong  p_exceed  pairs     "
          "trusted  wrong  wrong_trusted")
    for name, depth, phaser, beta, group, c in results + old_results:
        print(f"{name:<12} {depth:>5} {phaser:<10} {beta:<7g} {group:<10} {c[0]:>6} {c[1]:>6}"
              f"  {c[2]:<8.2e}  {c[3]:>8} {c[4]:>8} {c[5]:>6} {c[6]:>14}")

    def decide(name, chosen, registered):
        level = 0.05 / registered
        bad = [r for r in chosen if r[5][0] and r[5][2] < level]
        certifying = sum(1 for r in chosen if r[5][7])
        print(f"{name}: {len(chosen)} cells with test pairs of {registered} registered, "
              f"{certifying} certifying, level {level:.2e}, violating {len(bad)}, smallest p "
              f"{min((r[5][2] for r in chosen), default=1.0):.3g}"
              + "".join(f"\n  {r[1]}x {r[2]} beta {r[3]:g} {r[4]} p {r[5][2]:.2e}" for r in bad))
        return bad

    print()
    v1 = [r for r in results if r[4] == "easy"]
    v2 = [r for r in results if r[4] in DIFFICULT]
    bad1 = decide("V1 (easy)", v1, V1_CELLS)
    bad2 = decide("V2 (difficult and indel groups)", v2, V2_CELLS)
    v3 = [r for r in v1 if r[1] in V3_DEPTHS and r[3] == V3_BETA]
    ok3 = [r for r in v3 if r[5][3] and r[5][4] / r[5][3] >= V3_SHARE and fewer_errors(r[5])]
    print(f"V3 (easy, beta {V3_BETA:g}, depths {V3_DEPTHS}): {len(ok3)} of {len(v3)} cells succeed")
    ok4 = [r for r in v4 if r[2] > r[3]]
    print("V4 (per-group trusts more than pooled, beta 1e-3): "
          + "; ".join(f"{d}x {ph} {x:.4f} vs {y:.4f}" for d, ph, x, y in v4)
          + f"; {len(ok4)} of {len(v4)} cells")
    abstain = {}
    for name, depth, phaser, beta, group, c in results:
        a, n = abstain.get((depth, phaser, beta), (0, 0))
        abstain[(depth, phaser, beta)] = (a + (0 if c[7] else c[3]), n + c[3])
    print("V5 abstention (share of pairs in groups certifying nothing), beta 1e-3: "
          + "; ".join(f"{k[0]}x {k[1]} {a / n:.3f}" for k, (a, n) in abstain.items()
                      if k[2] == 1e-3))
    r6 = s7((test,), test_chroms)
    print("V6 (whatshap worse on shared indel pairs, 30x, level 0.05): "
          + "; ".join(f"{t} shared {n} only-whatshap {w} only-LongPhase {lp} p {p:.2e}"
                      for t, n, w, lp, p in r6))
    if old_results:
        o1 = [r for r in old_results if r[4] == "easy"]
        o2 = [r for r in old_results if r[4] in DIFFICULT]
        print("V7 (descriptive, thresholds from cal_r10):")
        decide("  easy", o1, V1_CELLS)
        decide("  difficult and indel", o2, V2_CELLS)

    def verdict(n_ok, n, partly):
        return ("supported" if n and n_ok == n else
                "partly supported" if n and n_ok >= partly else "not supported")

    print(f"\nV1 {'supported' if v1 and not bad1 else 'not supported'}; V2 "
          f"{'supported' if v2 and not bad2 else 'not supported'}; V3 "
          f"{verdict(len(ok3), len(v3), 3)}; V4 {verdict(len(ok4), len(v4), 3)}; V6 "
          + ("supported" if r6 and r6[0][4] < 0.05 else "not supported"))


if __name__ == "__main__":
    main()
