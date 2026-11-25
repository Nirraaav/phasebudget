"""Decide A1 to A4 of `docs/design/w1-ancestry-plan.md` (with amendment 1) and report A5 and A6.

Per depth and phaser, per-group thresholds come from `learn_then_test_by_group`
on HG002's HPRC HiFi pairs against v5.0q (`hprc_HG002`, 21 chromosomes); each
test genome's pairs are scored against the consensus truth of dipcall and
1000 Genomes transmission. A cell is (genome, depth, phaser, beta, group). A1
and A2 count only the (group, beta) cells that amendment 1's gate makes
testable (`w1_consensus_gate.decide`, recomputed here), at the registered
levels 0.05/792 and 0.05/3,960. A test genome is scored when its father and
mother are in the 1000 Genomes pedigree and all 21 chromosomes finished.

With `--dry-run` no test genome is read: HG002's odd chromosomes calibrate and
its even chromosomes, against v5.0q, stand in for every test genome. That run
checks the code, not the hypotheses.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/w1_score.py [--dry-run]`
"""

from __future__ import annotations

import gzip
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import confidence_score as cs  # noqa: E402
import genome_common as common  # noqa: E402
from stratified_score import DIFFICULT, cells, fewer_errors, grouped  # noqa: E402

from phasebudget.confidence import learn_then_test_by_group  # noqa: E402

CAL = "hprc_HG002"
A1_CELLS, A2_CELLS = 44 * 3 * 2 * 3, 44 * 3 * 2 * 3 * len(DIFFICULT)
A3_DEPTHS, A3_BETA, A3_SHARE, A3_GENOMES = ("30", "20"), 1e-3, 0.99, 0.9
A4_MARGIN, A4_ALPHA = 0.02, 0.05


def candidates() -> dict:
    """`{sample: superpopulation}` from `w1-candidates.tsv`, HG005 left out."""
    with open(common.ROOT / "docs" / "design" / "w1-candidates.tsv") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in handle]
    return {r["sample"]: r["superpopulation"] for r in rows if r["sample"] != "HG005"}


def finished(platform: str) -> bool:
    """Every confirmatory chromosome passed `confidence_run.py`'s last step."""
    done = common.ROOT / "data" / f"ext_{platform}" / "done"
    return all((done / f"{c}__cleanup.ok").exists() for c in common.CONFIRMATORY)


def prefetch(platforms, chroms) -> None:
    """Compute the row caches in parallel; each `load_rows` runs its own subprocess."""
    jobs = [(p, d, ph) for p in platforms for d in cs.DEPTHS for ph in cs.PHASERS]
    with ThreadPoolExecutor(24) as pool:
        list(pool.map(lambda j: cs.load_rows(j[0], j[1], j[2], chroms), jobs))


def welch_tost(a: list, b: list, margin: float) -> tuple:
    """`(difference, p)` of two one-sided Welch t tests that mean(a) - mean(b) lies within
    plus or minus `margin`; p is the larger of the two."""
    from statistics import mean, variance

    from scipy.stats import t as tdist

    if len(a) < 2 or len(b) < 2:
        return float("nan"), 1.0
    va, vb = variance(a) / len(a), variance(b) / len(b)
    diff, se = mean(a) - mean(b), (va + vb) ** 0.5
    if se == 0:
        return diff, 0.0 if abs(diff) < margin else 1.0
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    p_low = 1 - tdist.cdf((diff + margin) / se, df)
    p_high = tdist.cdf((diff - margin) / se, df)
    return diff, max(p_low, p_high)


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
    stratum = strata(chroms)
    population = candidates()
    print(f"W1 scoring{' (DRY RUN: HG002 only)' if args.dry_run else ''}; calibration {CAL} "
          f"(v5.0q); units {cs.UNIT // 1000} kb per group; betas {cs.BETAS}; delta {cs.DELTA}; "
          f"depths {cs.DEPTHS}; phasers {cs.PHASERS}", flush=True)
    testable = w1_consensus_gate.decide(chroms, stratum)
    if args.dry_run:
        tests, cal_chroms, test_chroms = {CAL: "dry run"}, odd, even
    else:
        tests = {f"hprc_{s}": sp for s, sp in population.items()
                 if common.kg_parents(s) and finished(f"hprc_{s}")}
        cal_chroms = test_chroms = chroms
        left = sorted(s for s in population if f"hprc_{s}" not in tests)
        print(f"\n{len(tests)} test genomes scored; not scored: {' '.join(left)}", flush=True)
        prefetch([CAL, *tests], chroms)
    results, a6 = [], []  # (test, depth, phaser, beta, group, cell)
    for depth in cs.DEPTHS:
        for phaser in cs.PHASERS:
            cal_rows = grouped(CAL, depth, phaser, cal_chroms, stratum)
            os.environ["PHASEBUDGET_TRUTH"] = "dipcall"
            dip_rows = grouped(CAL, depth, phaser, cal_chroms, stratum)
            del os.environ["PHASEBUDGET_TRUTH"]
            test_rows = {t: grouped(t, depth, phaser, test_chroms, stratum) for t in tests}
            for beta in cs.BETAS:
                for rows_from, out in ((cal_rows, results), (dip_rows, a6)):
                    thresholds = learn_then_test_by_group(rows_from, beta, cs.DELTA)
                    for t, rows in test_rows.items():
                        for group, cell in cells(thresholds, rows, beta).items():
                            out.append((t, depth, phaser, beta, group, cell))
            print(f"  {depth}x {phaser} scored", flush=True)
    print("\ntest              depth phaser     beta    group      units  wrong  p_exceed  pairs"
          "     trusted  wrong  wrong_trusted")
    for t, depth, phaser, beta, group, c in results:
        print(f"{t:<17} {depth:>5} {phaser:<10} {beta:<7g} {group:<10} {c[0]:>6} {c[1]:>6}"
              f"  {c[2]:<8.2e}  {c[3]:>8} {c[4]:>8} {c[5]:>6} {c[6]:>14}")

    def decide(name, chosen, registered):
        level = 0.05 / registered
        bad = [r for r in chosen if r[5][0] and r[5][2] < level]
        certifying = sum(1 for r in chosen if r[5][7])
        print(f"{name}: {len(chosen)} testable cells with test pairs of {registered} registered, "
              f"{certifying} certifying, level {level:.2e}, violating {len(bad)}, smallest p "
              f"{min((r[5][2] for r in chosen), default=1.0):.3g}"
              + "".join(f"\n  {r[0]} {r[1]}x {r[2]} beta {r[3]:g} {r[4]} p {r[5][2]:.2e} "
                        f"({r[5][1]}/{r[5][0]})" for r in bad))
        return bad

    def is_testable(r):
        return r[3] in testable.get(r[4], [])

    print()
    a1 = [r for r in results if r[4] == "easy" and is_testable(r)]
    a2 = [r for r in results if r[4] in DIFFICULT and is_testable(r)]
    bad1 = decide("A1 (easy)", a1, A1_CELLS)
    bad2 = decide("A2 (difficult and indel groups)", a2, A2_CELLS)
    verdicts = [f"A1 {'supported' if a1 and not bad1 else 'not supported' if a1 else 'no testable cell'}",
                f"A2 {'supported' if a2 and not bad2 else 'not supported' if a2 else 'no testable cell'}"]
    if A3_BETA in testable.get("easy", []):
        ok = {}
        for r in results:
            if r[4] == "easy" and r[1] in A3_DEPTHS and r[3] == A3_BETA:
                c = r[5]
                ok.setdefault(r[0], []).append(bool(c[3]) and c[4] / c[3] >= A3_SHARE
                                                and fewer_errors(c))
        good = [t for t, v in ok.items() if all(v) and len(v) == 4]
        print(f"A3 (easy, beta {A3_BETA:g}, depths {A3_DEPTHS}): {len(good)} of {len(ok)} genomes "
              f"succeed in all 4 cells; failing: {' '.join(sorted(set(ok) - set(good))) or 'none'}")
        verdicts.append("A3 " + ("supported" if ok and len(good) >= A3_GENOMES * len(ok)
                                 else "not supported"))
        for phaser in cs.PHASERS:
            share = {t: c[4] / max(1, c[3]) for t, d, ph, b, g, c in results
                     if g == "easy" and d == "30" and ph == phaser and b == A3_BETA}
            afr = [v for t, v in share.items() if tests[t] == "AFR"]
            other = [v for t, v in share.items() if tests[t] in ("AMR", "EAS")]
            diff, p = welch_tost(afr, other, A4_MARGIN)
            print(f"A4 ({phaser}, easy, 30x, beta {A3_BETA:g}): AFR {len(afr)} genomes, AMR and "
                  f"EAS {len(other)}; difference {diff:.4f}, TOST p {p:.3g}")
            verdicts.append(f"A4 {phaser} " + ("equivalent" if p < A4_ALPHA else "not shown"))
    else:
        verdicts.append("A3 and A4 not scored: (easy, 1e-3) is not testable")
    print("\nA5 (descriptive), beta 1e-3, 30x: per genome and group, certified share, error "
          "among trusted pairs, truth pairs; then the share of pairs in groups certifying nothing")
    for t in tests:
        for phaser in cs.PHASERS:
            mine = [(g, c) for tt, d, ph, b, g, c in results
                    if tt == t and d == "30" and ph == phaser and b == 1e-3]
            text = "; ".join(f"{g} {c[4] / max(1, c[3]):.4f} {c[6] / max(1, c[4]):.1e} {c[3]}"
                             for g, c in sorted(mine))
            none = sum(c[3] for _, c in mine if not c[7]) / max(1, sum(c[3] for _, c in mine))
            print(f"  {t} {tests[t]} {phaser}: {text}; abstaining {none:.3f}")
    from genome_run import LENGTH_MB

    density = {}
    for t, sp in tests.items():
        n = 0
        for c in test_chroms:
            with gzip.open(common.ROOT / "data" / f"ext_{t}" / c / "hetsnv.vcf.gz", "rt") as h:
                n += sum(1 for line in h if not line.startswith("#"))
        density.setdefault(sp, []).append(n / sum(LENGTH_MB[c] for c in test_chroms))
    print("  heterozygous variants in the phaser input per Mb of chromosome, mean by "
          "superpopulation: " + "; ".join(f"{sp} {sum(v) / len(v):.0f} ({len(v)} genomes)"
                                           for sp, v in density.items()))
    print("\nA6 (descriptive): thresholds from HG002 against its dipcall truth instead of v5.0q, "
          "testable cells: trusted share and violations at the A1 and A2 levels")
    for name, chosen, registered in (("easy", "easy", A1_CELLS), ("difficult and indel", None,
                                                                    A2_CELLS)):
        for label, source in (("v5.0q", results), ("dipcall", a6)):
            mine = [r for r in source if is_testable(r)
                    and (r[4] == chosen if chosen else r[4] in DIFFICULT)]
            bad = sum(1 for r in mine if r[5][0] and r[5][2] < 0.05 / registered)
            share = sum(r[5][4] for r in mine) / max(1, sum(r[5][3] for r in mine))
            print(f"  {name:<20} {label:<8} cells {len(mine):>5} trusted {share:.4f} violating {bad}")
    print("\n" + "; ".join(verdicts))


if __name__ == "__main__":
    main()
