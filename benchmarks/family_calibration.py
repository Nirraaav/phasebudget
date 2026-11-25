"""Y1 step 2, `docs/design/family-and-clinical-gates.md`: calibrate on transmission labels.

The registered HG002 ONT sup calibration rows get three losses: against
v5.0q (`truth`), against transmission labels from curated genotypes
(`labels-truth`) and from Clair3 calls (`labels-called`), on the pairs both
label sources cover. Each calibration's thresholds, pooled and per stratum,
are scored on HG002 ONT hac (v5.0q) and HG001 ONT sup (Platinum Pedigree).

The data were seen before: this is development, not a test.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/family_calibration.py`
"""

from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import confidence_score as cs  # noqa: E402
import genome_common as common  # noqa: E402
from family_labels import transmission  # noqa: E402

from phasebudget.confidence import Units, binomial_cdf, learn_then_test_by_group  # noqa: E402

SOURCES = ("truth", "labels-truth", "labels-called")
LEVEL = 0.05


def labelled(chrom: str, min_gq: int = 0):
    """`{pos1: (curated label, called label, v5.0q first allele)}` where both sources label."""
    curated, _, child = transmission(chrom)
    called, _, called_child = transmission(chrom, called=True, min_gq=min_gq)
    return chrom, {p: (curated[p][0], called[p][0], child[p].genotype[0]) for p in curated
                   if p in called and child[p].phased and called_child[p].alt == child[p].alt}


def score(thresholds: dict, test: list, beta: float) -> tuple[int, float, int]:
    """`(groups violating at LEVEL, trusted share, wrong trusted pairs)` on `test` rows."""
    violations = trusted = wrong = 0
    for group in {r[0] for r in test}:
        t = thresholds.get(group)
        if t is None:
            continue
        rows = [r for r in test if r[0] == group]
        n, k = Units.from_pairs([(u, key, bad) for _, u, key, bad in rows]).at(t)
        p = 1.0 if k == 0 else 1.0 - binomial_cdf(k - 1, n, beta)
        violations += p < LEVEL
        kept = [r for r in rows if r[2] >= t]
        trusted += len(kept)
        wrong += sum(r[3] for r in kept)
    return violations, trusted / len(test), wrong


def main() -> None:
    import argparse
    from functools import partial

    from stratified_gate import hg001_pedigree, strata

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-gq", type=int, default=0, help="step 3's filter on called labels")
    parser.add_argument("--even", action="store_true",
                        help="step 3: calibrate on the even chromosomes only")
    args = parser.parse_args()
    chroms = list(common.CONFIRMATORY)
    cal_chroms = [c for c in chroms if not args.even or int(c[3:]) % 2 == 0]
    print(f"min GQ {args.min_gq}; calibration chromosomes {' '.join(cal_chroms)}")
    print("Y1 step 2, docs/design/family-and-clinical-gates.md (development; data seen): HG002 "
          "ONT sup calibration rows restricted to pairs both label sources cover; losses truth "
          "(v5.0q), labels-truth (transmission, curated HG002 v5.0q and HG003/HG004 v4.2.1 in "
          "their BEDs), labels-called (transmission, Clair3 PASS calls of giab_2023.05 "
          "*_sup_all, absent parent site = homozygous reference); tests HG002 ONT hac (v5.0q, all "
          "pairs) and HG001 ONT sup (Platinum Pedigree, track C's adjudicable pairs); pooled and "
          f"per-stratum Learn then Test; betas {cs.BETAS}; delta {cs.DELTA}; units "
          f"{cs.UNIT // 1000} kb; violation p < {LEVEL} per group; chromosomes {' '.join(chroms)}")
    with ProcessPoolExecutor(6) as pool:
        labels = dict(pool.map(partial(labelled, min_gq=args.min_gq), cal_chroms))
    print(f"SNVs labelled by both sources and phased by v5.0q: "
          f"{sum(len(v) for v in labels.values())}")
    stratum = strata(chroms)
    hg001 = hg001_pedigree(chroms)
    failures, shares = [], {}
    for depth in cs.DEPTHS:
        for phaser in cs.PHASERS:
            cal = {s: [] for s in SOURCES}
            differ = {"labels-truth": [0, 0], "labels-called": [0, 0]}
            for c, p1, key, wrong, p2, _ in cs.load_rows("ont", depth, phaser, chroms):
                if c not in labels:
                    continue
                a, b = labels[c].get(p1), labels[c].get(p2)
                if a is None or b is None:
                    continue
                same_phaser = (a[2] == b[2]) != wrong
                losses = {"truth": wrong, "labels-truth": same_phaser != (a[0] == b[0]),
                          "labels-called": same_phaser != (a[1] == b[1])}
                for s in SOURCES:
                    cal[s].append((c, p1, p2, key, losses[s]))
                for s in differ:
                    if losses[s] != wrong:
                        differ[s][0 if wrong else 1] += 1  # [errors hidden, errors invented]
            hac = [(c, p1, p2, key, wrong)
                   for c, p1, key, wrong, p2, _ in cs.load_rows("ont_hac", depth, phaser, chroms)]
            ped = hg001(depth, phaser)
            ped = [(c, p1, p2, key, bad) for c, p1, key, bad, p2 in ped]
            print(f"\n{depth}x {phaser}: calibration pairs {len(cal['truth'])}, truth-wrong "
                  f"{sum(r[4] for r in cal['truth'])}; label against truth: "
                  + "; ".join(f"{s} hides {h}, invents {i}" for s, (h, i) in differ.items()))
            print("  mode        beta    source          HG002 hac: viol  share   wrong"
                  "    HG001: viol  share   wrong")
            for mode in ("pooled", "stratified"):
                def grouped(rows):
                    return [("all" if mode == "pooled" else stratum(c, p1, p2),
                             (c, p1 // cs.UNIT), key, loss) for c, p1, p2, key, loss in rows]

                tests = (grouped(hac), grouped(ped))
                for beta in cs.BETAS:
                    results = {}
                    for s in SOURCES:
                        thresholds = learn_then_test_by_group(grouped(cal[s]), beta, cs.DELTA)
                        results[s] = [score(thresholds, t, beta) for t in tests]
                        print(f"  {mode:<10}  {beta:<6g}  {s:<14}  "
                              f"{results[s][0][0]:>13}  {results[s][0][1]:.4f}  "
                              f"{results[s][0][2]:>6}  {results[s][1][0]:>11}  "
                              f"{results[s][1][1]:.4f}  {results[s][1][2]:>6}")
                    for i, name in enumerate(("HG002 hac", "HG001")):
                        if results["labels-called"][i][0] > results["truth"][i][0]:
                            failures.append(f"{depth}x {phaser} {mode} {beta:g} {name}")
                    if beta == 1e-3 and depth == "30":
                        shares[(phaser, mode)] = (results["labels-called"][0][1],
                                                  results["truth"][0][1])
    print(f"\n(a) cells where labels-called violates in more groups than truth: {len(failures)}"
          + "".join(f"\n  {f}" for f in failures))
    close = all(abs(a - b) <= 0.01 for a, b in shares.values())
    print("(b) HG002 hac trusted share at beta 1e-3, 30x, labels-called against truth: "
          + ", ".join(f"{p} {m} {a:.4f} vs {b:.4f}" for (p, m), (a, b) in shares.items()))
    print(f"Y1 step 2: {'continue' if not failures and close else 'stop'}")


if __name__ == "__main__":
    main()
