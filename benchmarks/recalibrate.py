"""Fit the base-quality recalibration table registered in `genome-run-plan.md`.

On the full-depth `basequality` fragments of the odd confirmatory
chromosomes, each read is assigned to the truth haplotype it agrees with at
the majority of its truth-phased sites; reads with fewer than 3 such sites,
or a tie, are skipped. Its alleles at truth-phased sites are counted as right
or wrong per base quality value. Values are pooled from the top down until
every bin holds at least 50 wrong alleles; the last bin absorbs any
remainder. A bin's error rate is (wrong + 1) / (alleles + 2), and its weight
is the log-odds of that rate. `apply_table` maps a matrix's base qualities to
those weights.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

from phasebudget._types import Fragment
from phasebudget.fragments import FragmentMatrix

MIN_SITES = 3
MIN_WRONG = 50


def count_errors(matrix: FragmentMatrix, truth) -> tuple[Counter, Counter]:
    """`(alleles, wrong)` per integer base quality, for reads with a clear haplotype."""
    alleles: Counter = Counter()
    wrong: Counter = Counter()
    for fragment in matrix:
        calls = [(j, a, int(round(w))) for j, a, w in fragment.observed() if truth[j] is not None]
        if len(calls) < MIN_SITES:
            continue
        agree = sum(1 for j, a, _ in calls if a == truth[j])
        disagree = len(calls) - agree
        if agree == disagree:
            continue
        side = 0 if agree > disagree else 1
        for j, a, q in calls:
            alleles[q] += 1
            if a != truth[j] ^ side:
                wrong[q] += 1
    return alleles, wrong


def fit_table(alleles: Counter, wrong: Counter) -> list[dict]:
    """Bins from the highest quality down, each with at least `MIN_WRONG` wrong alleles."""
    bins: list[dict] = []
    current = {"high": None, "low": None, "alleles": 0, "wrong": 0}
    for q in sorted(alleles, reverse=True):
        if current["high"] is None:
            current["high"] = q
        current["low"] = q
        current["alleles"] += alleles[q]
        current["wrong"] += wrong[q]
        if current["wrong"] >= MIN_WRONG:
            bins.append(current)
            current = {"high": None, "low": None, "alleles": 0, "wrong": 0}
    if current["high"] is not None:
        if bins:
            last = bins[-1]
            last["low"] = current["low"]
            last["alleles"] += current["alleles"]
            last["wrong"] += current["wrong"]
        else:
            bins.append(current)
    for b in bins:
        rate = (b["wrong"] + 1) / (b["alleles"] + 2)
        b["error"] = rate
        b["weight"] = math.log((1 - rate) / rate)
    return bins


def weight_of(table: list[dict], quality: float) -> float:
    q = int(round(quality))
    for b in table:
        if b["low"] <= q <= b["high"]:
            return b["weight"]
    # Outside the fitted range: the nearest bin.
    nearest = min(table, key=lambda b: min(abs(q - b["low"]), abs(q - b["high"])))
    return nearest["weight"]


def apply_table(matrix: FragmentMatrix, table: list[dict]) -> FragmentMatrix:
    fragments = [
        Fragment(f.read_id, f.start, f.alleles,
                 tuple(0.0 if a is None else weight_of(table, w)
                       for a, w in zip(f.alleles, f.weights)))
        for f in matrix
    ]
    return FragmentMatrix(fragments, matrix.n_variants)


def main() -> None:
    import genome_common as common

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(common.DATA / "recalibration.json"))
    args = parser.parse_args()
    extractor = common._module("extract_fragments.py")
    alleles: Counter = Counter()
    wrong: Counter = Counter()
    used = []
    for chrom in common.ODD:
        path = common.fragments_path(chrom, "full", "basequality")
        matrix, positions, _, _ = extractor.load_fragments(str(path))
        filtered, _, _ = common.load_truths(chrom, positions)
        a, w = count_errors(matrix, filtered)
        alleles.update(a)
        wrong.update(w)
        used.append(chrom)
    table = fit_table(alleles, wrong)
    Path(args.out).write_text(json.dumps({"chromosomes": used, "bins": table}, indent=1))
    print(f"fitted on {', '.join(used)}: {sum(alleles.values())} alleles, "
          f"{sum(wrong.values())} wrong")
    for b in table:
        print(f"quality {b['low']:>3}-{b['high']:<3} alleles {b['alleles']:>10} wrong "
              f"{b['wrong']:>7} error {b['error']:.6f} weight {b['weight']:.3f}")


if __name__ == "__main__":
    main()
