"""Measure the per-allele error `e` of the vote model on one chromosome.

As `docs/design/ont-prediction.md` defines it: each read with at least 3
called alleles at truth-phased sites is assigned to the haplotype it agrees
with at most of them (ties are skipped), and `e` is the share of those
alleles that disagree with the assigned haplotype. The truth is v5.0q for
HG002 and the platform's primary truth otherwise (transmission for HG005).
Chromosome 20 is a smoke-test chromosome on every platform and enters no test.

Usage: `PHASEBUDGET_PLATFORM=hg005 PYTHONPATH=src:benchmarks .venv/bin/python
benchmarks/measure_allele_error.py chr20`
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402


def allele_error(fragments, truth):
    """`(errors, alleles, reads, per-site error counts)` over assigned reads."""
    errors = alleles = reads = 0
    site_wrong, site_seen = {}, {}
    for f in fragments:
        calls = [(f.start + k, a) for k, a in enumerate(f.alleles)
                 if a is not None and truth[f.start + k] is not None]
        if len(calls) < 3:
            continue
        agree = sum(a == truth[j] for j, a in calls)
        if 2 * agree == len(calls):
            continue
        hap0 = 2 * agree > len(calls)
        reads += 1
        for j, a in calls:
            wrong = (a != truth[j]) if hap0 else (a == truth[j])
            errors += wrong
            alleles += 1
            site_seen[j] = site_seen.get(j, 0) + 1
            site_wrong[j] = site_wrong.get(j, 0) + wrong
    return errors, alleles, reads, site_wrong, site_seen


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chrom")
    parser.add_argument("--depth", default="full")
    args = parser.parse_args()
    extractor = common._module("extract_fragments.py")
    score = common._module("score_hg002.py")
    matrix, positions, _, _ = extractor.load_fragments(
        str(common.fragments_path(args.chrom, args.depth, "realign")))
    variants = score.variant_table(common.hetsnv_path(args.chrom))
    assert [p - 1 for p, _, _ in variants] == list(positions)
    if common.SAMPLE == "HG002":
        truth, name = common.v5_truth(args.chrom, variants), "GIAB v5.0q"
    else:
        truth, name = common.load_truths(args.chrom, positions)[1], "transmission"
    errors, alleles, reads, wrong, seen = allele_error(matrix.fragments, truth)
    noisy = [j for j in seen if seen[j] >= 5 and wrong[j] / seen[j] > 0.2]
    print(f"platform {common.PLATFORM}, sample {common.SAMPLE}, {args.chrom}, depth "
          f"{args.depth}, truth {name}: {reads} reads, {alleles} alleles, {errors} wrong, "
          f"e = {errors / alleles:.5f}; sites {len(seen)}, with error above 0.2 (5 or more "
          f"reads): {len(noisy)}, holding {sum(wrong[j] for j in noisy)} of the errors")


if __name__ == "__main__":
    main()
