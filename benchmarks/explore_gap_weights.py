"""Exploratory: weight each read's alleles by the error its gap fraction predicts.

WhatsHap's `realign` alleles all carry phred 30, but a read's share of
uncalled sites inside its span (its gap fraction) predicts its allele error:
on ONT chromosome 20, 0.0035 below 2 per cent gaps and 0.025 at 10 to 20 per
cent. This phases with constant weights and with gap weights, under
WhatsHap's selection and the exact one, and scores both against the truth.

Development set only: chromosome 20 fits the error table and is scored here.
Usage: `PHASEBUDGET_PLATFORM=ont PYTHONPATH=src:benchmarks .venv/bin/python
benchmarks/explore_gap_weights.py chr20 --caps 5 15`
"""

from __future__ import annotations

import argparse
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402

from phasebudget._types import Fragment  # noqa: E402
from phasebudget.dp import CoverageDP  # noqa: E402
from phasebudget.fragments import FragmentMatrix  # noqa: E402
from phasebudget.metrics import comparable_pairs, phasing_from_solution, switch_positions  # noqa: E402
from phasebudget.paired import common_pairs  # noqa: E402
from phasebudget.selection import read_weights, select_optimal  # noqa: E402

pipelines = common._module("hg002_pipelines.py")
extractor = common._module("extract_fragments.py")
score = common._module("score_hg002.py")

# Allele error by gap fraction, ONT chromosome 20 against v5.0q
# (`benchmarks/results/runs/2026-10-01_gap_error_ont_chr20.txt`).
GAP_ERROR = ((0.02, 0.00346), (0.05, 0.00557), (0.10, 0.01197), (0.20, 0.02455),
             (0.40, 0.03578), (1.01, 0.13889))


def gap_fraction(f) -> float:
    return 1 - sum(a is not None for a in f.alleles) / len(f.alleles)


def gap_weighted(matrix: FragmentMatrix) -> FragmentMatrix:
    """The same reads, each called allele weighted `log((1 - e) / e)` for its read's gap `e`."""
    out = []
    for f in matrix.fragments:
        g = gap_fraction(f)
        e = next(err for upper, err in GAP_ERROR if g < upper)
        w = math.log((1 - e) / e)
        out.append(Fragment(f.read_id, f.start, f.alleles,
                            tuple(w if a is not None else 0.0 for a in f.alleles)))
    return FragmentMatrix(out, matrix.n_variants)


def one(job):
    chrom, cap, select_name, weight_name = job
    matrix, positions, _, _ = extractor.load_fragments(
        str(common.fragments_path(chrom, "full", "realign")))
    matrix = extractor.phred_to_log_odds(matrix)
    weighted = gap_weighted(matrix) if weight_name == "gap" else matrix
    if select_name == "whatshap":
        selected = pipelines.whatshap_selection(matrix, cap, positions).selected
    else:
        selected = select_optimal(weighted, cap, read_weights(weighted)).selected
    subset = weighted.subset(selected)
    phasing = phasing_from_solution(subset, CoverageDP(subset).solve())
    return job, phasing


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chrom")
    parser.add_argument("--caps", nargs="+", type=int, default=[5, 15])
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    variants = score.variant_table(common.hetsnv_path(args.chrom))
    positions = [p - 1 for p, _, _ in variants]
    truth = common.v5_truth(args.chrom, variants)
    jobs = [(args.chrom, cap, s, w) for cap in args.caps for s in ("whatshap", "optimal")
            for w in ("constant", "gap")]
    with ProcessPoolExecutor(args.jobs) as pool:
        phasings = dict(pool.map(one, jobs))
    print(f"exploratory gap weights, platform {common.PLATFORM}, {args.chrom}, full depth, "
          "truth GIAB v5.0q")
    for cap in args.caps:
        for s in ("whatshap", "optimal"):
            a, b = phasings[(args.chrom, cap, s, "constant")], phasings[(args.chrom, cap, s, "gap")]
            line = []
            for name, p in (("constant", a), ("gap", b)):
                sw, pairs = len(switch_positions(truth, p)), comparable_pairs(p, truth)
                line.append(f"{name} {sw}/{pairs} = {1e4 * sw / pairs:.2f}")
            rows = common_pairs(truth, a, b)
            only_a = sum(1 for _, x, y in rows if x and not y)
            only_b = sum(1 for _, x, y in rows if y and not x)
            print(f"cap {cap:>2} {s:>8}: {'; '.join(line)} per 10k; common pairs {len(rows)}, "
                  f"only constant wrong {only_a}, only gap wrong {only_b}")


if __name__ == "__main__":
    main()
