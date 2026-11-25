"""Predict WhatsHap's selection minus the exact one from a vote model.

`docs/design/ont-prediction.md` fixes the model and the rule. For each
chromosome, depth and cap, both selections are recomputed from the
`realign` fragments, every boundary gets `i`, the number of selected reads
calling both of its variants, and the predicted difference is the mean of
`f(i_whatshap) - f(i_optimal)` over boundaries both selections call, per
10,000. No truth is read.

Usage: `PHASEBUDGET_PLATFORM=ont PYTHONPATH=src:benchmarks .venv/bin/python
benchmarks/predict_regime.py --error 0.00503`
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

from phasebudget.selection import read_weights, select_optimal  # noqa: E402

pipelines = common._module("hg002_pipelines.py")
extractor = common._module("extract_fragments.py")
CAPS = (3, 5, 10, 15)


def vote_error(i: int, q: float) -> float:
    """P(a majority of `i` independent readings is wrong), ties split evenly."""
    total = 0.0
    for k in range(i + 1):
        p = math.comb(i, k) * q ** k * (1 - q) ** (i - k)
        if 2 * k > i:
            total += p
        elif 2 * k == i:
            total += p / 2
    return total


def informative(matrix, selected):
    counts = [0] * (matrix.n_variants - 1)
    for f in matrix.subset(selected):
        for b in range(f.start, f.end - 1):
            if f.alleles[b - f.start] is not None and f.alleles[b + 1 - f.start] is not None:
                counts[b] += 1
    return counts


def joined(chrom, depth, name, cap):
    """Per boundary, whether the stored phasing puts both variants in one block."""
    import gzip
    import json

    path = common.chrom_dir(chrom) / "runs" / f"lib_{depth}_{name}_c{cap}.json.gz"
    with gzip.open(path, "rt") as handle:
        record = json.load(handle)
    a, b = record["alleles"], record["blocks"]
    return [a[j] >= 0 and a[j + 1] >= 0 and b[j] >= 0 and b[j] == b[j + 1]
            for j in range(len(a) - 1)]


def one(job):
    chrom, depth, cap, q, g0 = job
    matrix, positions, _, _ = extractor.load_fragments(
        str(common.fragments_path(chrom, depth, "realign")))
    matrix = extractor.phred_to_log_odds(matrix)
    w = informative(matrix, pipelines.whatshap_selection(matrix, cap, positions).selected)
    o = informative(matrix, select_optimal(matrix, cap, read_weights(matrix)).selected)
    f = [vote_error(i, q) for i in range(cap + 1)]
    if g0 is None:
        both = [(a, b) for a, b in zip(w, o) if a >= 1 and b >= 1]
    else:
        # Amendment 1: every boundary both phasings join; gap-linked ones (i = 0) at g0.
        f[0] = g0
        jw, jo = joined(chrom, depth, "whatshap", cap), joined(chrom, depth, "optimal", cap)
        both = [(a, b) for a, b, x, y in zip(w, o, jw, jo) if x and y]
    return chrom, depth, cap, len(both), sum(f[a] - f[b] for a, b in both)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--error", type=float, required=True, help="per-allele error e")
    parser.add_argument("--chromosomes", nargs="+", default=list(common.CONFIRMATORY))
    parser.add_argument("--caps", nargs="+", type=int, default=list(CAPS))
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--tsv", help="also write chrom, depth, cap, boundaries, summed difference")
    parser.add_argument("--gap-rate", type=float,
                        help="amendment 1: score every boundary both phasings join, with this "
                             "error rate where no selected read calls both variants")
    args = parser.parse_args()
    q = 2 * args.error * (1 - args.error)
    jobs = [(c, d, cap, q, args.gap_rate) for c in args.chromosomes for d in common.DEPTHS for cap in args.caps]
    pooled: dict = {}
    rows = []
    with ProcessPoolExecutor(args.jobs) as pool:
        for chrom, depth, cap, n, diff in pool.map(one, jobs):
            rows.append((chrom, depth, cap, n, diff))
            acc = pooled.setdefault((depth, cap), [0, 0.0])
            acc[0] += n
            acc[1] += diff
    if args.tsv:
        with open(args.tsv, "w") as handle:
            handle.write("chrom\tdepth\tcap\tboundaries\tsummed_difference\n")
            for row in rows:
                handle.write("\t".join(map(str, row)) + "\n")
    model = "registered" if args.gap_rate is None else f"amended, g0 = {args.gap_rate}"
    print(f"vote-model prediction ({model}), platform {common.PLATFORM}, e = {args.error}, "
          f"q = {q:.6f}, "
          f"chromosomes {len(args.chromosomes)}; whatshap minus optimal per 10,000 boundaries "
          "both selections call; no truth read")
    for cap in args.caps:
        print(f"  cap {cap:>2}: " + "; ".join(
            f"{d}x {1e4 * pooled[(d, cap)][1] / max(pooled[(d, cap)][0], 1):+.2f}"
            for d in common.DEPTHS))


if __name__ == "__main__":
    main()
