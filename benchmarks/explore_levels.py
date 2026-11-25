"""Exploratory: can a count model predict each selection's own error rate?

The registered vote model predicts the difference between two selections
well and each selection's own rate badly (1.3 to 4 times too high on ONT).
For every boundary a phasing joins and the truth phases, this compares the
observed wrong calls with two truth-free predictions:
- `f(i)`, with `i` the selected reads calling both variants (the registered
  count);
- `f(c)`, with `c` the selected reads having calls on both sides of the
  boundary, which the DP also uses.

This reads truths of seen chromosomes; it supports no claim. Usage:
`PHASEBUDGET_PLATFORM=ont PYTHONPATH=src:benchmarks .venv/bin/python
benchmarks/explore_levels.py --error 0.00503 --chromosomes chr21 chr22`
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402
from predict_regime import vote_error  # noqa: E402

from phasebudget.selection import read_weights, select_optimal  # noqa: E402

pipelines = common._module("hg002_pipelines.py")

extractor = common._module("extract_fragments.py")
score = common._module("score_hg002.py")


def counts(fragments, selected, n):
    """`(i, c)` per boundary: reads calling both variants, reads with calls on both sides."""
    i = [0] * (n - 1)
    diff = [0] * n
    for f in (fragments[k] for k in selected):
        called = [f.start + k for k, a in enumerate(f.alleles) if a is not None]
        if len(called) < 2:
            continue
        diff[called[0]] += 1
        diff[called[-1]] -= 1
        for a, b in zip(called, called[1:]):
            if b == a + 1:
                i[a] += 1
    c, run = [0] * (n - 1), 0
    for b in range(n - 1):
        run += diff[b]
        c[b] = run
    return i, c


def one(job):
    chrom, depth, cap, q = job
    matrix, positions, _, _ = extractor.load_fragments(
        str(common.fragments_path(chrom, depth, "realign")))
    variants = score.variant_table(common.hetsnv_path(chrom))
    if common.SAMPLE == "HG002":
        truth = common.v5_truth(chrom, variants)
    else:
        truth = common.load_truths(chrom, positions)[1]
    matrix = extractor.phred_to_log_odds(matrix)
    selections = {
        "whatshap": pipelines.whatshap_selection(matrix, cap, positions).selected,
        "optimal": select_optimal(matrix, cap, read_weights(matrix)).selected,
    }
    f = [vote_error(k, q) for k in range(cap + 1)]
    out = []
    for arm in ("whatshap", "optimal"):
        path = common.chrom_dir(chrom) / "runs" / f"lib_{depth}_{arm}_c{cap}.json.gz"
        with gzip.open(path, "rt") as handle:
            r = json.load(handle)
        a, blocks = r["alleles"], r["blocks"]
        i, c = counts(matrix.fragments, selections[arm], matrix.n_variants)
        n = wrong = 0
        pi = pc = 0.0
        for b in range(matrix.n_variants - 1):
            if not (a[b] >= 0 and a[b + 1] >= 0 and blocks[b] >= 0 and blocks[b] == blocks[b + 1]):
                continue
            if truth[b] is None or truth[b + 1] is None or i[b] == 0:
                continue
            n += 1
            wrong += (a[b] == a[b + 1]) != (truth[b] == truth[b + 1])
            pi += f[min(i[b], cap)]
            pc += f[min(c[b], cap)]
        out.append((depth, cap, arm, n, wrong, pi, pc))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--error", type=float, required=True)
    parser.add_argument("--chromosomes", nargs="+", required=True)
    parser.add_argument("--caps", nargs="+", type=int, default=[3, 5, 15])
    parser.add_argument("--jobs", type=int, default=3)
    args = parser.parse_args()
    q = 2 * args.error * (1 - args.error)
    jobs = [(ch, d, cap, q) for ch in args.chromosomes for d in common.DEPTHS for cap in args.caps]
    pooled = {}
    with ProcessPoolExecutor(args.jobs) as pool:
        for rows in pool.map(one, jobs):
            for depth, cap, arm, n, wrong, pi, pc in rows:
                acc = pooled.setdefault((depth, cap, arm), [0, 0, 0.0, 0.0])
                for k, v in enumerate((n, wrong, pi, pc)):
                    acc[k] += v
    print(f"exploratory levels, platform {common.PLATFORM}, sample {common.SAMPLE}, e = "
          f"{args.error}, chromosomes {' '.join(args.chromosomes)}; per 10,000 boundaries "
          "the phasing joins, the truth phases and i >= 1")
    print(f"{'depth':>5} {'cap':>3} {'arm':>8} {'boundaries':>10} {'observed':>9} "
          f"{'f(i)':>7} {'f(c)':>7}")
    for (depth, cap, arm), (n, wrong, pi, pc) in sorted(pooled.items(),
                                                        key=lambda kv: (kv[0][0], kv[0][1])):
        print(f"{depth:>5} {cap:>3} {arm:>8} {n:>10} {1e4 * wrong / n:>9.2f} "
              f"{1e4 * pi / n:>7.2f} {1e4 * pc / n:>7.2f}")


if __name__ == "__main__":
    main()
