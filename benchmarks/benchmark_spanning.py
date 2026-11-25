"""Why selection matters only at small caps: errors by spanning-read count.

Exploratory, on the genome-wide HiFi run, after its results were seen. For
each chromosome, depth and cap, the selections of `whatshap` and `optimal`
are recomputed (both are deterministic, so they are the selections behind
the stored phasings) and every boundary between index-adjacent variants is
given `s`, the number of selected reads spanning it. A boundary is scored
when the truth phases both of its variants and the phasing puts them in one
block; it is wrong when the phasing's relation differs from the truth's.

The table gives, per selector, the number of scored boundaries and the
error rate at each `s`, and the share of boundaries at each `s`.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/benchmark_spanning.py`
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

from phasebudget.selection import read_weights, select_optimal  # noqa: E402

pipelines = common._module("hg002_pipelines.py")
extractor = common._module("extract_fragments.py")
score = common._module("score_hg002.py")


def one(job):
    chrom, depth, cap, truth_name, count = job
    matrix, positions, _, _ = extractor.load_fragments(
        str(common.fragments_path(chrom, depth, "realign")))
    matrix = extractor.phred_to_log_odds(matrix)
    variants = score.variant_table(common.hetsnv_path(chrom))
    if truth_name == "v5":
        truth = common.v5_truth(chrom, variants)
    else:
        truth = common.load_truths(chrom, positions)[0]
    selections = {
        "whatshap": pipelines.whatshap_selection(matrix, cap, positions).selected,
        "optimal": select_optimal(matrix, cap, read_weights(matrix)).selected,
    }
    out = {}
    for name, selected in selections.items():
        spans = [0] * (matrix.n_variants - 1)
        for fragment in matrix.subset(selected):
            for b in range(fragment.start, fragment.end - 1):
                if count == "span" or (fragment.alleles[b - fragment.start] is not None
                                       and fragment.alleles[b + 1 - fragment.start] is not None):
                    spans[b] += 1
        path = common.chrom_dir(chrom) / "runs" / f"lib_{depth}_{name}_c{cap}.json.gz"
        with gzip.open(path, "rt") as handle:
            record = json.load(handle)
        alleles, blocks = record["alleles"], record["blocks"]
        table = {}
        for b, s in enumerate(spans):
            row = table.setdefault(s, [0, 0, 0])
            row[0] += 1
            if truth[b] is None or truth[b + 1] is None or alleles[b] < 0 or alleles[b + 1] < 0:
                continue
            if blocks[b] != blocks[b + 1]:
                continue
            row[1] += 1
            row[2] += (alleles[b] ^ alleles[b + 1]) != (truth[b] ^ truth[b + 1])
        out[name] = table
    return (depth, cap), out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", choices=("v5", "filtered"), default="v5")
    parser.add_argument("--depths", nargs="+", default=["full", "10"])
    parser.add_argument("--caps", nargs="+", type=int, default=[3, 5])
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--count", choices=("span", "informative"), default="span",
                        help="span: reads whose span covers the boundary; informative: reads "
                             "that call both of its variants")
    args = parser.parse_args()
    jobs = [(c, d, cap, args.truth, args.count) for c in common.CONFIRMATORY for d in args.depths
            for cap in args.caps]
    pooled: dict = {}
    with ProcessPoolExecutor(args.jobs) as pool:
        for key, out in pool.map(one, jobs):
            for name, table in out.items():
                target = pooled.setdefault((*key, name), {})
                for s, row in table.items():
                    acc = target.setdefault(s, [0, 0, 0])
                    for i in range(3):
                        acc[i] += row[i]
    print(f"errors by spanning-read count ({args.count}); platform {common.PLATFORM}, truth {args.truth}, "
          f"{len(common.CONFIRMATORY)} chromosomes; EXPLORATORY, after the results were seen")
    print(f"{'depth':>5} {'cap':>3} {'selector':>9} {'s':>3} {'boundaries':>11} {'share':>7} "
          f"{'scored':>9} {'wrong':>6} {'per 10k':>8}")
    for (depth, cap, name), table in sorted(pooled.items()):
        total = sum(row[0] for row in table.values())
        for s in sorted(table):
            n, scored, wrong = table[s]
            rate = 1e4 * wrong / scored if scored else 0.0
            print(f"{depth:>5} {cap:>3} {name:>9} {s:>3} {n:>11} {n / total:>7.4f} "
                  f"{scored:>9} {wrong:>6} {rate:>8.2f}")
        print()


if __name__ == "__main__":
    main()
