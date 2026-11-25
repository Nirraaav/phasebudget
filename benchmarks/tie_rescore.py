"""Post hoc: re-score library phasings with tied sites left unphased, as WhatsHap does.

At a site where the selected reads, each oriented to fit its block, support
both alleles with equal weight, flipping the site's allele leaves the wMEC
cost unchanged. An exact solver must guess there. `whatshap phase` leaves
such a site unphased; the library phases it. On ONT chromosome 20 at 10-fold
and cap 15, 21 such sites held 17 of the library's 31 switch errors.

For each chromosome, depth and cap this recomputes both selections, marks
tied sites in each stored phasing, removes them, and compares `whatshap`
with `optimal` on common pairs against the truth: v5.0q for HG002,
transmission for HG005. It prints the pooled difference per 10,000 common
pairs with its 1 Mb window-bootstrap interval, the sign test, and the tied
sites of each selection. Usage: `PHASEBUDGET_PLATFORM=ont
PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/tie_rescore.py`
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

from phasebudget._types import Phasing  # noqa: E402
from phasebudget.paired import common_pairs, sign_test, window_discordance  # noqa: E402
from phasebudget.selection import read_weights, select_optimal  # noqa: E402

pipelines = common._module("hg002_pipelines.py")
extractor = common._module("extract_fragments.py")
score = common._module("score_hg002.py")


def tied_sites(fragments, alleles, blocks):
    """Sites where the reads, each oriented to fit its block, weigh equally for both alleles."""
    support = {}
    for f in fragments:
        groups = {}
        for k, allele in enumerate(f.alleles):
            j = f.start + k
            if allele is None or alleles[j] < 0:
                continue
            groups.setdefault(blocks[j], []).append((j, allele, f.weights[k]))
        for calls in groups.values():
            agree = sum(w for j, a, w in calls if a == alleles[j])
            flip = agree < sum(w for _, _, w in calls) - agree
            for j, a, w in calls:
                s = support.setdefault(j, [0.0, 0.0])
                s[(a == alleles[j]) == flip] += w
    return {j for j, (pro, con) in support.items() if abs(pro - con) < 1e-9}


def one(job):
    chrom, depth, cap = job
    matrix, positions, _, _ = extractor.load_fragments(
        str(common.fragments_path(chrom, depth, "realign")))
    matrix = extractor.phred_to_log_odds(matrix)
    variants = score.variant_table(common.hetsnv_path(chrom))
    truth = (common.v5_truth(chrom, variants) if common.SAMPLE == "HG002"
             else common.load_truths(chrom, positions)[1])
    selections = {"whatshap": pipelines.whatshap_selection(matrix, cap, positions).selected,
                  "optimal": select_optimal(matrix, cap, read_weights(matrix)).selected}
    phasings, ties = {}, {}
    for arm, selected in selections.items():
        path = common.chrom_dir(chrom) / "runs" / f"lib_{depth}_{arm}_c{cap}.json.gz"
        with gzip.open(path, "rt") as handle:
            record = json.load(handle)
        a, b = record["alleles"], record["blocks"]
        tied = tied_sites([matrix.fragments[k] for k in selected], a, b)
        ties[arm] = len(tied)
        phasings[arm] = Phasing(tuple(None if x < 0 or j in tied else x for j, x in enumerate(a)),
                                tuple(None if x < 0 or j in tied else y
                                      for j, (x, y) in enumerate(zip(a, b))))
    rows = common_pairs(truth, phasings["whatshap"], phasings["optimal"])
    return (depth, cap, len(rows), sum(1 for _, x, y in rows if x and not y),
            sum(1 for _, x, y in rows if y and not x), sum(1 for _, x, y in rows if x and y),
            window_discordance(rows, positions, 1_000_000), ties["whatshap"], ties["optimal"])


def main() -> None:
    from phasebudget.bootstrap import statistic_interval

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chromosomes", nargs="+", default=list(common.CONFIRMATORY))
    parser.add_argument("--depths", nargs="+", default=list(common.DEPTHS))
    parser.add_argument("--caps", nargs="+", type=int, default=[3, 5, 15])
    parser.add_argument("--jobs", type=int, default=3)
    args = parser.parse_args()
    jobs = [(c, d, cap) for c in args.chromosomes for d in args.depths for cap in args.caps]
    pooled = {}
    with ProcessPoolExecutor(args.jobs) as pool:
        for depth, cap, n, first, second, both, windows, tw, to in pool.map(one, jobs):
            acc = pooled.setdefault((cap, depth), [0, 0, 0, 0, [], 0, 0])
            for k, v in enumerate((n, first, second, both)):
                acc[k] += v
            acc[4].extend(windows)
            acc[5] += tw
            acc[6] += to
    truth = "GIAB v5.0q" if common.SAMPLE == "HG002" else "transmission"
    print(f"tie-aware re-score (post hoc), platform {common.PLATFORM}, sample {common.SAMPLE}, "
          f"truth {truth}, {len(args.chromosomes)} chromosomes; tied sites unphased in both "
          "selections; whatshap minus optimal per 10,000 common pairs")
    print(f"{'cap':>3} {'depth':>5} {'pairs':>9} {'only W':>7} {'only O':>7} {'both':>5} "
          f"{'p':>9} {'diff [95%]':>22} {'tied W':>7} {'tied O':>7}")
    for (cap, depth) in sorted(pooled, key=lambda k: (k[0], args.depths.index(k[1]))):
        n, first, second, both, windows, tw, to = pooled[(cap, depth)]

        def stat(picks, windows=windows):
            m = sum(windows[k][0] for k in picks)
            return None if m == 0 else 1e4 * (sum(windows[k][1] for k in picks)
                                              - sum(windows[k][2] for k in picks)) / m

        iv = statistic_interval(stat, len(windows), resamples=2000, seed=0)
        print(f"{cap:>3} {depth:>5} {n:>9} {first:>7} {second:>7} {both:>5} "
              f"{sign_test(first, first + second):>9.2e} {iv.estimate:+7.2f} "
              f"[{iv.low:+.2f}, {iv.high:+.2f}] {tw:>7} {to:>7}")


if __name__ == "__main__":
    main()
