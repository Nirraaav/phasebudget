"""The objective-truth gap and margin informativeness along the density axis.

`benchmark_objective_truth.py` varies coverage, read span, and error rate, and
records as a stated limit that it cannot vary heterozygosity: its reads are
measured in variants, so density does not exist in it. This script uses
`simulate_genomic`, where reads have a physical length and variants are spaced
by density, so the quantity being swept is phase linkage:

    E[variants per read] = heterozygosity * read_length.

Connectivity is an outcome here, not a control. Lowering the density is exactly
what breaks a phase block, so forcing one block per instance would discard the
main effect and keep an unrepresentative remainder. Blocks, phased fraction,
and comparable pairs are therefore reported as columns, and the switch error
rate is a rate over whatever pairs remained comparable. That makes the switch
column across density rows a comparison over different denominators, which is
the honest form of this measurement rather than a defect of it.
"""

from __future__ import annotations

import argparse
import statistics

from phasebudget._validation import CoverageLimitExceeded
from phasebudget.dp import CoverageDP
from phasebudget.metrics import (
    comparable_pairs,
    phased_fraction,
    phasing_from_solution,
    switch_positions,
)
from phasebudget.posterior import posterior_margins
from phasebudget.simulate import simulate_genomic


def rank_probability(bad: list[float], good: list[float]) -> float:
    """P(margin at a switch < margin elsewhere), ties counted as one half."""
    if not bad or not good:
        return float("nan")
    wins = 0.0
    for left in bad:
        for right in good:
            if left < right:
                wins += 1.0
            elif left == right:
                wins += 0.5
    return wins / (len(bad) * len(good))


def _mean(values: list[float]) -> float:
    return statistics.fmean(values) if values else float("nan")


def run(
    heterozygosities,
    error_rates,
    replicates,
    n_variants,
    coverage,
    read_length,
    max_active_reads,
    include_unspanned=False,
):
    """One row per (heterozygosity, error rate). Instances over the limit skip."""
    rows = []
    for heterozygosity in heterozygosities:
        for error_rate in error_rates:
            per_read: list[float] = []
            depth: list[float] = []
            blocks: list[float] = []
            phased: list[float] = []
            pairs: list[float] = []
            gaps: list[float] = []
            switch_rates: list[float] = []
            probabilities: list[float] = []
            posterior_probabilities: list[float] = []
            truth_optimal = 0
            used = 0
            skipped = 0
            for replicate in range(replicates):
                instance = simulate_genomic(
                    n_variants,
                    coverage=coverage,
                    read_length=read_length,
                    heterozygosity=heterozygosity,
                    error_rate=error_rate,
                    seed=2100 + 19 * replicate,
                )
                matrix = instance.matrix
                if not len(matrix):
                    skipped += 1
                    continue
                try:
                    dp = CoverageDP(matrix, max_active_reads=max_active_reads)
                except CoverageLimitExceeded:
                    skipped += 1
                    continue
                used += 1
                solution = dp.solve()
                truth_cost = matrix.cost(instance.haplotype)
                gaps.append(truth_cost - solution.cost)
                truth_optimal += truth_cost <= solution.cost + 1e-9

                per_read.append(statistics.fmean(len(f) for f in matrix))
                depth.append(statistics.fmean(matrix.coverage()))

                phasing = phasing_from_solution(matrix, solution)
                phased.append(phased_fraction(phasing))
                comparable = comparable_pairs(phasing)
                pairs.append(float(comparable))
                blocks.append(
                    float(len({b for b in phasing.blocks if b is not None}))
                )
                switches = set(switch_positions(instance.haplotype, phasing))
                switch_rates.append(
                    len(switches) / comparable if comparable else 0.0
                )
                # A boundary no read spans has margin zero by construction and
                # can never carry a switch, since switches are counted inside
                # blocks. Ranking it among the "no switch" boundaries charges
                # the signal for a block break, which was the error in the 10
                # Sep table; only spanned boundaries are ranked unless
                # `include_unspanned` reproduces that table.
                margins = [m for m in dp.margins() if m.spanned or include_unspanned]
                bad = [m.margin for m in margins if m.boundary + 1 in switches]
                good = [m.margin for m in margins if m.boundary + 1 not in switches]
                probability = rank_probability(bad, good)
                if probability == probability:
                    probabilities.append(probability)
                soft = [m for m in posterior_margins(matrix) if m.spanned or include_unspanned]
                bad = [m.margin for m in soft if m.boundary + 1 in switches]
                good = [m.margin for m in soft if m.boundary + 1 not in switches]
                probability = rank_probability(bad, good)
                if probability == probability:
                    posterior_probabilities.append(probability)
            rows.append(
                {
                    "heterozygosity": heterozygosity,
                    "error_rate": error_rate,
                    "per_read": _mean(per_read),
                    "depth": _mean(depth),
                    "blocks": _mean(blocks),
                    "phased": _mean(phased),
                    "pairs": _mean(pairs),
                    "gap": _mean(gaps),
                    "truth_optimal": truth_optimal / used if used else float("nan"),
                    "switch": _mean(switch_rates),
                    "informative": _mean(probabilities),
                    "posterior_informative": _mean(posterior_probabilities),
                    "instances": used,
                    "skipped": skipped,
                    "with_a_switch": len(probabilities),
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=int, default=40)
    parser.add_argument("--coverage", type=float, default=6.0)
    parser.add_argument("--read-length", type=int, default=10000)
    parser.add_argument("--replicates", type=int, default=25)
    parser.add_argument("--max-active-reads", type=int, default=20)
    parser.add_argument(
        "--heterozygosities",
        type=float,
        nargs="+",
        default=[2e-4, 5e-4, 1e-3, 2e-3],
    )
    parser.add_argument(
        "--error-rates", type=float, nargs="+", default=[0.05, 0.10, 0.20]
    )
    parser.add_argument(
        "--include-unspanned", action="store_true",
        help="rank block breaks too, as the 10 Sep table did; wrong, kept to reproduce it",
    )
    args = parser.parse_args()

    rows = run(
        args.heterozygosities,
        args.error_rates,
        args.replicates,
        args.variants,
        args.coverage,
        args.read_length,
        args.max_active_reads,
        include_unspanned=args.include_unspanned,
    )
    print(
        f"{args.variants} variants, physical coverage {args.coverage:g}, read "
        f"length {args.read_length} bp, {args.replicates} replicates, "
        f"full coverage to the solver"
    )
    print(
        "connectivity is an outcome, not a control: the switch rate is over "
        "the comparable pairs that survived"
    )
    print("simulated data; variant sites are all truly heterozygous")
    print("informativeness ranks " + ("every boundary, block breaks included (the 10 Sep error)"
          if args.include_unspanned else "spanned boundaries only"))
    header = (
        f"{'het':>8} {'err':>6} {'var/read':>9} {'depth':>6} {'blocks':>7} "
        f"{'phased':>7} {'pairs':>7} {'gap':>8} {'truth opt':>10} {'switch':>8} "
        f"{'informative':>12} {'posterior':>10} {'n':>4} {'skip':>5}"
    )
    print()
    print(header)
    print("-" * len(header))
    previous = None
    for row in rows:
        if previous is not None and row["heterozygosity"] != previous:
            print()
        previous = row["heterozygosity"]
        print(
            f"{row['heterozygosity']:>8.1e} {row['error_rate']:>6.2f} "
            f"{row['per_read']:>9.2f} {row['depth']:>6.2f} {row['blocks']:>7.2f} "
            f"{row['phased']:>7.3f} {row['pairs']:>7.1f} {row['gap']:>8.3f} "
            f"{row['truth_optimal']:>10.2f} {row['switch']:>8.4f} "
            f"{row['informative']:>12.3f} {row['posterior_informative']:>10.3f} "
            f"{row['with_a_switch']:>4} "
            f"{row['skipped']:>5}"
        )


if __name__ == "__main__":
    main()
