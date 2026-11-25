"""Does the wMEC optimum track the truth, and do margins know where it fails?

Two quantities per instance, both at full coverage so that read selection is
not in the way:

  * the objective-truth gap: the wMEC cost of the true haplotype minus the
    optimum. It is zero when the truth is an optimal solution and grows when
    the objective prefers something else.
  * margin informativeness: the rank-sum probability that a boundary carrying
    a switch error has a smaller margin than a boundary that does not. A value
    of 0.5 means margins say nothing about where the phasing is wrong.

The project's central hypothesis is that the second decays as the first grows.
"""

from __future__ import annotations

import argparse
import statistics

from phasebudget.dp import CoverageDP
from phasebudget.metrics import phasing_from_solution, switch_positions
from phasebudget.simulate import simulate_connected


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


def run(error_rates, replicates, n_variants, coverages, spans):
    rows = []
    for coverage in coverages:
        for span in spans:
            for error_rate in error_rates:
                gaps = []
                switch_rates = []
                probabilities = []
                truth_optimal = 0
                for replicate in range(replicates):
                    instance = simulate_connected(
                        n_variants,
                        coverage=coverage,
                        read_span=span,
                        error_rate=error_rate,
                        seed=31 * replicate + 5,
                    )
                    matrix = instance.matrix
                    dp = CoverageDP(matrix)
                    solution = dp.solve()
                    truth_cost = matrix.cost(instance.haplotype)
                    gaps.append(truth_cost - solution.cost)
                    truth_optimal += truth_cost <= solution.cost + 1e-9

                    phasing = phasing_from_solution(matrix, solution)
                    switches = set(switch_positions(instance.haplotype, phasing))
                    margins = dp.margins()
                    comparable = len(margins)
                    switch_rates.append(len(switches) / comparable if comparable else 0.0)
                    bad = [m.margin for m in margins if m.boundary + 1 in switches]
                    good = [m.margin for m in margins if m.boundary + 1 not in switches]
                    probability = rank_probability(bad, good)
                    if probability == probability:
                        probabilities.append(probability)
                rows.append(
                    {
                        "coverage": coverage,
                        "span": span,
                        "error_rate": error_rate,
                        "gap": statistics.mean(gaps),
                        "truth_optimal": truth_optimal / replicates,
                        "switch": statistics.mean(switch_rates),
                        "informative": (
                            statistics.mean(probabilities)
                            if probabilities
                            else float("nan")
                        ),
                        "instances_with_a_switch": len(probabilities),
                    }
                )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=int, default=40)
    parser.add_argument("--coverage", type=float, help="one source coverage")
    parser.add_argument("--coverages", type=float, nargs="+", default=[10.0])
    parser.add_argument("--span", type=int, help="one read span")
    parser.add_argument("--spans", type=int, nargs="+", default=[6])
    parser.add_argument("--replicates", type=int, default=30)
    parser.add_argument(
        "--error-rates",
        type=float,
        nargs="+",
        default=[0.01, 0.03, 0.05, 0.1, 0.15, 0.2, 0.3],
    )
    args = parser.parse_args()

    coverages = [args.coverage] if args.coverage is not None else args.coverages
    spans = [args.span] if args.span is not None else args.spans
    rows = run(args.error_rates, args.replicates, args.variants, coverages, spans)
    header = (
        f"{'cov':>4} {'span':>5} {'err':>6} {'gap':>9} {'truth opt':>10} {'switch':>8} "
        f"{'informative':>12} {'n':>4}"
    )
    print(
        f"{args.variants} variants, coverages {coverages}, spans {spans}, "
        f"{args.replicates} replicates, full coverage to the solver"
    )
    print("all sites are preselected heterozygous variants; heterozygosity is not varied")
    print(header)
    print("-" * len(header))
    previous = None
    for row in rows:
        key = (row["coverage"], row["span"])
        if previous is not None and key != previous:
            print()
        previous = key
        print(
            f"{row['coverage']:>4g} {row['span']:>5} {row['error_rate']:>6.2f} {row['gap']:>9.3f} "
            f"{row['truth_optimal']:>10.2f} {row['switch']:>8.4f} "
            f"{row['informative']:>12.3f} {row['instances_with_a_switch']:>4}"
        )


if __name__ == "__main__":
    main()
