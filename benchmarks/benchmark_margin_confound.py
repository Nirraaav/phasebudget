"""Is the adaptive margin signal the phasing's fragility or the subset's thinness?

This is OQ-10. The adaptive loop measures margins on the subset that selection
just chose, so a small margin can mean the phase relation is genuinely
uncertain or merely that the subset is thin there. The two are separated by
measuring the margins a second time on the whole read set, which is affordable
only in the regime where the full matrix is under the dynamic program's state
limit. Every number here therefore comes from instances that did not need
selection, and it does not transfer to the coverage-constrained regime by
itself. `docs/design/margin-source.md` states the argument.

Three quantities per cell:

  * agreement: the fraction of full-matrix spanned boundaries the subset
    reports unspanned, and the rank correlation of the two margin vectors over
    the boundaries both report spanned.
  * informativeness: the rank-sum probability that a boundary carrying a switch
    error has the smaller margin, computed for each source against the same
    switch positions.
  * accuracy: switch error of the adaptive loop under each margin source, with
    static selection on the same instances as the baseline. The two sources do
    not run the same number of rounds, because a constant margin vector still
    moves the selection until the boosted weights stop changing it, so the
    round count is reported and a matched pair at one reweight each is reported
    beside it.
"""

from __future__ import annotations

import argparse
import statistics

from phasebudget._validation import CoverageLimitExceeded
from phasebudget.adaptive import adaptive_selection, static_selection
from phasebudget.dp import CoverageDP
from phasebudget.metrics import phasing_from_solution, switch_error_rate, switch_positions
from phasebudget.selection import read_weights, select_optimal
from phasebudget.simulate import simulate_connected


def rank_probability(bad: list[float], good: list[float]) -> float:
    """P(margin at a switch < margin elsewhere), ties counted as one half.

    The same statistic as `benchmark_objective_truth.rank_probability`, applied
    here to two different signals on one set of switch positions.
    """
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


def _ranks(values: list[float]) -> list[float]:
    """Average ranks, so tied margins do not invent an ordering."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    index = 0
    while index < len(order):
        stop = index
        while stop + 1 < len(order) and values[order[stop + 1]] == values[order[index]]:
            stop += 1
        shared = (index + stop) / 2.0 + 1.0
        for position in range(index, stop + 1):
            ranks[order[position]] = shared
        index = stop + 1
    return ranks


def spearman(left: list[float], right: list[float]) -> float:
    """Pearson correlation of the average ranks.

    Returns `nan` for fewer than three points or when either side is constant,
    which happens whenever every margin in the comparable set is equal.
    """
    if len(left) != len(right):
        raise ValueError("the two margin vectors must have the same length")
    if len(left) < 3:
        return float("nan")
    a, b = _ranks(left), _ranks(right)
    mean_a, mean_b = statistics.fmean(a), statistics.fmean(b)
    da = [value - mean_a for value in a]
    db = [value - mean_b for value in b]
    denominator = (sum(x * x for x in da) * sum(x * x for x in db)) ** 0.5
    if denominator == 0.0:
        return float("nan")
    return sum(x * y for x, y in zip(da, db)) / denominator


def _switch_of(instance, matrix, result) -> float:
    """Switch error of one adaptive result, against the instance's truth."""
    return switch_error_rate(
        instance.haplotype,
        phasing_from_solution(
            matrix.subset(result.selection.selected), result.solution
        ),
    )


def _mean(values: list[float]) -> float:
    return statistics.fmean(values) if values else float("nan")


def _worst(values: list[float], *, low: bool = True) -> float:
    if not values:
        return float("nan")
    return min(values) if low else max(values)


def run(error_rates, caps, replicates, n_variants, coverage, span, iterations):
    """One row per (error rate, cap). Instances over the state limit are skipped."""
    rows = []
    for error_rate in error_rates:
        for cap in caps:
            hidden: list[float] = []
            rhos: list[float] = []
            informative_subset: list[float] = []
            informative_full: list[float] = []
            switch_static: list[float] = []
            switch_subset: list[float] = []
            switch_full: list[float] = []
            matched_subset: list[float] = []
            matched_full: list[float] = []
            rounds_subset: list[float] = []
            rounds_full: list[float] = []
            skipped = 0
            used = 0
            for replicate in range(replicates):
                instance = simulate_connected(
                    n_variants,
                    coverage=coverage,
                    read_span=span,
                    error_rate=error_rate,
                    seed=9000 + 17 * replicate,
                )
                matrix = instance.matrix
                try:
                    full_dp = CoverageDP(matrix)
                except CoverageLimitExceeded:
                    skipped += 1
                    continue
                used += 1
                full_margins = full_dp.margins()

                selection = select_optimal(matrix, cap, read_weights(matrix))
                subset = matrix.subset(selection.selected)
                subset_dp = CoverageDP(subset)
                subset_solution = subset_dp.solve()
                subset_margins = subset_dp.margins()

                # Agreement. A boundary the full matrix spans and the subset
                # does not is the structural half of the confound: the subset
                # loop cannot see it at all.
                full_spanned = [m for m in full_margins if m.spanned]
                if full_spanned:
                    lost = sum(
                        1 for m in full_spanned if not subset_margins[m.boundary].spanned
                    )
                    hidden.append(lost / len(full_spanned))
                both = [
                    m.boundary
                    for m in full_margins
                    if m.spanned and subset_margins[m.boundary].spanned
                ]
                rho = spearman(
                    [full_margins[j].margin for j in both],
                    [subset_margins[j].margin for j in both],
                )
                if rho == rho:
                    rhos.append(rho)

                # Informativeness. Both signals are scored against the switch
                # positions of the subset phasing, which is what the loop is
                # trying to repair.
                phasing = phasing_from_solution(subset, subset_solution)
                switches = set(switch_positions(instance.haplotype, phasing))
                for margins, sink in (
                    (subset_margins, informative_subset),
                    (full_margins, informative_full),
                ):
                    # Only boundaries a read spans are ranked: an unspanned one
                    # has margin zero by construction and no switch can be
                    # counted there (corrected 26 Sep, see the development notebook).
                    bad = [m.margin for m in margins
                           if m.spanned and m.boundary + 1 in switches]
                    good = [m.margin for m in margins
                            if m.spanned and m.boundary + 1 not in switches]
                    probability = rank_probability(bad, good)
                    if probability == probability:
                        sink.append(probability)

                # Accuracy.
                static = static_selection(matrix, cap)
                switch_static.append(
                    switch_error_rate(
                        instance.haplotype,
                        phasing_from_solution(
                            matrix.subset(static.selection.selected), static.solution
                        ),
                    )
                )
                for source, sink, rounds, matched in (
                    ("subset", switch_subset, rounds_subset, matched_subset),
                    ("full", switch_full, rounds_full, matched_full),
                ):
                    loop = adaptive_selection(
                        matrix, cap, iterations=iterations, margin_source=source
                    )
                    sink.append(_switch_of(instance, matrix, loop))
                    rounds.append(float(loop.iterations))
                    # One reweight each, so the round budget cannot explain a
                    # difference between the two sources.
                    matched.append(
                        _switch_of(
                            instance,
                            matrix,
                            adaptive_selection(
                                matrix, cap, iterations=2, margin_source=source
                            ),
                        )
                    )
            rows.append(
                {
                    "error_rate": error_rate,
                    "cap": cap,
                    "hidden": _mean(hidden),
                    "hidden_worst": _worst(hidden, low=False),
                    "rho": _mean(rhos),
                    "rho_worst": _worst(rhos),
                    "informative_subset": _mean(informative_subset),
                    "informative_subset_worst": _worst(informative_subset),
                    "informative_full": _mean(informative_full),
                    "informative_full_worst": _worst(informative_full),
                    "switch_static": _mean(switch_static),
                    "switch_subset": _mean(switch_subset),
                    "switch_full": _mean(switch_full),
                    "matched_subset": _mean(matched_subset),
                    "matched_full": _mean(matched_full),
                    "rounds_subset": _mean(rounds_subset),
                    "rounds_full": _mean(rounds_full),
                    "instances": used,
                    "skipped": skipped,
                    "with_a_switch": len(informative_subset),
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=int, default=30)
    parser.add_argument("--coverage", type=float, default=8.0)
    parser.add_argument("--span", type=int, default=5)
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--caps", type=int, nargs="+", default=[3, 5])
    parser.add_argument(
        "--error-rates", type=float, nargs="+", default=[0.05, 0.10, 0.20]
    )
    parser.add_argument("--iterations", type=int, default=4)
    args = parser.parse_args()

    rows = run(
        args.error_rates,
        args.caps,
        args.replicates,
        args.variants,
        args.coverage,
        args.span,
        args.iterations,
    )
    print(
        f"{args.variants} variants, source coverage {args.coverage:g}, span "
        f"{args.span}, {args.replicates} replicates, caps {args.caps}, "
        f"{args.iterations} iterations"
    )
    print(
        "the full-margin source solves the whole matrix, so every instance here "
        "is one selection was not needed for"
    )
    print("simulated data; all sites are preselected heterozygous variants")

    agreement = (
        f"{'err':>6} {'cap':>4} {'hidden':>8} {'worst':>7} {'rho':>7} {'worst':>7} "
        f"{'inf sub':>8} {'worst':>7} {'inf full':>9} {'worst':>7} {'n':>4}"
    )
    print()
    print("agreement and informativeness")
    print(agreement)
    print("-" * len(agreement))
    for row in rows:
        print(
            f"{row['error_rate']:>6.2f} {row['cap']:>4} {row['hidden']:>8.3f} "
            f"{row['hidden_worst']:>7.3f} {row['rho']:>7.3f} {row['rho_worst']:>7.3f} "
            f"{row['informative_subset']:>8.3f} {row['informative_subset_worst']:>7.3f} "
            f"{row['informative_full']:>9.3f} {row['informative_full_worst']:>7.3f} "
            f"{row['with_a_switch']:>4}"
        )

    accuracy = (
        f"{'err':>6} {'cap':>4} {'static':>8} {'sub':>8} {'full':>8} "
        f"{'sub@2':>8} {'full@2':>8} {'r sub':>6} {'r full':>7} {'n':>4} {'skip':>5}"
    )
    print()
    print("switch error rate by margin source; @2 is one reweight each")
    print(accuracy)
    print("-" * len(accuracy))
    for row in rows:
        print(
            f"{row['error_rate']:>6.2f} {row['cap']:>4} {row['switch_static']:>8.4f} "
            f"{row['switch_subset']:>8.4f} {row['switch_full']:>8.4f} "
            f"{row['matched_subset']:>8.4f} {row['matched_full']:>8.4f} "
            f"{row['rounds_subset']:>6.2f} {row['rounds_full']:>7.2f} "
            f"{row['instances']:>4} {row['skipped']:>5}"
        )


if __name__ == "__main__":
    main()
