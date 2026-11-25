"""Accuracy of the whole pipeline under each read-selection strategy.

Every strategy is given the same coverage cap and the same instances, so the
only difference is which reads reach the solver. Switch error is the accuracy
number; the selection objective is reported next to it so the two can be
compared directly.
"""

from __future__ import annotations

import argparse
import statistics
import time

from phasebudget.adaptive import adaptive_selection, static_selection
from phasebudget.dp import CoverageDP
from phasebudget.metrics import phasing_from_solution, summarize
from phasebudget.selection import (
    read_weights,
    selection_objective,
    select_fischer_marschall,
    select_greedy,
    select_max_min_coverage,
    select_random,
    select_whatshap,
)
from phasebudget.simulate import simulate_connected

STRATEGIES = (
    "random",
    "greedy",
    "fischer_marschall",
    "whatshap",
    "max_min_coverage",
    "exact",
    "adaptive",
    "adaptive_posterior",
    "adaptive_shuffled",
)


def phase(matrix, cap, strategy, seed, iterations):
    weights = read_weights(matrix)
    started = time.perf_counter()
    if strategy in ("adaptive", "adaptive_posterior", "adaptive_shuffled"):
        signal = {"adaptive": "margin", "adaptive_posterior": "posterior",
                  "adaptive_shuffled": "shuffled"}[strategy]
        result = adaptive_selection(matrix, cap, iterations=iterations, signal=signal,
                                    shuffle_seed=seed)
        selection, solution = result.selection, result.solution
    elif strategy == "exact":
        result = static_selection(matrix, cap)
        selection, solution = result.selection, result.solution
    else:
        if strategy == "random":
            selection = select_random(matrix, cap, seed=seed, weights=weights)
        elif strategy == "greedy":
            selection = select_greedy(matrix, cap, weights)
        elif strategy == "fischer_marschall":
            selection = select_fischer_marschall(matrix, cap, weights)
        elif strategy == "whatshap":
            selection = select_whatshap(matrix, cap, weights)
        else:
            selection = select_max_min_coverage(matrix, cap, weights)
        solution = CoverageDP(matrix.subset(selection.selected)).solve()
    elapsed = time.perf_counter() - started
    subset = matrix.subset(selection.selected)
    return selection, solution, subset, elapsed


def run(caps, error_rates, replicates, n_variants, coverage, span, iterations):
    rows = []
    for error_rate in error_rates:
        for cap in caps:
            per_strategy = {name: {"switch": [], "mec": [], "obj": [], "sec": []}
                            for name in STRATEGIES}
            for replicate in range(replicates):
                instance = simulate_connected(
                    n_variants,
                    coverage=coverage,
                    read_span=span,
                    error_rate=error_rate,
                    seed=7000 + 13 * replicate,
                )
                for name in STRATEGIES:
                    selection, solution, subset, elapsed = phase(
                        instance.matrix, cap, name, replicate, iterations
                    )
                    phasing = phasing_from_solution(subset, solution)
                    report = summarize(instance.haplotype, phasing, instance.positions)
                    per_strategy[name]["switch"].append(report["switch_error_rate"])
                    per_strategy[name]["mec"].append(
                        instance.matrix.cost(solution.haplotype)
                    )
                    # Score every selection under the base weights. The
                    # adaptive loop reports its own objective under boosted
                    # weights, which is not comparable across strategies.
                    per_strategy[name]["obj"].append(
                        selection_objective(
                            read_weights(instance.matrix), selection.selected
                        )
                    )
                    per_strategy[name]["sec"].append(elapsed)
            for name in STRATEGIES:
                data = per_strategy[name]
                rows.append(
                    {
                        "error_rate": error_rate,
                        "cap": cap,
                        "strategy": name,
                        "switch": statistics.mean(data["switch"]),
                        "mec": statistics.mean(data["mec"]),
                        "objective": statistics.mean(data["obj"]),
                        "ms": 1000 * statistics.mean(data["sec"]),
                    }
                )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=int, default=60)
    parser.add_argument("--coverage", type=float, default=20.0)
    parser.add_argument("--span", type=int, default=6)
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--caps", type=int, nargs="+", default=[3, 5, 8])
    parser.add_argument(
        "--error-rates", type=float, nargs="+", default=[0.01, 0.05, 0.1, 0.2]
    )
    parser.add_argument("--iterations", type=int, default=4)
    args = parser.parse_args()

    rows = run(
        args.caps,
        args.error_rates,
        args.replicates,
        args.variants,
        args.coverage,
        args.span,
        args.iterations,
    )
    header = (
        f"{'err':>6} {'cap':>4} {'strategy':>18} {'switch':>8} "
        f"{'MEC':>9} {'sel obj':>9} {'ms':>8}"
    )
    print(
        f"{args.variants} variants, coverage {args.coverage:g}, span {args.span}, "
        f"{args.replicates} replicates"
    )
    print("switch is the switch error rate over comparable adjacent pairs")
    print("the selection objective is additive weight; Fischer--Marschall and max-min optimize other outcomes")
    print(header)
    print("-" * len(header))
    previous = None
    for row in rows:
        key = (row["error_rate"], row["cap"])
        if previous is not None and key != previous:
            print()
        previous = key
        print(
            f"{row['error_rate']:>6.2f} {row['cap']:>4} {row['strategy']:>18} "
            f"{row['switch']:>8.4f} {row['mec']:>9.2f} "
            f"{row['objective']:>9.1f} {row['ms']:>8.1f}"
        )


if __name__ == "__main__":
    main()
