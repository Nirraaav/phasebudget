"""Every selector scored against every selection objective.

Three objectives are in play and they are not the same problem:

  * `weight`, the additive sum of read scores, maximized exactly by
    `select_optimal` through min-cost flow;
  * `floor`, the smallest per-variant coverage, maximized exactly by
    `select_max_min_coverage` (Mäkinen, Staneva, Tomescu and Valenzuela 2015);
  * `covered` and `blocks`, the number of distinct variants covered and the
    number of connected components left, which is what the slice-and-bridge
    procedure of Fischer and Marschall 2016 is built to improve.

Reporting one selector against another selector's objective is the only way to
see how far apart the objectives are. Phasing accuracy is not measured here;
`benchmark_phasing.py` does that.
"""

from __future__ import annotations

import argparse
import statistics
import time

from phasebudget.selection import (
    component_count,
    distinct_variants_covered,
    min_coverage,
    read_weights,
    select_fischer_marschall,
    select_greedy,
    select_max_min_coverage,
    select_optimal,
    select_random,
    select_whatshap,
)
from phasebudget.simulate import simulate

SELECTORS = (
    "random",
    "greedy",
    "fischer_marschall",
    "whatshap",
    "max_min_coverage",
    "optimal",
)


def choose(matrix, coverage, name, weights, seed):
    if name == "random":
        return select_random(matrix, coverage, seed=seed, weights=weights)
    if name == "greedy":
        return select_greedy(matrix, coverage, weights)
    if name == "fischer_marschall":
        return select_fischer_marschall(matrix, coverage, weights)
    if name == "whatshap":
        return select_whatshap(matrix, coverage, weights)
    if name == "max_min_coverage":
        return select_max_min_coverage(matrix, coverage, weights)
    return select_optimal(matrix, coverage, weights)


def run(caps, coverages, spans, replicates, n_variants, error_rate):
    rows = []
    for coverage in coverages:
        for span in spans:
            for cap in caps:
                totals = {
                    name: {
                        "weight": [],
                        "floor": [],
                        "covered": [],
                        "blocks": [],
                        "reads": [],
                        "ms": [],
                    }
                    for name in SELECTORS
                }
                for replicate in range(replicates):
                    instance = simulate(
                        n_variants,
                        coverage=coverage,
                        read_span=span,
                        error_rate=error_rate,
                        seed=1000 * replicate + span,
                    )
                    matrix = instance.matrix
                    weights = read_weights(matrix)
                    for name in SELECTORS:
                        started = time.perf_counter()
                        result = choose(matrix, cap, name, weights, replicate)
                        elapsed = time.perf_counter() - started
                        chosen = result.selected
                        totals[name]["weight"].append(result.objective)
                        totals[name]["floor"].append(min_coverage(matrix, chosen))
                        totals[name]["covered"].append(
                            distinct_variants_covered(matrix, chosen)
                        )
                        totals[name]["blocks"].append(component_count(matrix, chosen))
                        totals[name]["reads"].append(len(chosen))
                        totals[name]["ms"].append(1000 * elapsed)
                best_weight = max(
                    statistics.mean(totals[name]["weight"]) for name in SELECTORS
                )
                for name in SELECTORS:
                    data = totals[name]
                    rows.append(
                        {
                            "coverage": coverage,
                            "span": span,
                            "cap": cap,
                            "selector": name,
                            "weight": statistics.mean(data["weight"]),
                            "weight_ratio": (
                                statistics.mean(data["weight"]) / best_weight
                                if best_weight
                                else 1.0
                            ),
                            "floor": statistics.mean(data["floor"]),
                            "covered": statistics.mean(data["covered"]),
                            "blocks": statistics.mean(data["blocks"]),
                            "reads": statistics.mean(data["reads"]),
                            "ms": statistics.mean(data["ms"]),
                        }
                    )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=int, default=200)
    parser.add_argument("--replicates", type=int, default=8)
    parser.add_argument("--error-rate", type=float, default=0.03)
    parser.add_argument("--caps", type=int, nargs="+", default=[2, 4, 8])
    parser.add_argument("--coverages", type=float, nargs="+", default=[10.0, 30.0])
    parser.add_argument("--spans", type=int, nargs="+", default=[4, 12])
    args = parser.parse_args()

    rows = run(
        args.caps,
        args.coverages,
        args.spans,
        args.replicates,
        args.variants,
        args.error_rate,
    )
    header = (
        f"{'cov':>4} {'span':>5} {'cap':>4} {'selector':>18} {'weight':>9} "
        f"{'w/best':>7} {'floor':>6} {'covered':>8} {'blocks':>7} "
        f"{'reads':>7} {'ms':>8}"
    )
    print(
        f"{args.variants} variants, error rate {args.error_rate:g}, "
        f"{args.replicates} replicates per cell"
    )
    print("each selector optimizes a different column; none optimizes them all")
    print(header)
    print("-" * len(header))
    previous = None
    for row in rows:
        key = (row["coverage"], row["span"], row["cap"])
        if previous is not None and key != previous:
            print()
        previous = key
        print(
            f"{row['coverage']:>4.0f} {row['span']:>5} {row['cap']:>4} "
            f"{row['selector']:>18} {row['weight']:>9.1f} "
            f"{row['weight_ratio']:>7.3f} {row['floor']:>6.2f} "
            f"{row['covered']:>8.1f} {row['blocks']:>7.2f} "
            f"{row['reads']:>7.1f} {row['ms']:>8.2f}"
        )


if __name__ == "__main__":
    main()
