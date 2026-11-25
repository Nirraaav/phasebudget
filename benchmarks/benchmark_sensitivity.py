"""Is the adaptive gain a result or one lucky setting of three free parameters?

This is OQ-6. The development notebook reports adaptive selection beating static selection
at boost 2.0, quantile 0.25, and four iterations, which is one point in a
three-dimensional parameter space that was never swept. A method that wins only
there is not a method.

The sweep runs the whole grid on one set of instances and reports, per cell,
the switch error against the static baseline on those same instances. The
summary that matters is the fraction of settings that beat static and the worst
setting in the grid, not the best one: the best setting of a swept parameter is
a selection effect, and quoting it is how a null result gets published as a
positive one.
"""

from __future__ import annotations

import argparse
import statistics

from phasebudget.adaptive import adaptive_selection, static_selection
from phasebudget.metrics import phasing_from_solution, switch_error_rate
from phasebudget.simulate import simulate_connected


def _switch(instance, matrix, result) -> float:
    return switch_error_rate(
        instance.haplotype,
        phasing_from_solution(
            matrix.subset(result.selection.selected), result.solution
        ),
    )


def run(
    boosts,
    quantiles,
    iteration_caps,
    error_rates,
    caps,
    replicates,
    n_variants,
    coverage,
    span,
    margin_source,
):
    """One row per (error rate, cap, boost, quantile, iteration cap).

    The static baseline is computed once per instance and shared by every cell,
    so a difference between two rows is never a difference in instances.
    """
    rows = []
    for error_rate in error_rates:
        for cap in caps:
            instances = [
                simulate_connected(
                    n_variants,
                    coverage=coverage,
                    read_span=span,
                    error_rate=error_rate,
                    seed=4400 + 23 * replicate,
                )
                for replicate in range(replicates)
            ]
            baseline = [
                _switch(instance, instance.matrix, static_selection(instance.matrix, cap))
                for instance in instances
            ]
            static_mean = statistics.fmean(baseline)
            for boost in boosts:
                for quantile in quantiles:
                    for iterations in iteration_caps:
                        switches = []
                        wins = 0
                        rounds = []
                        for instance, reference in zip(instances, baseline):
                            result = adaptive_selection(
                                instance.matrix,
                                cap,
                                iterations=iterations,
                                boost=boost,
                                quantile=quantile,
                                margin_source=margin_source,
                            )
                            value = _switch(instance, instance.matrix, result)
                            switches.append(value)
                            wins += value < reference - 1e-12
                            rounds.append(float(result.iterations))
                        rows.append(
                            {
                                "error_rate": error_rate,
                                "cap": cap,
                                "boost": boost,
                                "quantile": quantile,
                                "iterations": iterations,
                                "static": static_mean,
                                "switch": statistics.fmean(switches),
                                "delta": statistics.fmean(switches) - static_mean,
                                "instances_improved": wins,
                                "rounds": statistics.fmean(rounds),
                                "replicates": replicates,
                            }
                        )
    return rows


def summarize(rows):
    """Per (error rate, cap): how much of the grid beats static, and the worst."""
    groups: dict[tuple[float, int], list[dict]] = {}
    for row in rows:
        groups.setdefault((row["error_rate"], row["cap"]), []).append(row)
    summary = []
    for (error_rate, cap), cells in sorted(groups.items()):
        deltas = [cell["delta"] for cell in cells]
        better = sum(1 for value in deltas if value < -1e-12)
        worse = sum(1 for value in deltas if value > 1e-12)
        summary.append(
            {
                "error_rate": error_rate,
                "cap": cap,
                "cells": len(cells),
                "static": cells[0]["static"],
                "better": better,
                "worse": worse,
                "best": min(deltas),
                "worst": max(deltas),
                "median": statistics.median(deltas),
            }
        )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=int, default=60)
    parser.add_argument("--coverage", type=float, default=20.0)
    parser.add_argument("--span", type=int, default=6)
    parser.add_argument("--replicates", type=int, default=12)
    parser.add_argument("--caps", type=int, nargs="+", default=[3, 5])
    parser.add_argument("--error-rates", type=float, nargs="+", default=[0.05, 0.10, 0.20])
    parser.add_argument("--boosts", type=float, nargs="+", default=[1.25, 1.5, 2.0, 4.0])
    parser.add_argument(
        "--quantiles", type=float, nargs="+", default=[0.1, 0.25, 0.5]
    )
    parser.add_argument("--iteration-caps", type=int, nargs="+", default=[2, 4, 8])
    parser.add_argument("--margin-source", default="subset", choices=("subset", "full"))
    parser.add_argument("--full-grid", action="store_true", help="print every cell")
    args = parser.parse_args()

    rows = run(
        args.boosts,
        args.quantiles,
        args.iteration_caps,
        args.error_rates,
        args.caps,
        args.replicates,
        args.variants,
        args.coverage,
        args.span,
        args.margin_source,
    )
    cells = len(args.boosts) * len(args.quantiles) * len(args.iteration_caps)
    print(
        f"{args.variants} variants, coverage {args.coverage:g}, span {args.span}, "
        f"{args.replicates} replicates, margin source {args.margin_source}"
    )
    print(
        f"boosts {args.boosts}, quantiles {args.quantiles}, iteration caps "
        f"{args.iteration_caps}: {cells} settings per error rate and cap"
    )
    print("delta is adaptive switch error minus static on the same instances; negative is better")
    print("simulated data")

    header = (
        f"{'err':>6} {'cap':>4} {'static':>8} {'cells':>6} {'better':>7} "
        f"{'worse':>6} {'best d':>8} {'median d':>9} {'worst d':>8}"
    )
    print()
    print("grid summary")
    print(header)
    print("-" * len(header))
    for row in summarize(rows):
        print(
            f"{row['error_rate']:>6.2f} {row['cap']:>4} {row['static']:>8.4f} "
            f"{row['cells']:>6} {row['better']:>7} {row['worse']:>6} "
            f"{row['best']:>8.4f} {row['median']:>9.4f} {row['worst']:>8.4f}"
        )

    if args.full_grid:
        detail = (
            f"{'err':>6} {'cap':>4} {'boost':>6} {'quant':>6} {'iter':>5} "
            f"{'switch':>8} {'delta':>8} {'improved':>9} {'rounds':>7}"
        )
        print()
        print("every setting")
        print(detail)
        print("-" * len(detail))
        for row in rows:
            print(
                f"{row['error_rate']:>6.2f} {row['cap']:>4} {row['boost']:>6.2f} "
                f"{row['quantile']:>6.2f} {row['iterations']:>5} "
                f"{row['switch']:>8.4f} {row['delta']:>8.4f} "
                f"{row['instances_improved']:>4}/{row['replicates']:<4} "
                f"{row['rounds']:>7.2f}"
            )


if __name__ == "__main__":
    main()
