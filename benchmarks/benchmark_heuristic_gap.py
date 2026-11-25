"""How far deployed read selection sits from the certified optimum.

This is the headline measurement. `select_optimal` certifies the maximum of
the additive selection objective under the coverage cap, so every other
selector's objective can be divided by it and the quotient is an optimality
gap rather than a comparison. The number exists whatever the outcome, and
nothing in the project depends on it coming out large.

Four levels are reported and they are kept apart, because a selector that wins
one can lose another:

  1. **selection objective**, as a fraction of the certified optimum. Mean and
     worst case over replicates, since a mean ratio hides the instance that
     went wrong.
  2. **runtime** of the selection step alone, in milliseconds, separately from
     the solve it enables.
  3. **wMEC cost** of the reported haplotype, scored against the full read
     set. Also printed against the retained read set, which is not comparable
     across selectors and is here to show that it is not.
  4. **switch error** against the true haplotype.

Two versions of the Fischer and Marschall selector are measured.
`fischer_marschall` follows the paper's pseudocode; `whatshap` is a port of
WhatsHap 2.8's `readselect.pyx`, which selects the same reads as WhatsHap on
every instance checked (`tests/test_whatshap_port.py`) and is the deployed
selector. The Fischer and Marschall selectors and the maximum-minimum-coverage selector
maximize other objectives, stated in
`docs/design/published-selection-objectives.md`. Their ratio column is the
distance from the additive optimum, not a claim that they tried to reach it.
The adaptive loop is absent on purpose: it reweights, so it is answering a
different question, and `benchmark_phasing.py` is where it belongs.

Coverage is span coverage throughout, which is the convention fixed in
`docs/design/coverage-convention.md`.
"""

from __future__ import annotations

import argparse
import platform
import statistics
import subprocess
import time

from phasebudget.dp import CoverageDP
from phasebudget.metrics import phasing_from_solution, summarize
from phasebudget.selection import (
    distinct_variants_covered,
    min_coverage,
    read_weights,
    select_fischer_marschall,
    select_greedy,
    select_max_min_coverage,
    select_optimal,
    select_random,
    select_whatshap,
    selection_objective,
)
from phasebudget.simulate import simulate_connected

SELECTORS = (
    "random",
    "greedy",
    "fischer_marschall",
    "whatshap",
    "max_min_coverage",
    "optimal",
)


def choose(matrix, cap, name, weights, seed, positions=None):
    if name == "random":
        return select_random(matrix, cap, seed=seed, weights=weights)
    if name == "greedy":
        return select_greedy(matrix, cap, weights)
    if name == "fischer_marschall":
        return select_fischer_marschall(matrix, cap, weights)
    if name == "whatshap":
        return select_whatshap(matrix, cap, weights, positions=positions)
    if name == "max_min_coverage":
        return select_max_min_coverage(matrix, cap, weights)
    return select_optimal(matrix, cap, weights)


def measure(instance, cap, name, seed):
    """One selector on one instance, at all four levels."""
    matrix = instance.matrix
    weights = read_weights(matrix)
    started = time.perf_counter()
    selection = choose(matrix, cap, name, weights, seed, instance.positions)
    select_ms = 1000 * (time.perf_counter() - started)

    subset = matrix.subset(selection.selected)
    started = time.perf_counter()
    solution = CoverageDP(subset).solve()
    solve_ms = 1000 * (time.perf_counter() - started)

    phasing = phasing_from_solution(subset, solution)
    report = summarize(instance.haplotype, phasing, instance.positions)
    return {
        # Every selection is scored under the base weights, never under a
        # selector's own notion of value.
        "objective": selection_objective(weights, selection.selected),
        "select_ms": select_ms,
        "solve_ms": solve_ms,
        "mec_full": matrix.cost(solution.haplotype),
        "mec_kept": solution.cost,
        "switch": report["switch_error_rate"],
        "flips": report["flip_errors"],
        "n50": report["block_n50"],
        "qan50": report["qan50"],
        "phased": report["phased_fraction"],
        # Each selector's own objective gets a column, so the ratio above is
        # read next to what the selector was actually trying to do: `floor` is
        # the maximum-minimum-coverage objective and `covered` is the quantity
        # the Fischer and Marschall slice pays for.
        "floor": float(min_coverage(matrix, selection.selected)),
        "covered": float(distinct_variants_covered(matrix, selection.selected)),
        "reads": len(selection.selected),
    }


def run(
    caps,
    error_rates,
    replicates,
    n_variants,
    coverage,
    span,
    seed_base=7000,
    quality_spread=0.0,
):
    rows = []
    for error_rate in error_rates:
        for cap in caps:
            per_selector = {
                name: {key: [] for key in
                       ("objective", "ratio", "select_ms", "solve_ms", "mec_full",
                        "mec_kept", "switch", "flips", "n50", "qan50", "phased",
                        "floor", "covered", "reads")}
                for name in SELECTORS
            }
            for replicate in range(replicates):
                instance = simulate_connected(
                    n_variants,
                    coverage=coverage,
                    read_span=span,
                    error_rate=error_rate,
                    seed=seed_base + 13 * replicate,
                    quality_spread=quality_spread,
                )
                measured = {
                    name: measure(instance, cap, name, replicate)
                    for name in SELECTORS
                }
                certified = measured["optimal"]["objective"]
                for name in SELECTORS:
                    record = measured[name]
                    # Per-instance ratio against the optimum on that same
                    # instance. Averaging the numerator and denominator
                    # separately would report a ratio of means instead.
                    ratio = (
                        record["objective"] / certified if certified > 0 else 1.0
                    )
                    for key, value in record.items():
                        per_selector[name][key].append(value)
                    per_selector[name]["ratio"].append(ratio)
            for name in SELECTORS:
                data = per_selector[name]
                rows.append(
                    {
                        "error_rate": error_rate,
                        "cap": cap,
                        "selector": name,
                        "objective": statistics.mean(data["objective"]),
                        "ratio": statistics.mean(data["ratio"]),
                        "ratio_worst": min(data["ratio"]),
                        "select_ms": statistics.mean(data["select_ms"]),
                        "solve_ms": statistics.mean(data["solve_ms"]),
                        "mec_full": statistics.mean(data["mec_full"]),
                        "mec_kept": statistics.mean(data["mec_kept"]),
                        "switch": statistics.mean(data["switch"]),
                        "switch_worst": max(data["switch"]),
                        "flips": statistics.mean(data["flips"]),
                        "n50": statistics.mean(data["n50"]),
                        "qan50": statistics.mean(data["qan50"]),
                        "phased": statistics.mean(data["phased"]),
                        "floor": statistics.mean(data["floor"]),
                        "covered": statistics.mean(data["covered"]),
                        "reads": statistics.mean(data["reads"]),
                    }
                )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=int, default=60)
    parser.add_argument("--coverage", type=float, default=20.0)
    parser.add_argument("--span", type=int, default=6)
    parser.add_argument("--replicates", type=int, default=25)
    parser.add_argument("--caps", type=int, nargs="+", default=[3, 5, 8])
    parser.add_argument(
        "--error-rates", type=float, nargs="+", default=[0.01, 0.05, 0.10, 0.20]
    )
    parser.add_argument("--seed-base", type=int, default=7000)
    parser.add_argument(
        "--quality-spread",
        type=float,
        default=0.0,
        help=(
            "phred jitter per base. At 0 every base carries the same log-odds "
            "weight, so a read's score is its called-allele count and the "
            "additive objective has few distinct values"
        ),
    )
    args = parser.parse_args()

    rows = run(
        args.caps,
        args.error_rates,
        args.replicates,
        args.variants,
        args.coverage,
        args.span,
        args.seed_base,
        args.quality_spread,
    )

    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    print(f"commit {commit}; {platform.system()} {platform.release()} {platform.machine()}; "
          f"Python {platform.python_version()}")
    print(
        f"{args.variants} variants, depth {args.coverage:g}, span {args.span}, "
        f"caps {' '.join(str(c) for c in args.caps)}, "
        f"error rates {' '.join(f'{e:g}' for e in args.error_rates)}, "
        f"{args.replicates} replicates, quality spread {args.quality_spread:g}, "
        f"seed {args.seed_base} + 13 * replicate"
    )
    print("simulated data; span coverage; every selection scored under base weights")
    print(
        "obj/opt is the additive objective over the certified optimum on the "
        "same instance; worst is the worst replicate"
    )
    print(
        "MEC full is comparable across selectors, MEC kept is not: it is the "
        "cost on the retained reads, and keeping fewer reads lowers it"
    )
    print(
        "floor and covered are the objectives of the max-min and "
        "Fischer--Marschall selectors, printed so the ratio is read against "
        "what each one optimizes"
    )
    print(
        "flips is the mean count of flip errors; N50 and QAN50 are in base pairs "
        "of the simulated positions"
    )
    header = (
        f"{'err':>5} {'cap':>4} {'selector':>18} {'obj/opt':>8} {'worst':>7} "
        f"{'sel ms':>7} {'dp ms':>7} {'MEC full':>9} {'MEC kept':>9} "
        f"{'switch':>8} {'sw worst':>9} {'flips':>6} {'N50':>8} {'QAN50':>8} "
        f"{'phased':>7} {'floor':>6} {'covered':>8} {'reads':>6}"
    )
    print(header)
    print("-" * len(header))
    previous = None
    for row in rows:
        key = (row["error_rate"], row["cap"])
        if previous is not None and key != previous:
            print()
        previous = key
        print(
            f"{row['error_rate']:>5.2f} {row['cap']:>4} {row['selector']:>18} "
            f"{row['ratio']:>8.4f} {row['ratio_worst']:>7.4f} "
            f"{row['select_ms']:>7.2f} {row['solve_ms']:>7.1f} "
            f"{row['mec_full']:>9.2f} {row['mec_kept']:>9.2f} "
            f"{row['switch']:>8.4f} {row['switch_worst']:>9.4f} "
            f"{row['flips']:>6.2f} {row['n50']:>8.0f} {row['qan50']:>8.0f} "
            f"{row['phased']:>7.3f} {row['floor']:>6.2f} "
            f"{row['covered']:>8.1f} {row['reads']:>6.1f}"
        )


if __name__ == "__main__":
    main()
