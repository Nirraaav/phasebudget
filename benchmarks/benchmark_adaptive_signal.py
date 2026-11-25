"""Does the posterior margin steer the adaptive loop better than the MEC margin?

Three strategies on the same instances, through `benchmark_phasing.phase`:
exact static selection, the loop steered by MEC margins, and the loop steered
by posterior margins (`docs/design/posterior-margins.md`), and a control that
steers by the posterior margins shuffled across boundaries. Differences in
switch error are paired per instance and printed with a 95 per cent interval;
negative means the first strategy of the pair phases better. The loop's
boost, quantile and iteration cap are the library defaults for both signals.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
from pathlib import Path

from phasebudget.metrics import phasing_from_solution, summarize
from phasebudget.simulate import simulate_connected

PAIRS = (
    ("adaptive_posterior", "adaptive"),
    ("adaptive", "exact"),
    ("adaptive_posterior", "exact"),
    ("adaptive", "adaptive_shuffled"),
    ("adaptive_shuffled", "exact"),
)


def _phasing_module():
    path = Path(__file__).with_name("benchmark_phasing.py")
    spec = importlib.util.spec_from_file_location("phasing", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _interval(values):
    n = len(values)
    mean = sum(values) / n
    variance = sum((v - mean) ** 2 for v in values) / (n - 1) if n > 1 else 0.0
    error = math.sqrt(variance / n)
    return mean, mean - 1.96 * error, mean + 1.96 * error


def run(caps, error_rates, replicates, n_variants, coverage, span, iterations):
    phasing = _phasing_module()
    rows = []
    for error_rate in error_rates:
        for cap in caps:
            switch = {name: [] for name in
                      ("exact", "adaptive", "adaptive_posterior", "adaptive_shuffled")}
            for replicate in range(replicates):
                instance = simulate_connected(
                    n_variants,
                    coverage=coverage,
                    read_span=span,
                    error_rate=error_rate,
                    seed=7000 + 13 * replicate,
                )
                for name in switch:
                    selection, solution, subset, _ = phasing.phase(
                        instance.matrix, cap, name, replicate, iterations
                    )
                    report = summarize(
                        instance.haplotype, phasing_from_solution(subset, solution)
                    )
                    switch[name].append(report["switch_error_rate"])
            row = {
                "error_rate": error_rate,
                "cap": cap,
                **{f"mean_{name}": sum(v) / len(v) for name, v in switch.items()},
            }
            for first, second in PAIRS:
                differences = [a - b for a, b in zip(switch[first], switch[second])]
                row[(first, second)] = (
                    *_interval(differences),
                    sum(d < 0 for d in differences),
                    sum(d > 0 for d in differences),
                )
            rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=int, default=60)
    parser.add_argument("--coverage", type=float, default=20.0)
    parser.add_argument("--span", type=int, default=6)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--caps", type=int, nargs="+", default=[3, 5, 8])
    parser.add_argument("--error-rates", type=float, nargs="+", default=[0.05, 0.10, 0.20])
    parser.add_argument("--iterations", type=int, default=4)
    args = parser.parse_args()
    print(
        f"{args.variants} variants, depth {args.coverage:g}, span {args.span}, "
        f"{args.replicates} replicates, {args.iterations} iterations, seed 7000 + 13 * "
        "replicate; boost and quantile at library defaults; simulated data"
    )
    print("switch error means, then paired differences with 95% intervals and wins")
    for row in run(args.caps, args.error_rates, args.replicates, args.variants,
                   args.coverage, args.span, args.iterations):
        print(
            f"err {row['error_rate']:.2f} cap {row['cap']}: exact {row['mean_exact']:.4f}, "
            f"adaptive {row['mean_adaptive']:.4f}, posterior {row['mean_adaptive_posterior']:.4f}"
        )
        for first, second in PAIRS:
            mean, low, high, first_wins, second_wins = row[(first, second)]
            print(
                f"    {first} - {second}: {mean:+.4f} [{low:+.4f}, {high:+.4f}] "
                f"wins {first_wins}:{second_wins}"
            )


if __name__ == "__main__":
    main()
