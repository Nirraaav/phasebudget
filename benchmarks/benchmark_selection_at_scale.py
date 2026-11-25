"""Milestone 2: the optimality gap of read selection at chromosome scale.

`benchmark_heuristic_gap.py` measures the gap on 60-variant instances. This is
the same measurement, through its `measure` function, on simulated chromosome
stretches of thousands of variants with physical read lengths and densities
(`docs/design/genomic-simulator.md`), at the deployed cap and below. Two read
profiles stand in for long-read technologies: `hifi` is 15 kb at error 0.01
and `ont` is 10 kb at error 0.05; both carry a base-quality spread. The
objective ratio is against the certified optimum; mean and worst replicate are
both reported, and accuracy is a separate column.

HapCHAT's reduction (`phasebudget.hapchat`) merges reads before it selects, so
its reads are not the original ones and it has no objective ratio. It is
compared on switch error and on wMEC scored against the full matrix.
`merge_optimal` merges the same way and then selects the merged reads by the
exact min-cost flow instead of HapCHAT's greedy step. For both, `mixed` is the
share of merged reads that hold reads from both haplotypes.

`--systematic-rate`, `--systematic-strength`, `--burst-rate`, `--burst-length`
and `--burst-error` add errors the base qualities do not report
(`docs/design/genomic-simulator.md`); they are off by default.
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import os
import platform
import statistics
import subprocess
from multiprocessing import Pool
from pathlib import Path

import time

from phasebudget.dp import CoverageDP
from phasebudget.hapchat import merge_reads, select_hapchat
from phasebudget.selection import read_weights, select_optimal
from phasebudget.metrics import phasing_from_solution, summarize
from phasebudget.simulate import simulate_genomic

PROFILES = {
    "hifi": {"read_length": 15_000, "error_rate": 0.01},
    "ont": {"read_length": 10_000, "error_rate": 0.05},
    "noisy": {"read_length": 10_000, "error_rate": 0.10},
}
SELECTORS = ("random", "greedy", "fischer_marschall", "whatshap", "max_min_coverage", "optimal")
ALL = SELECTORS + ("hapchat", "merge_optimal")


def mixed_share(merged, instance):
    origin = dict(zip(instance.matrix.read_ids, instance.origin))
    clusters = [c for c in merged.members if len(c) > 1]
    if not clusters:
        return 0.0
    return sum(len({origin[r] for r in c}) > 1 for c in clusters) / len(clusters)


def measure_hapchat(instance, cap, seed, exact=False):
    started = time.perf_counter()
    if exact:
        merged = merge_reads(instance.matrix)
        selection = select_optimal(merged.matrix, cap, read_weights(merged.matrix))
    else:
        merged, selection = select_hapchat(instance.matrix, cap, seed=seed)
    select_ms = 1000 * (time.perf_counter() - started)
    subset = merged.matrix.subset(selection.selected)
    started = time.perf_counter()
    solution = CoverageDP(subset).solve()
    solve_ms = 1000 * (time.perf_counter() - started)
    report = summarize(instance.haplotype, phasing_from_solution(subset, solution),
                       instance.positions)
    return {
        "objective": float("nan"),
        "select_ms": select_ms,
        "solve_ms": solve_ms,
        "mec_full": instance.matrix.cost(solution.haplotype),
        "switch": report["switch_error_rate"],
        "flips": report["flip_errors"],
        "n50": report["block_n50"],
        "qan50": report["qan50"],
        "phased": report["phased_fraction"],
        "mixed": mixed_share(merged, instance),
    }


def _gap_module():
    path = Path(__file__).with_name("benchmark_heuristic_gap.py")
    spec = importlib.util.spec_from_file_location("heuristic_gap", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def unit(arguments):
    profile, cap, replicate, n_variants, depth, heterozygosity, spread, errors = arguments
    gap = _gap_module()
    instance = simulate_genomic(
        n_variants,
        coverage=depth,
        heterozygosity=heterozygosity,
        quality_spread=spread,
        seed=8000 + 17 * replicate,
        **PROFILES[profile],
        **errors,
    )
    measured = {name: gap.measure(instance, cap, name, replicate) for name in SELECTORS}
    measured["hapchat"] = measure_hapchat(instance, cap, replicate)
    measured["merge_optimal"] = measure_hapchat(instance, cap, replicate, exact=True)
    return profile, cap, replicate, len(instance.matrix), measured


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", nargs="+", default=["hifi", "ont"])
    parser.add_argument("--variants", type=int, default=2000)
    parser.add_argument("--depth", type=float, default=30.0)
    parser.add_argument("--heterozygosity", type=float, default=0.001)
    parser.add_argument("--quality-spread", type=float, default=2.0)
    parser.add_argument("--caps", type=int, nargs="+", default=[5, 10, 15])
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--processes", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    parser.add_argument("--systematic-rate", type=float, default=0.0)
    parser.add_argument("--systematic-strength", type=float, default=0.5)
    parser.add_argument("--burst-rate", type=float, default=0.0)
    parser.add_argument("--burst-length", type=int, default=5)
    parser.add_argument("--burst-error", type=float, default=0.35)
    args = parser.parse_args()
    errors = {
        "systematic_rate": args.systematic_rate,
        "systematic_strength": args.systematic_strength,
        "burst_rate": args.burst_rate,
        "burst_length": args.burst_length,
        "burst_error": args.burst_error,
    }

    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    print(f"commit {commit}; {platform.system()} {platform.release()} {platform.machine()}; "
          f"Python {platform.python_version()}")
    print(f"{args.variants} variants, depth {args.depth:g}, heterozygosity {args.heterozygosity:g}, "
          f"quality spread {args.quality_spread:g}, {args.replicates} replicates, seed 8000 + 17 * "
          f"replicate; profiles " + ", ".join(
              f"{p} ({PROFILES[p]['read_length']} bp, error {PROFILES[p]['error_rate']})"
              for p in args.profiles) + "; simulated data")
    print("unreported errors: " + ", ".join(f"{k} {v:g}" for k, v in errors.items()))
    print("obj/opt is the additive objective over the certified optimum; worst is the worst "
          "replicate; switch, flips (mean count), N50 and QAN50 (kb) come from the phasing of "
          "the selected reads")

    units = [(p, c, r, args.variants, args.depth, args.heterozygosity, args.quality_spread,
              errors)
             for p, c, r in itertools.product(args.profiles, args.caps, range(args.replicates))]
    with Pool(args.processes) as pool:
        results = pool.map(unit, units, chunksize=1)

    header = (f"{'profile':>7} {'cap':>4} {'selector':>18} {'obj/opt':>8} {'worst':>7} "
              f"{'MEC full':>10} {'switch':>8} {'sw worst':>9} {'flips':>7} {'N50':>7} "
              f"{'QAN50':>7} {'phased':>7} {'sel s':>7} "
              f"{'dp s':>7} {'mixed':>6}")
    print(header)
    print("-" * len(header))
    for profile, cap in itertools.product(args.profiles, args.caps):
        group = [r for r in results if r[0] == profile and r[1] == cap]
        for name in ALL:
            ratios = [r[4][name]["objective"] / r[4]["optimal"]["objective"] for r in group]
            switches = [r[4][name]["switch"] for r in group]
            ratio_text = (f"{'n/a':>8} {'n/a':>7}" if name in ("hapchat", "merge_optimal") else
                          f"{statistics.mean(ratios):>8.4f} {min(ratios):>7.4f}")
            print(
                f"{profile:>7} {cap:>4} {name:>18} {ratio_text} "
                f"{statistics.mean(r[4][name]['mec_full'] for r in group):>10.1f} "
                f"{statistics.mean(switches):>8.5f} {max(switches):>9.5f} "
                f"{statistics.mean(r[4][name]['flips'] for r in group):>7.1f} "
                f"{statistics.mean(r[4][name]['n50'] for r in group) / 1000:>7.0f} "
                f"{statistics.mean(r[4][name]['qan50'] for r in group) / 1000:>7.0f} "
                f"{statistics.mean(r[4][name]['phased'] for r in group):>7.3f} "
                f"{statistics.mean(r[4][name]['select_ms'] for r in group) / 1000:>7.2f} "
                f"{statistics.mean(r[4][name]['solve_ms'] for r in group) / 1000:>7.2f} "
                + (f"{statistics.mean(r[4][name]['mixed'] for r in group):>6.3f}"
                   if name in ("hapchat", "merge_optimal") else f"{'':>6}")
            )
        print()
    reads = [r[3] for r in results]
    print(f"reads per instance: mean {statistics.mean(reads):.0f}, range {min(reads)} to {max(reads)}")


if __name__ == "__main__":
    main()
