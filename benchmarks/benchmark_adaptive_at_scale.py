"""The adaptive loop at chromosome scale, paired against exact static selection.

`benchmark_adaptive_signal.py` measures the loop on 60-variant instances.
This runs it on the chromosome-scale instances of
`benchmark_selection_at_scale.py`: exact static selection and the loop
steered by MEC margins, on the same instances, switch error paired per
instance with a 95 per cent interval. Negative means the loop phases better.
A third profile, `noisy`, is 10 kb at error 0.10, where switch errors are
frequent enough at chromosome scale to compare.
The loop's boost, quantile and iteration cap are the library defaults.
"""

from __future__ import annotations

import argparse
import itertools
import math
import os
import platform
import subprocess
from multiprocessing import Pool

from phasebudget.adaptive import adaptive_selection, static_selection
from phasebudget.hapchat import merge_reads
from phasebudget.metrics import phasing_from_solution, summarize
from phasebudget.simulate import simulate_genomic

PROFILES = {
    "hifi": {"read_length": 15_000, "error_rate": 0.01},
    "ont": {"read_length": 10_000, "error_rate": 0.05},
    "noisy": {"read_length": 10_000, "error_rate": 0.10},
}


def unit(arguments):
    (profile, cap, replicate, n_variants, depth, heterozygosity, spread, signal,
     probability, acceptance, window, control, iterations, boost, merge, seed_base,
     errors, merge_options) = arguments
    instance = simulate_genomic(
        n_variants, coverage=depth, heterozygosity=heterozygosity, quality_spread=spread,
        seed=seed_base + 17 * replicate, **PROFILES[profile], **errors,
    )
    matrix = merge_reads(instance.matrix, **merge_options).matrix if merge else instance.matrix
    loop = dict(iterations=iterations, boost=boost)
    rates = {}
    blocks = {}
    arms = [
        ("exact", static_selection(matrix, cap)),
        ("adaptive", adaptive_selection(matrix, cap, signal=signal,
                                        error_probability=probability,
                                        acceptance=acceptance, window=window,
                                        shuffle_seed=replicate, **loop)),
    ]
    if control:
        arms.append(("control", adaptive_selection(
            matrix, cap, signal="shuffled", error_probability=probability,
            acceptance=acceptance, window=window, shuffle_seed=replicate, **loop)))
    for name, result in arms:
        subset = matrix.subset(result.support or result.selected)
        report = summarize(instance.haplotype, phasing_from_solution(subset, result.solution),
                           instance.positions)
        rates[name] = report["switch_error_rate"]
        blocks[name] = (report["block_n50"], report["qan50"])
    return profile, cap, rates["exact"], rates["adaptive"], rates.get("control"), blocks


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", nargs="+", default=list(PROFILES))
    parser.add_argument("--variants", type=int, default=2000)
    parser.add_argument("--depth", type=float, default=30.0)
    parser.add_argument("--heterozygosity", type=float, default=0.001)
    parser.add_argument("--quality-spread", type=float, default=2.0)
    parser.add_argument("--caps", type=int, nargs="+", default=[3, 5, 10])
    parser.add_argument("--replicates", type=int, default=16)
    parser.add_argument("--signal", default="margin")
    parser.add_argument("--error-probability", type=float, default=None,
                        help="absolute fragile rule; needs --signal posterior")
    parser.add_argument("--acceptance", default="global", choices=["global", "stitched", "stitched_oriented"])
    parser.add_argument("--window", type=int, default=25)
    parser.add_argument("--seed-base", type=int, default=8000,
                        help="instances use seed-base + 17 * replicate; change it for held-out runs")
    parser.add_argument("--systematic-rate", type=float, default=0.0)
    parser.add_argument("--systematic-strength", type=float, default=0.5)
    parser.add_argument("--burst-rate", type=float, default=0.0)
    parser.add_argument("--burst-length", type=int, default=5)
    parser.add_argument("--burst-error", type=float, default=0.35)
    parser.add_argument("--iterations", type=int, default=4)
    parser.add_argument("--boost", type=float, default=2.0)
    parser.add_argument("--merge", action="store_true",
                        help="merge reads by HapCHAT's rule first; both arms then run on merged reads")
    parser.add_argument("--merge-p", type=float, default=0.15,
                        help="HapCHAT's error probability in the merge rule")
    parser.add_argument("--merge-threshold", type=float, default=1e6,
                        help="HapCHAT's likelihood-ratio threshold for merging")
    parser.add_argument("--control", action="store_true",
                        help="also run the shuffled-signal control and pair the signal against it")
    parser.add_argument("--processes", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    args = parser.parse_args()
    errors = {"systematic_rate": args.systematic_rate,
              "systematic_strength": args.systematic_strength, "burst_rate": args.burst_rate,
              "burst_length": args.burst_length, "burst_error": args.burst_error}
    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    print(f"commit {commit}; {platform.system()} {platform.release()} {platform.machine()}; "
          f"Python {platform.python_version()}")
    print(f"{args.variants} variants, depth {args.depth:g}, heterozygosity {args.heterozygosity:g}, "
          f"quality spread {args.quality_spread:g}, {args.replicates} replicates, seed "
          f"{args.seed_base} + 17 * replicate, signal {args.signal}, fragile rule "
          f"{'P(wrong) > ' + str(args.error_probability) if args.error_probability else 'quantile 0.25'}"
          f", acceptance {args.acceptance}"
          + (f" with window {args.window}" if args.acceptance == "stitched" else "")
          + f", {args.iterations} iterations, boost {args.boost:g}"
          + (f", reads merged by HapCHAT's rule first (p {args.merge_p:g}, threshold "
             f"{args.merge_threshold:g})" if args.merge else "")
          + "; unreported errors " + ", ".join(f"{k} {v:g}" for k, v in errors.items())
          + "; simulated data")
    units = [(p, c, r, args.variants, args.depth, args.heterozygosity, args.quality_spread,
              args.signal, args.error_probability, args.acceptance, args.window, args.control,
              args.iterations, args.boost, args.merge, args.seed_base, errors,
              {"p": args.merge_p, "threshold": args.merge_threshold})
             for p, c, r in itertools.product(args.profiles, args.caps, range(args.replicates))]
    with Pool(args.processes) as pool:
        results = pool.map(unit, units, chunksize=1)
    print("switch error, paired; N50 and QAN50 are means in kb, exact arm then adaptive arm")
    header = (f"{'profile':>7} {'cap':>4} {'exact':>9} {'adaptive':>9} {'diff':>9} "
              f"{'95% interval':>21} {'wins':>7} {'N50 ex':>7} {'N50 ad':>7} "
              f"{'QAN50 ex':>8} {'QAN50 ad':>8}")
    print(header)
    print("-" * len(header))
    for profile, cap in itertools.product(args.profiles, args.caps):
        group = [r for r in results if r[0] == profile and r[1] == cap]
        pairs = [(r[2], r[3]) for r in group]
        kb = lambda arm, i: sum(r[5][arm][i] for r in group) / len(group) / 1000
        diffs = [b - a for a, b in pairs]
        n = len(diffs)
        mean = sum(diffs) / n
        se = math.sqrt(sum((d - mean) ** 2 for d in diffs) / (n - 1) / n) if n > 1 else 0.0
        print(f"{profile:>7} {cap:>4} {sum(a for a, _ in pairs) / n:>9.5f} "
              f"{sum(b for _, b in pairs) / n:>9.5f} {mean:>+9.5f} "
              f"[{mean - 1.96 * se:+.5f}, {mean + 1.96 * se:+.5f}] "
              f"{sum(d < 0 for d in diffs):>3}:{sum(d > 0 for d in diffs):<3} "
              f"{kb('exact', 0):>7.0f} {kb('adaptive', 0):>7.0f} "
              f"{kb('exact', 1):>8.0f} {kb('adaptive', 1):>8.0f}")

    if args.control:
        print()
        print("signal minus shuffled control, paired; negative means the signal steers better")
        for profile, cap in itertools.product(args.profiles, args.caps):
            diffs = [r[3] - r[4] for r in results if r[0] == profile and r[1] == cap]
            n = len(diffs)
            mean = sum(diffs) / n
            se = math.sqrt(sum((d - mean) ** 2 for d in diffs) / (n - 1) / n) if n > 1 else 0.0
            print(f"{profile:>7} {cap:>4} {mean:>+9.5f} [{mean - 1.96 * se:+.5f}, "
                  f"{mean + 1.96 * se:+.5f}] {sum(d < 0 for d in diffs):>3}:"
                  f"{sum(d > 0 for d in diffs):<3}")


if __name__ == "__main__":
    main()
