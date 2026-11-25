"""What WhatsHap's score update costs: the port against the paper's rule.

`selection.whatshap_score_drop` reproduces WhatsHap 2.8, which lowers a
queued read's score by the number of its variants the admitted read did not
newly cover. The paper's `scoredyn` lowers it by the number that were newly
covered. This diagnostic runs the port twice, once as shipped and once with
the paper's rule substituted and nothing else changed, and phases each
selection with the exact DP. It is a diagnostic of one line of WhatsHap, not
a selector the library offers.
"""

from __future__ import annotations

import argparse
import inspect
import itertools
import platform
import statistics
import subprocess
from multiprocessing import Pool

import phasebudget.selection as selection
from phasebudget.dp import CoverageDP
from phasebudget.metrics import phasing_from_solution, summarize
from phasebudget.simulate import simulate_genomic

PROFILES = {
    "ont": {"read_length": 10_000, "error_rate": 0.05},
    "noisy": {"read_length": 10_000, "error_rate": 0.10},
}


def _paper_rule_selector():
    source = inspect.getsource(selection.select_whatshap)
    shipped = "drop = whatshap_score_drop(called[other], fresh)"
    paper = "drop = sum(1 for variant in called[other] if variant in fresh)"
    if shipped not in source:
        raise RuntimeError("select_whatshap no longer has the line this diagnostic replaces")
    namespace = dict(vars(selection))
    exec(source.replace(shipped, paper).replace("def select_whatshap(", "def paper_rule("),
         namespace)
    return namespace["paper_rule"]


def unit(arguments):
    profile, cap, replicate, n_variants = arguments
    instance = simulate_genomic(n_variants, coverage=30.0, heterozygosity=0.001,
                                quality_spread=2.0, seed=8000 + 17 * replicate,
                                **PROFILES[profile])
    matrix = instance.matrix
    weights = selection.read_weights(matrix)
    optimum = selection.select_optimal(matrix, cap, weights).objective
    paper_rule = _paper_rule_selector()
    out = {}
    for name, chosen in (
        ("shipped", selection.select_whatshap(matrix, cap, weights, positions=instance.positions)),
        ("paper rule", paper_rule(matrix, cap, weights, positions=instance.positions)),
    ):
        subset = matrix.subset(chosen.selected)
        report = summarize(instance.haplotype,
                           phasing_from_solution(subset, CoverageDP(subset).solve()),
                           instance.positions)
        out[name] = {
            "ratio": selection.selection_objective(weights, chosen.selected) / optimum,
            "blocks": selection.component_count(matrix, chosen.selected),
            "switch": report["switch_error_rate"],
            "n50": report["block_n50"] / 1000,
            "qan50": report["qan50"] / 1000,
        }
    return profile, cap, out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", nargs="+", default=["ont", "noisy"])
    parser.add_argument("--caps", type=int, nargs="+", default=[3, 5, 10, 15])
    parser.add_argument("--variants", type=int, default=2000)
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--processes", type=int, default=8)
    args = parser.parse_args()
    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    print(f"commit {commit}; {platform.system()} {platform.release()} {platform.machine()}; "
          f"Python {platform.python_version()}")
    print(f"{args.variants} variants, depth 30, heterozygosity 0.001, quality spread 2, "
          f"{args.replicates} replicates, seed 8000 + 17 * replicate; profiles "
          + ", ".join(f"{p} ({PROFILES[p]['read_length']} bp, error {PROFILES[p]['error_rate']})"
                      for p in args.profiles) + "; simulated data")
    print("ratio is the additive objective over the certified optimum (mean, worst); blocks, "
          "switch (mean, worst), N50 and QAN50 (kb) are means over replicates")
    units = list(itertools.product(args.profiles, args.caps, range(args.replicates),
                                   [args.variants]))
    with Pool(args.processes) as pool:
        results = pool.map(unit, units, chunksize=1)
    header = (f"{'profile':>7} {'cap':>4} {'rule':>11} {'ratio':>7} {'worst':>7} {'blocks':>7} "
              f"{'switch':>8} {'sw worst':>9} {'N50':>6} {'QAN50':>6}")
    print(header)
    print("-" * len(header))
    for profile, cap in itertools.product(args.profiles, args.caps):
        group = [r[2] for r in results if r[0] == profile and r[1] == cap]
        for name in ("shipped", "paper rule"):
            rows = [g[name] for g in group]
            print(f"{profile:>7} {cap:>4} {name:>11} "
                  f"{statistics.mean(r['ratio'] for r in rows):>7.4f} "
                  f"{min(r['ratio'] for r in rows):>7.4f} "
                  f"{statistics.mean(r['blocks'] for r in rows):>7.1f} "
                  f"{statistics.mean(r['switch'] for r in rows):>8.5f} "
                  f"{max(r['switch'] for r in rows):>9.5f} "
                  f"{statistics.mean(r['n50'] for r in rows):>6.0f} "
                  f"{statistics.mean(r['qan50'] for r in rows):>6.0f}")
        print()


if __name__ == "__main__":
    main()
