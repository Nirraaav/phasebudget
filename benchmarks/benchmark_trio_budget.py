"""How a trio should spend WhatsHap's coverage budget (research plan track T3).

WhatsHap's pedigree DP tracks every active read of every member, so its state
space grows with the total number of selected reads covering a column across
the trio. `whatshap phase` splits a budget `C` into `C // 3` per member. This
compares three allocations of the same total `C`:

* `thirds`: each member's reads selected exactly at `C // 3` (WhatsHap's split);
* `shared`: one exact selection over the union of the three members' reads,
  under a single cap `C` on their summed span coverage;
* `shared_weighted`: the same, with each allele's selection weight 1 where
  the child is heterozygous and both parents are heterozygous (the sites
  transmission cannot phase) and 0.1 elsewhere.

The union of interval matrices is an interval matrix, so `shared` is still
exact. Phasing always uses the original qualities, by WhatsHap's
`PedigreeDPTable`; the child's switch error, the DP's wall time and the
largest column state are reported.

Decision rule, fixed in the research plan (development repository): the shared budget is useful
if `shared_weighted` at a total of 9 or less matches the child's switch error
of `thirds` at 15 (paired 95 per cent interval of the difference containing
zero or below it), which would make the DP state space at least 64 times
smaller.
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import math
import os
import platform
import statistics
import subprocess
import time
from multiprocessing import Pool
from pathlib import Path

from phasebudget._types import Fragment
from phasebudget.fragments import FragmentMatrix
from phasebudget.selection import read_weights, select_optimal

ALLOCATIONS = ("thirds", "shared", "shared_weighted")


def _pedigree():
    path = Path(__file__).with_name("benchmark_pedigree.py")
    spec = importlib.util.spec_from_file_location("benchmark_pedigree", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _union(reads, members, weight_of=None):
    """One matrix of every member's reads, with a map back to (member, read id)."""
    fragments, origin = [], []
    for member in members:
        for f in reads[member]:
            weights = f.weights if weight_of is None else tuple(
                0.0 if a is None else weight_of(f.start + i)
                for i, a in enumerate(f.alleles))
            fragments.append(Fragment(len(fragments), f.start, f.alleles, weights))
            origin.append((member, f.read_id))
    n = next(iter(reads.values())).n_variants
    return FragmentMatrix(fragments, n), origin


def select(allocation, reads, genomes, total, members):
    if allocation == "thirds":
        per = max(1, total // len(members))
        return {m: reads[m].subset(select_optimal(reads[m], per, read_weights(reads[m])).selected)
                for m in members}
    weight_of = None
    if allocation == "shared_weighted":
        het = lambda m, j: genomes[m][0][j] != genomes[m][1][j]
        ambiguous = [het("child", j) and het("father", j) and het("mother", j)
                     for j in range(len(genomes["child"][0]))]
        weight_of = lambda j: 1.0 if ambiguous[j] else 0.1
    union, origin = _union(reads, members, weight_of)
    chosen = select_optimal(union, total, read_weights(union)).selected
    picked = {m: [] for m in members}
    for index in chosen:
        member, read_id = origin[index]
        picked[member].append(read_id)
    return {m: reads[m].subset(picked[m]) for m in members}


def unit(arguments):
    profile, total, replicate, n_sites, depth, spread, seed_base = arguments
    pedigree = _pedigree()
    positions, genomes, reads = pedigree.simulate_trio(
        n_sites, depth=depth, spread=spread, seed=seed_base + 17 * replicate,
        **pedigree.PROFILES[profile])
    out = {}
    for allocation in ALLOCATIONS:
        selected = select(allocation, reads, genomes, total, pedigree.MEMBERS)
        started = time.perf_counter()
        predicted = pedigree._phase(selected, genomes, positions, pedigree.MEMBERS)
        seconds = time.perf_counter() - started
        switches, pairs = pedigree._switch_rate(genomes["child"], predicted)
        coverage = [0] * len(positions)
        for matrix in selected.values():
            for j, c in enumerate(matrix.coverage()):
                coverage[j] += c
        out[allocation] = (switches, pairs, seconds, max(coverage, default=0))
    return profile, total, out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", nargs="+", default=["ont", "noisy"])
    parser.add_argument("--totals", type=int, nargs="+", default=[6, 9, 12, 15])
    parser.add_argument("--sites", type=int, default=2000)
    parser.add_argument("--depth", type=float, default=30.0)
    parser.add_argument("--quality-spread", type=float, default=2.0)
    parser.add_argument("--replicates", type=int, default=16)
    parser.add_argument("--seed-base", type=int, default=30000)
    parser.add_argument("--processes", type=int, default=max(1, (os.cpu_count() or 2) - 4))
    args = parser.parse_args()
    import whatshap

    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    print(f"commit {commit}; {platform.system()} {platform.release()} {platform.machine()}; "
          f"Python {platform.python_version()}; WhatsHap {whatshap.__version__}")
    print(f"trio, {args.sites} segregating sites, depth {args.depth:g} per member, quality "
          f"spread {args.quality_spread:g}, {args.replicates} replicates, seed {args.seed_base} + "
          "17 * replicate; simulated data")
    print("child switch error pooled over replicates, worst replicate, mean DP seconds, mean "
          "largest summed coverage; last column: shared_weighted minus thirds at total 15, paired")
    units = list(itertools.product(args.profiles, args.totals, range(args.replicates),
                                   [args.sites], [args.depth], [args.quality_spread],
                                   [args.seed_base]))
    with Pool(args.processes) as pool:
        results = pool.map(unit, units, chunksize=1)
    header = (f"{'profile':>7} {'total':>6} {'allocation':>16} {'switch':>9} {'worst':>8} "
              f"{'DP s':>7} {'max cov':>8} {'vs thirds@max [95%]':>34}")
    print(header)
    print("-" * len(header))
    for profile, total in itertools.product(args.profiles, args.totals):
        group = [r[2] for r in results if r[0] == profile and r[1] == total]
        base = [r[2]["thirds"] for r in results
                if r[0] == profile and r[1] == max(args.totals)]
        for allocation in ALLOCATIONS:
            rows = [g[allocation] for g in group]
            pooled = sum(r[0] for r in rows) / max(1, sum(r[1] for r in rows))
            worst = max(r[0] / r[1] if r[1] else 0.0 for r in rows)
            diffs = [a[0] / a[1] - b[0] / b[1] for a, b in zip(rows, base) if a[1] and b[1]]
            mean = sum(diffs) / len(diffs)
            se = (math.sqrt(sum((d - mean) ** 2 for d in diffs) / (len(diffs) - 1) / len(diffs))
                  if len(diffs) > 1 else 0.0)
            print(f"{profile:>7} {total:>6} {allocation:>16} {pooled:>9.5f} {worst:>8.5f} "
                  f"{statistics.mean(r[2] for r in rows):>7.2f} "
                  f"{statistics.mean(r[3] for r in rows):>8.1f} "
                  f"{mean:>+10.5f} [{mean - 1.96 * se:+.5f}, {mean + 1.96 * se:+.5f}]")
        print()


if __name__ == "__main__":
    main()
