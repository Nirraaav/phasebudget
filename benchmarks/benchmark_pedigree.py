"""Read selection in pedigree mode: a trio shares WhatsHap's coverage budget.

`whatshap phase` with a pedigree gives each member `max(1, cap // members)`
reads of coverage, so a trio at the default cap of 15 selects each member's
reads at cap 5, which is where selection changed accuracy on single samples.
This asks whether it still does when transmission from the parents also
carries phase.

The trio is simulated here, since the library's simulators make one diploid.
Four parental haplotypes draw alleles independently with probability 1/2 at
sites spaced by the genomic simulator's geometric gaps; sites where no member
is heterozygous are dropped. Each parent transmits one haplotype to the child,
recombining at the uniform rate WhatsHap assumes (1.26 cM/Mb). Each member
gets reads of one length at a given depth, each read from one of its two
haplotypes, with per-base phred qualities jittered around the profile's error
rate, and each allele flipped with the probability its quality states. A
read's integer quality at a site is both WhatsHap's quality and the selectors'
weight.

Each member's reads are selected at the per-member cap by one selector, and
the trio is phased by WhatsHap 2.8's `PedigreeDPTable` with trusted genotypes
and the uniform recombination map, as `whatshap phase --ped` does. The child's
switch error is counted over consecutive child-heterozygous sites as one
block: with genotypes trusted, WhatsHap's own components join every site where
some member is homozygous, so the whole stretch is one block in its output.

As a reference the child alone is phased by the same DP at the same cap,
where switch error is counted inside blocks of the selected reads.

WhatsHap is an optional dependency (`pip install -e ".[whatshap]"`).
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import math
import os
import platform
import random
import statistics
import subprocess
from multiprocessing import Pool
from pathlib import Path

from phasebudget._types import Fragment
from phasebudget.fragments import FragmentMatrix
from phasebudget.hapchat import merge_reads
from phasebudget.selection import read_weights, select_optimal, select_random, select_whatshap
from phasebudget.simulate import geometric_gap

PROFILES = {
    "hifi": {"read_length": 15_000, "error_rate": 0.01},
    "ont": {"read_length": 10_000, "error_rate": 0.05},
    "noisy": {"read_length": 10_000, "error_rate": 0.10},
}
SELECTORS = ("random", "whatshap", "optimal", "merge_optimal")
MEMBERS = ("father", "mother", "child")
RECOMBINATION_RATE = 1.26  # cM/Mb, WhatsHap's default


def simulate_trio(n_sites, *, depth, read_length, error_rate, spread, seed,
                  heterozygosity=1e-3):
    rng = random.Random(seed)
    haplotypes = [[] for _ in range(4)]  # father 0, 1; mother 2, 3
    positions = []
    position = 1
    while len(positions) < n_sites:
        position += geometric_gap(rng, heterozygosity)
        alleles = [rng.randint(0, 1) for _ in range(4)]
        if len(set(alleles)) == 1 and rng.random() < 0.5:
            continue
        positions.append(position)
        for h, allele in enumerate(alleles):
            haplotypes[h].append(allele)
    # Transmission with crossovers at the uniform rate.
    transmitted = []
    for parent in (0, 1):
        current = rng.randint(0, 1)
        choice = []
        for j in range(n_sites):
            if j:
                morgans = (positions[j] - positions[j - 1]) * 1e-8 * RECOMBINATION_RATE
                if rng.random() < morgans:
                    current = 1 - current
            choice.append(current)
        transmitted.append(choice)
    child = [
        [haplotypes[c][j] for j, c in enumerate(transmitted[0])],
        [haplotypes[2 + c][j] for j, c in enumerate(transmitted[1])],
    ]
    genomes = {"father": haplotypes[0:2], "mother": haplotypes[2:4], "child": child}
    keep = [j for j in range(n_sites)
            if any(genomes[m][0][j] != genomes[m][1][j] for m in MEMBERS)]
    positions = [positions[j] for j in keep]
    genomes = {m: [[h[j] for j in keep] for h in pair] for m, pair in genomes.items()}
    reads = {m: _reads(rng, genomes[m], positions, depth, read_length, error_rate, spread)
             for m in MEMBERS}
    return positions, genomes, reads


def _reads(rng, pair, positions, depth, read_length, error_rate, spread):
    lowest = positions[0] - read_length + 1
    span = positions[-1] - lowest + 1
    n_reads = round(depth * span / read_length)
    base_phred = -10 * math.log10(error_rate)
    fragments = []
    for read_id in range(n_reads):
        start = rng.randint(lowest, positions[-1])
        hap = pair[rng.randint(0, 1)]
        calls = []
        for j, position in enumerate(positions):
            if position < start:
                continue
            if position >= start + read_length:
                break
            phred = max(1.0, base_phred + rng.gauss(0.0, spread))
            quality = int(round(phred))
            wrong = rng.random() < 10 ** (-quality / 10)
            calls.append((j, hap[j] ^ int(wrong), float(quality)))
        if not calls:
            continue
        first = calls[0][0]
        alleles = tuple(a for _, a, _ in calls)
        weights = tuple(w for _, _, w in calls)
        assert [j for j, _, _ in calls] == list(range(first, first + len(calls)))
        fragments.append(Fragment(len(fragments), first, alleles, weights))
    return FragmentMatrix(fragments, len(positions))


def _select(name, matrix, cap, positions, seed):
    if name == "merge_optimal":
        merged = merge_reads(matrix).matrix
        chosen = select_optimal(merged, cap, read_weights(merged)).selected
        return merged.subset(chosen)
    weights = read_weights(matrix)
    if name == "random":
        chosen = select_random(matrix, cap, seed=seed, weights=weights).selected
    elif name == "whatshap":
        chosen = select_whatshap(matrix, cap, weights, positions=positions).selected
    else:
        chosen = select_optimal(matrix, cap, weights).selected
    return matrix.subset(chosen)


def _add_reads(readset, matrix, positions, sample_id, prefix):
    from whatshap.core import Read

    for fragment in matrix:
        if sum(1 for a in fragment.alleles if a is not None) < 2:
            continue
        read = Read(f"{prefix}{fragment.read_id}", 60, 0, sample_id)
        for j, allele, weight in fragment.observed():
            read.add_variant(position=positions[j], allele=int(allele),
                             quality=max(1, int(round(weight))))
        read.sort()
        readset.add(read)


def _phase(selected, genomes, positions, members):
    from whatshap.core import (Genotype, NumericSampleIds, Pedigree, PedigreeDPTable,
                               ReadSet)
    from whatshap.pedigree import UniformRecombinationCostComputer

    ids = NumericSampleIds()
    pedigree = Pedigree(ids)
    readset = ReadSet()
    for member in members:
        pair = genomes[member]
        pedigree.add_individual(member, [Genotype(sorted([a, b])) for a, b in zip(*pair)])
        _add_reads(readset, selected[member], positions, ids[member], member[0])
    if len(members) == 3:
        pedigree.add_relationship(father_id="father", mother_id="mother", child_id="child")
    readset.sort()
    costs = UniformRecombinationCostComputer.uniform_recombination_map(
        RECOMBINATION_RATE, positions)
    table = PedigreeDPTable(readset, costs, pedigree, False, positions)
    superreads, _ = table.get_super_reads()
    child = superreads[members.index("child")]
    alleles = {variant.position: variant.allele for variant in child[0]}
    return [alleles.get(p) for p in positions]


def _switch_rate(truth, predicted, blocks=None):
    het = [j for j in range(len(truth[0])) if truth[0][j] != truth[1][j]
           and predicted[j] is not None]
    switches = pairs = 0
    for left, right in zip(het, het[1:]):
        if blocks is not None and (blocks[left] is None or blocks[left] != blocks[right]):
            continue
        pairs += 1
        if (predicted[left] ^ predicted[right]) != (truth[0][left] ^ truth[0][right]):
            switches += 1
    return switches, pairs


def unit(arguments):
    profile, cap, replicate, n_sites, depth, spread, seed_base = arguments
    positions, genomes, reads = simulate_trio(
        n_sites, depth=depth, spread=spread, seed=seed_base + 17 * replicate,
        **PROFILES[profile])
    out = {}
    for name in SELECTORS:
        selected = {m: _select(name, reads[m], cap, positions, replicate) for m in MEMBERS}
        trio = _switch_rate(genomes["child"], _phase(selected, genomes, positions, MEMBERS))
        child = selected["child"]
        alone = _switch_rate(genomes["child"], _phase(selected, genomes, positions, ("child",)),
                             blocks=child.blocks())
        out[name] = {"trio": trio, "alone": alone}
    return profile, cap, out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", nargs="+", default=["ont", "noisy"])
    parser.add_argument("--caps", type=int, nargs="+", default=[3, 5],
                        help="per-member caps; 5 is a trio at WhatsHap's default of 15")
    parser.add_argument("--sites", type=int, default=2000)
    parser.add_argument("--depth", type=float, default=30.0, help="per member")
    parser.add_argument("--quality-spread", type=float, default=2.0)
    parser.add_argument("--replicates", type=int, default=16)
    parser.add_argument("--seed-base", type=int, default=30000)
    parser.add_argument("--processes", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    args = parser.parse_args()
    import whatshap

    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    print(f"commit {commit}; {platform.system()} {platform.release()} {platform.machine()}; "
          f"Python {platform.python_version()}; WhatsHap {whatshap.__version__}")
    print(f"trio, {args.sites} segregating sites at 1e-3, depth {args.depth:g} per member, "
          f"quality spread {args.quality_spread:g}, {args.replicates} replicates, seed "
          f"{args.seed_base} + 17 * replicate, recombination {RECOMBINATION_RATE} cM/Mb; "
          f"profiles " + ", ".join(f"{p} ({PROFILES[p]['read_length']} bp, error "
                                   f"{PROFILES[p]['error_rate']})" for p in args.profiles)
          + "; simulated data")
    print("switch error of the child, pooled over replicates (switches / pairs), with the worst "
          "replicate; trio = PedigreeDPTable on all three members at the per-member cap, "
          "alone = the child by itself at the same cap, counted inside read blocks")
    units = list(itertools.product(args.profiles, args.caps, range(args.replicates),
                                   [args.sites], [args.depth], [args.quality_spread],
                                   [args.seed_base]))
    with Pool(args.processes) as pool:
        results = pool.map(unit, units, chunksize=1)
    header = (f"{'profile':>7} {'cap':>4} {'selector':>14} {'trio':>9} {'worst':>8} "
              f"{'alone':>9} {'worst':>8} {'trio pairs':>10}")
    print(header)
    print("-" * len(header))
    for profile, cap in itertools.product(args.profiles, args.caps):
        group = [r[2] for r in results if r[0] == profile and r[1] == cap]
        for name in SELECTORS:
            cells = []
            for mode in ("trio", "alone"):
                counts = [g[name][mode] for g in group]
                pooled = sum(s for s, _ in counts) / max(1, sum(p for _, p in counts))
                worst = max(s / p if p else 0.0 for s, p in counts)
                cells.append((pooled, worst))
            pairs = statistics.mean(g[name]["trio"][1] for g in group)
            print(f"{profile:>7} {cap:>4} {name:>14} {cells[0][0]:>9.5f} {cells[0][1]:>8.5f} "
                  f"{cells[1][0]:>9.5f} {cells[1][1]:>8.5f} {pairs:>10.0f}")
        print()


if __name__ == "__main__":
    main()
