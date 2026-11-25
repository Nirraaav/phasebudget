"""What does assuming every called site is heterozygous cost? This is OQ-7.

`phasebudget.dp` takes the second haplotype to be the complement of the first,
so it cannot represent a site that was called heterozygous and is not. The
general MEC model can. WhatsHap does not need the restriction and says so, so
the restriction is this project's own and its price has to be measured.

Both models are scored on the read bipartition, which is the only object both
produce, and on phase relations between sites that really are heterozygous.
`docs/design/genotype-errors.md` derives both metrics and states what each
outcome would mean.

The general model is solved by `phasebudget.general`, the coverage-
parameterized dynamic program, so the comparison runs at the scale of the other
benchmarks. `--solver bruteforce` switches back to the `2^n` oracle, which is
limited to about sixteen reads and was the only option until 26 September.
"""

from __future__ import annotations

import argparse
import statistics

from phasebudget.dp import CoverageDP
from phasebudget.exact import solve_bipartition_bruteforce
from phasebudget.general import solve_general
from phasebudget.simulate import simulate_with_genotype_errors


def bipartition_error(assignment, origin) -> float:
    """Hamming distance to the true origin, minimized over the global flip."""
    if not origin:
        return float("nan")
    same = sum(1 for a, b in zip(assignment, origin) if a != b)
    return min(same, len(origin) - same) / len(origin)


def general_phase(matrix, assignment):
    """`(haplotype of part 0, which variants the two parts disagree at)`.

    The general model is free to give both parts the same consensus at a
    variant, and when it does it has declined to phase that variant: the column
    carries an allele and no phase. `decided[j]` records whether it kept the two
    parts distinct, which is the general-model analogue of a variant being
    inside a phase block.
    """
    haplotype = []
    decided = []
    for position in range(matrix.n_variants):
        totals = [[0.0, 0.0], [0.0, 0.0]]
        for index in matrix.active(position):
            allele = matrix[index].allele_at(position)
            if allele is None:
                continue
            totals[assignment[index]][allele] += matrix[index].weight_at(position)
        consensus = [0 if side[0] >= side[1] else 1 for side in totals]
        covered = [side[0] + side[1] > 0.0 for side in totals]
        haplotype.append(consensus[0])
        decided.append(all(covered) and consensus[0] != consensus[1])
    return tuple(haplotype), tuple(decided)


def switch_error_on_sites(truth, predicted, sites, decided=None):
    """`(switch error over the scored pairs, fraction of pairs scored)`.

    With `decided` given, a pair is scored only when both of its sites were
    phased, so a model is not charged a switch error at a site where it said
    there is no phase. The complementary model phases every site by
    construction and scores 1.0 on the second value, which is the number that
    has to be read next to the first.
    """
    pairs = list(zip(sites, sites[1:]))
    if not pairs:
        return float("nan"), float("nan")
    scored = [
        (left, right)
        for left, right in pairs
        if decided is None or (decided[left] and decided[right])
    ]
    if not scored:
        return float("nan"), 0.0
    wrong = sum(
        1
        for left, right in scored
        if (predicted[left] ^ predicted[right]) != (truth[left] ^ truth[right])
    )
    return wrong / len(scored), len(scored) / len(pairs)


def run(homozygous_rates, error_rates, replicates, n_variants, coverage, span,
        solver="dp"):
    """One row per (homozygous rate, error rate)."""
    rows = []
    for homozygous_rate in homozygous_rates:
        for error_rate in error_rates:
            complementary_bip: list[float] = []
            general_bip: list[float] = []
            complementary_switch: list[float] = []
            general_switch: list[float] = []
            complementary_cost: list[float] = []
            general_cost: list[float] = []
            general_scored: list[float] = []
            general_phased: list[float] = []
            miscalled: list[float] = []
            reads: list[float] = []
            used = 0
            for replicate in range(replicates):
                instance = simulate_with_genotype_errors(
                    n_variants,
                    coverage=coverage,
                    read_span=span,
                    error_rate=error_rate,
                    homozygous_rate=homozygous_rate,
                    seed=5500 + 29 * replicate,
                )
                matrix = instance.matrix
                sites = instance.heterozygous_sites()
                if len(matrix) < 2 or len(sites) < 3:
                    continue
                used += 1
                reads.append(float(len(matrix)))
                miscalled.append(float(len(instance.homozygous)))

                solution = CoverageDP(matrix).solve()
                complementary_cost.append(solution.cost)
                complementary_bip.append(
                    bipartition_error(solution.assignment, instance.origin)
                )
                complementary_switch.append(
                    switch_error_on_sites(
                        instance.haplotype, solution.haplotype, sites
                    )[0]
                )

                if solver == "bruteforce":
                    cost, assignment = solve_bipartition_bruteforce(matrix)
                else:
                    cost, assignment = solve_general(matrix)
                general_cost.append(cost)
                general_bip.append(bipartition_error(assignment, instance.origin))
                predicted, decided = general_phase(matrix, assignment)
                general_switch.append(
                    switch_error_on_sites(instance.haplotype, predicted, sites)[0]
                )
                scored, fraction = switch_error_on_sites(
                    instance.haplotype, predicted, sites, decided
                )
                if scored == scored:
                    general_scored.append(scored)
                general_phased.append(fraction)
            rows.append(
                {
                    "homozygous_rate": homozygous_rate,
                    "error_rate": error_rate,
                    "reads": statistics.fmean(reads) if reads else float("nan"),
                    "miscalled": (
                        statistics.fmean(miscalled) if miscalled else float("nan")
                    ),
                    "complementary_cost": (
                        statistics.fmean(complementary_cost)
                        if complementary_cost
                        else float("nan")
                    ),
                    "general_cost": (
                        statistics.fmean(general_cost) if general_cost else float("nan")
                    ),
                    "complementary_bip": (
                        statistics.fmean(complementary_bip)
                        if complementary_bip
                        else float("nan")
                    ),
                    "general_bip": (
                        statistics.fmean(general_bip) if general_bip else float("nan")
                    ),
                    "complementary_switch": (
                        statistics.fmean(complementary_switch)
                        if complementary_switch
                        else float("nan")
                    ),
                    "general_switch": (
                        statistics.fmean(general_switch)
                        if general_switch
                        else float("nan")
                    ),
                    "general_scored": (
                        statistics.fmean(general_scored)
                        if general_scored
                        else float("nan")
                    ),
                    "general_phased": (
                        statistics.fmean(general_phased)
                        if general_phased
                        else float("nan")
                    ),
                    "instances": used,
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=int, default=12)
    parser.add_argument("--coverage", type=float, default=6.0)
    parser.add_argument("--span", type=int, default=5)
    parser.add_argument("--replicates", type=int, default=40)
    parser.add_argument(
        "--homozygous-rates", type=float, nargs="+", default=[0.0, 0.1, 0.2, 0.3]
    )
    parser.add_argument("--error-rates", type=float, nargs="+", default=[0.02, 0.10])
    parser.add_argument("--solver", choices=["dp", "bruteforce"], default="dp")
    args = parser.parse_args()

    rows = run(
        args.homozygous_rates,
        args.error_rates,
        args.replicates,
        args.variants,
        args.coverage,
        args.span,
        solver=args.solver,
    )
    print(
        f"{args.variants} variants, coverage {args.coverage:g}, span {args.span}, "
        f"{args.replicates} replicates, seed 5500 + 29 * replicate, "
        f"general model solved by {args.solver}; simulated data"
    )
    print(
        "bip is read-bipartition error against the true origin, flip-invariant; "
        "switch is over consecutive truly heterozygous sites"
    )
    print(
        "sw gen charges the general model at sites it declined to phase; sw gen* "
        "scores only the pairs it phased, and gen phas is that fraction"
    )
    print("the complementary model phases every site, so its fraction is 1.000")
    print("the general model is a relaxation, so its cost is a lower bound")
    print("simulated data; the general model is exact under either solver")
    header = (
        f"{'hom':>5} {'err':>6} {'miscall':>8} {'cost cmp':>9} "
        f"{'cost gen':>9} {'bip cmp':>8} {'bip gen':>8} {'sw cmp':>8} "
        f"{'sw gen':>8} {'sw gen*':>8} {'gen phas':>9} {'n':>4}"
    )
    print()
    print(header)
    print("-" * len(header))
    previous = None
    for row in rows:
        if previous is not None and row["homozygous_rate"] != previous:
            print()
        previous = row["homozygous_rate"]
        print(
            f"{row['homozygous_rate']:>5.2f} {row['error_rate']:>6.2f} "
            f"{row['miscalled']:>8.2f} "
            f"{row['complementary_cost']:>9.3f} {row['general_cost']:>9.3f} "
            f"{row['complementary_bip']:>8.4f} {row['general_bip']:>8.4f} "
            f"{row['complementary_switch']:>8.4f} {row['general_switch']:>8.4f} "
            f"{row['general_scored']:>8.4f} {row['general_phased']:>9.3f} "
            f"{row['instances']:>4}"
        )


if __name__ == "__main__":
    main()
