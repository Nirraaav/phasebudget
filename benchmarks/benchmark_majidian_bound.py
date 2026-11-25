"""Is Theorem 2 of Majidian, Kahaei and de Ridder 2020 the probability it claims?

Their Theorem 1 is a sufficient condition for the true haplotype to lose under
unweighted MEC: some column has more than half of its entries wrong. Their
Theorem 2 then computes P{c-MEC}, the probability that MEC recovers the truth,
as the probability that no column meets that condition. The step between the
two, "if the antecedent is not satisfied, the MEC approach works properly", is
asserted and not proved.

The condition is not necessary. `COUNTEREXAMPLE` has no column with half or
more of its entries wrong, and a haplotype other than the truth still has lower
MEC. So Theorem 2 is an upper bound on the probability that the truth is an
optimum, and this script measures how loose it is with the exact dynamic
program under uniform weights, which is their unweighted objective.

Columns per cell:

  * thm2: their formula evaluated on each instance's actual column coverages,
    averaged over replicates.
  * no-bad-col: the fraction of instances in which no column has more than half
    of its entries wrong, counted against each read's true origin. This is the
    event Theorem 2 computes, and agreeing with `thm2` checks the formula.
  * truth-opt: the fraction in which the truth is an optimum of exact MEC.
  * gap-cases: no bad column, and the truth is still not optimal.
"""

from __future__ import annotations

import argparse
from math import comb

from phasebudget import dp
from phasebudget.fragments import FragmentMatrix
from phasebudget.simulate import SimulatedInstance, simulate

#: Five reads over five variants. The truth is `00000`; origins are the side
#: each read was drawn from, 0 for the truth and 1 for its complement. Column
#: errors are (0, 1, 1, 1, 0) at coverages (1, 3, 4, 3, 2), so no column is
#: even half wrong, yet `00010` costs 2 against the truth's 3.
COUNTEREXAMPLE = {
    "rows": ("--10-", "-1101", "--111", "11---", "-10--"),
    "truth": (0, 0, 0, 0, 0),
    "origin": (0, 1, 1, 1, 0),
    "better": (0, 0, 0, 1, 0),
}


def column_errors(
    matrix: FragmentMatrix, haplotype, origin
) -> tuple[list[int], list[int]]:
    """Wrong entries and called entries per column, against true read origins."""
    errors = [0] * matrix.n_variants
    called = [0] * matrix.n_variants
    for fragment, side in zip(matrix.fragments, origin):
        for offset, allele in enumerate(fragment.alleles):
            if allele is None:
                continue
            column = fragment.start + offset
            called[column] += 1
            expected = haplotype[column] if side == 0 else 1 - haplotype[column]
            errors[column] += allele != expected
    return errors, called


def theorem2_probability(called: list[int], error_rate: float) -> float:
    """Majidian et al. equation 7 at the given per-column coverages."""
    probability = 1.0
    for coverage in called:
        probability *= sum(
            comb(coverage, n) * error_rate**n * (1 - error_rate) ** (coverage - n)
            for n in range(coverage // 2 + 1)
        )
    return probability


def has_bad_column(errors: list[int], called: list[int]) -> bool:
    """The antecedent of their Theorem 1: some column more than half wrong."""
    return any(c and e / c > 0.5 for e, c in zip(errors, called))


def truth_is_optimal(instance: SimulatedInstance) -> bool:
    matrix = instance.matrix
    return matrix.cost(instance.haplotype) <= dp.solve(matrix).cost + 1e-9


def run(n_variants, span, coverages, error_rates, replicates, seed_base):
    rows = []
    for coverage in coverages:
        for error_rate in error_rates:
            predicted = absent = optimal = gap_cases = 0.0
            for replicate in range(replicates):
                instance = simulate(
                    n_variants,
                    coverage=coverage,
                    read_span=span,
                    error_rate=error_rate,
                    seed=seed_base + replicate,
                    weight_model="uniform",
                )
                errors, called = column_errors(
                    instance.matrix, instance.haplotype, instance.origin
                )
                predicted += theorem2_probability(called, error_rate)
                no_bad_column = not has_bad_column(errors, called)
                truth_optimal = truth_is_optimal(instance)
                absent += no_bad_column
                optimal += truth_optimal
                gap_cases += no_bad_column and not truth_optimal
            rows.append(
                {
                    "coverage": coverage,
                    "error_rate": error_rate,
                    "thm2": predicted / replicates,
                    "no_bad_column": absent / replicates,
                    "truth_optimal": optimal / replicates,
                    "gap_cases": gap_cases / replicates,
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=int, default=30)
    parser.add_argument("--span", type=int, default=5)
    parser.add_argument("--coverages", type=float, nargs="+", default=[4.0, 6.0, 8.0])
    parser.add_argument(
        "--error-rates", type=float, nargs="+", default=[0.05, 0.10, 0.20]
    )
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--seed-base", type=int, default=9000)
    args = parser.parse_args()

    example = COUNTEREXAMPLE
    matrix = FragmentMatrix.from_rows(example["rows"])
    errors, called = column_errors(matrix, example["truth"], example["origin"])
    print(
        f"counterexample: column errors {errors} at coverages {called}; "
        f"MEC(truth) = {matrix.cost(example['truth']):g}, "
        f"MEC({''.join(map(str, example['better']))}) = "
        f"{matrix.cost(example['better']):g}"
    )
    print(
        f"{args.variants} variants, span {args.span}, coverages {args.coverages}, "
        f"{args.replicates} replicates, seed {args.seed_base} + replicate, "
        "uniform weights (unweighted MEC), exact DP"
    )
    header = (
        f"{'cov':>4} {'err':>5} {'thm2':>6} {'no-bad-col':>11} "
        f"{'truth-opt':>10} {'gap-cases':>10}"
    )
    print(header)
    print("-" * len(header))
    rows = run(
        args.variants,
        args.span,
        args.coverages,
        args.error_rates,
        args.replicates,
        args.seed_base,
    )
    for row in rows:
        print(
            f"{row['coverage']:>4g} {row['error_rate']:>5.2f} {row['thm2']:>6.3f} "
            f"{row['no_bad_column']:>11.3f} {row['truth_optimal']:>10.3f} "
            f"{row['gap_cases']:>10.3f}"
        )


if __name__ == "__main__":
    main()
