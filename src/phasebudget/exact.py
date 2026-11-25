"""Exhaustive wMEC solvers. These define correctness for everything faster."""

from __future__ import annotations

from itertools import product
from typing import Iterator

from ._types import WMECSolution
from .fragments import FragmentMatrix

#: Above this many variants, `solve_bruteforce` refuses rather than hangs.
MAX_BRUTEFORCE_VARIANTS = 22

#: Above this many fragments, `solve_bipartition_bruteforce` refuses.
MAX_BRUTEFORCE_FRAGMENTS = 20


def enumerate_haplotypes(n_variants: int) -> Iterator[tuple[int, ...]]:
    """Every haplotype with `h[0] = 0`, which breaks the complement symmetry.

    Under the complementary model `h` and its complement have equal cost, so
    fixing the first allele halves the search space and loses no optimum.
    """
    if n_variants == 0:
        yield ()
        return
    for tail in product((0, 1), repeat=n_variants - 1):
        yield (0,) + tail


def solve_bruteforce(
    matrix: FragmentMatrix,
    *,
    max_variants: int = MAX_BRUTEFORCE_VARIANTS,
) -> WMECSolution:
    """Optimal wMEC by enumerating haplotypes.

    The second haplotype is the complement of the first, so a solution is a
    haplotype plus the cheapest side for each read. Cost ties keep the first
    haplotype found in enumeration order, which makes the result reproducible
    but does not make it unique.
    """
    if matrix.n_variants > max_variants:
        raise ValueError(
            f"{matrix.n_variants} variants exceeds the brute-force limit of "
            f"{max_variants}; use phasebudget.dp.CoverageDP instead"
        )
    best_cost = float("inf")
    best_haplotype: tuple[int, ...] = ()
    for haplotype in enumerate_haplotypes(matrix.n_variants):
        cost = matrix.cost(haplotype)
        if cost < best_cost:
            best_cost, best_haplotype = cost, haplotype
    if best_cost == float("inf"):
        best_cost, best_haplotype = 0.0, tuple(0 for _ in range(matrix.n_variants))
    return WMECSolution(
        cost=best_cost,
        haplotype=best_haplotype,
        assignment=matrix.best_assignment(best_haplotype),
    )


def constrained_optimum(
    matrix: FragmentMatrix,
    boundary: int,
    relation: int,
    *,
    max_variants: int = MAX_BRUTEFORCE_VARIANTS,
) -> float:
    """Best wMEC cost subject to `h[boundary] ^ h[boundary+1] == relation`.

    This is the definition the block-wise margin extraction in
    `phasebudget.dp` is checked against. Returns infinity when no haplotype
    satisfies the constraint, which cannot happen for a valid boundary.
    """
    if not 0 <= boundary < matrix.n_variants - 1:
        raise ValueError(f"boundary {boundary} is outside the variant range")
    if relation not in (0, 1):
        raise ValueError(f"relation must be 0 or 1, got {relation}")
    if matrix.n_variants > max_variants:
        raise ValueError(f"{matrix.n_variants} variants exceeds the brute-force limit")
    best = float("inf")
    for haplotype in enumerate_haplotypes(matrix.n_variants):
        if haplotype[boundary] ^ haplotype[boundary + 1] != relation:
            continue
        best = min(best, matrix.cost(haplotype))
    return best


def bruteforce_margins(
    matrix: FragmentMatrix,
    *,
    max_variants: int = MAX_BRUTEFORCE_VARIANTS,
) -> tuple[float, ...]:
    """The MEC margin of every adjacent boundary, by constrained re-solving."""
    return tuple(
        abs(
            constrained_optimum(matrix, j, 0, max_variants=max_variants)
            - constrained_optimum(matrix, j, 1, max_variants=max_variants)
        )
        for j in range(matrix.n_variants - 1)
    )


def solve_bipartition_bruteforce(
    matrix: FragmentMatrix,
    *,
    max_fragments: int = MAX_BRUTEFORCE_FRAGMENTS,
) -> tuple[float, tuple[int, ...]]:
    """Optimal MEC under the general model, where the two haplotypes are free.

    Reads are split into two parts and each part is corrected to its own
    weighted column consensus, with no constraint tying the two consensus
    values at a variant. This is the classical MEC formulation. It lower bounds
    `solve_bruteforce`, and the gap is the price of assuming every variant is
    heterozygous. `tests/test_exact.py` carries an instance where the two
    models disagree.
    """
    if len(matrix) > max_fragments:
        raise ValueError(
            f"{len(matrix)} fragments exceeds the brute-force limit of {max_fragments}"
        )
    n = len(matrix)
    best_cost = float("inf")
    best_assignment: tuple[int, ...] = ()
    for code in range(1 << n) if n else (0,):
        assignment = tuple((code >> i) & 1 for i in range(n))
        cost = 0.0
        for position in range(matrix.n_variants):
            for side in (0, 1):
                to_zero = 0.0
                to_one = 0.0
                for index in matrix.active(position):
                    if assignment[index] != side:
                        continue
                    allele = matrix[index].allele_at(position)
                    if allele is None:
                        continue
                    weight = matrix[index].weight_at(position)
                    if allele == 0:
                        to_one += weight
                    else:
                        to_zero += weight
                cost += min(to_zero, to_one)
        if cost < best_cost:
            best_cost, best_assignment = cost, assignment
    if best_cost == float("inf"):
        return 0.0, ()
    return best_cost, best_assignment


def bruteforce_posterior_margins(
    matrix: FragmentMatrix,
    *,
    max_variants: int = MAX_BRUTEFORCE_VARIANTS,
) -> list[tuple[float, float]]:
    """`(log Z_0(b), log Z_1(b))` for every boundary, by enumerating haplotypes.

    `Z_r(b)` sums `exp(-cost(h, a))` over haplotypes `h` with `h[b] xor h[b+1]
    = r` and over every read assignment `a`. The sum over assignments factors
    read by read, so each haplotype contributes the product over reads of
    `exp(-cost(side 0)) + exp(-cost(side 1))`. Both a haplotype and its
    complement are enumerated, as the dynamic program counts both.
    `docs/design/posterior-margins.md` states the definition.
    """
    import math

    n = matrix.n_variants
    if n > max_variants:
        raise ValueError(f"{n} variants exceeds the brute-force limit of {max_variants}")
    if n < 2:
        return []
    terms: list[tuple[list[float], list[float]]] = [([], []) for _ in range(n - 1)]
    for code in range(1 << n):
        haplotype = tuple((code >> j) & 1 for j in range(n))
        log_weight = 0.0
        for index in range(len(matrix)):
            a = -matrix.fragment_cost(index, haplotype, 0)
            b = -matrix.fragment_cost(index, haplotype, 1)
            top = max(a, b)
            log_weight += top + math.log(math.exp(a - top) + math.exp(b - top))
        for boundary in range(n - 1):
            relation = haplotype[boundary] ^ haplotype[boundary + 1]
            terms[boundary][relation].append(log_weight)

    def logsumexp(values: list[float]) -> float:
        if not values:
            return float("-inf")
        top = max(values)
        return top + math.log(sum(math.exp(v - top) for v in values))

    return [(logsumexp(zero), logsumexp(one)) for zero, one in terms]


def solution_is_optimal(
    matrix: FragmentMatrix,
    solution: WMECSolution,
    *,
    tolerance: float = 1e-9,
    max_variants: int = MAX_BRUTEFORCE_VARIANTS,
) -> bool:
    """Check a solution against exhaustive enumeration, cost and consistency.

    A solution passes only when its recorded cost matches the cost of its own
    haplotype and assignment, and that cost is the enumerated optimum. A
    cheaper recorded cost fails the same way an expensive one does.
    """
    optimum = solve_bruteforce(matrix, max_variants=max_variants).cost
    stated = matrix.cost(solution.haplotype, solution.assignment)
    return (
        abs(stated - solution.cost) <= tolerance
        and abs(solution.cost - optimum) <= tolerance
    )
