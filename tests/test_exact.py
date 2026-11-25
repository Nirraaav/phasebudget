import random

import pytest

from phasebudget import FragmentMatrix
from phasebudget.exact import (
    bruteforce_margins,
    constrained_optimum,
    enumerate_haplotypes,
    solve_bipartition_bruteforce,
    solve_bruteforce,
    solution_is_optimal,
)

from conftest import random_matrix


def test_enumeration_fixes_the_first_allele():
    haplotypes = list(enumerate_haplotypes(3))
    assert len(haplotypes) == 4
    assert all(haplotype[0] == 0 for haplotype in haplotypes)


def test_error_free_instance_has_zero_cost():
    matrix = FragmentMatrix.from_rows(["0110", "0110", "1001"])
    solution = solve_bruteforce(matrix)
    assert solution.cost == 0.0
    assert solution.haplotype == (0, 1, 1, 0)


def test_one_flipped_allele_costs_its_weight():
    matrix = FragmentMatrix.from_rows(["0110", "0110", "0100"], [[2.0] * 4] * 3)
    assert solve_bruteforce(matrix).cost == 2.0


def test_a_single_read_is_always_free():
    matrix = FragmentMatrix.from_rows(["0101"])
    assert solve_bruteforce(matrix).cost == 0.0


def test_reported_cost_matches_the_reported_assignment():
    matrix = FragmentMatrix.from_rows(["0110", "0100", "1011"])
    solution = solve_bruteforce(matrix)
    assert matrix.cost(solution.haplotype, solution.assignment) == solution.cost
    assert solution_is_optimal(matrix, solution)


def test_solution_is_optimal_rejects_an_understated_cost():
    matrix = FragmentMatrix.from_rows(["0110", "0100", "1011"])
    solution = solve_bruteforce(matrix)
    understated = type(solution)(
        cost=solution.cost - 1.0,
        haplotype=solution.haplotype,
        assignment=solution.assignment,
    )
    assert not solution_is_optimal(matrix, understated)


def test_constrained_optimum_bounds_the_unconstrained_one():
    matrix = FragmentMatrix.from_rows(["011", "010", "101"])
    optimum = solve_bruteforce(matrix).cost
    for boundary in range(matrix.n_variants - 1):
        costs = [constrained_optimum(matrix, boundary, side) for side in (0, 1)]
        assert min(costs) == pytest.approx(optimum)
        assert max(costs) >= optimum


def test_margins_are_zero_at_an_unspanned_boundary():
    matrix = FragmentMatrix.from_rows(["01--", "--10"])
    assert bruteforce_margins(matrix)[1] == 0.0


def test_the_two_models_differ_on_a_homozygous_looking_column():
    # Every read calls allele 1 at the middle variant, so the general model
    # gives both haplotypes a 1 there and pays nothing. The complementary
    # model cannot, and pays the weight of the minority side.
    matrix = FragmentMatrix.from_rows(["011", "111", "110"])
    complementary = solve_bruteforce(matrix).cost
    general, _ = solve_bipartition_bruteforce(matrix)
    assert general < complementary


def test_the_general_model_never_costs_more():
    generator = random.Random(11)
    for _ in range(60):
        matrix = random_matrix(generator, max_variants=6, max_reads=5)
        general, _ = solve_bipartition_bruteforce(matrix)
        assert general <= solve_bruteforce(matrix).cost + 1e-9


def test_brute_force_refuses_an_instance_it_cannot_finish():
    matrix = FragmentMatrix.from_rows(["0" * 30])
    with pytest.raises(ValueError, match="brute-force limit"):
        solve_bruteforce(matrix)
