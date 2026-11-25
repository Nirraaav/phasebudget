import random

import pytest

from phasebudget import CoverageDP, FragmentMatrix
from phasebudget._validation import CoverageLimitExceeded
from phasebudget.exact import bruteforce_margins, solve_bruteforce
from phasebudget.simulate import simulate

from conftest import random_matrix


def test_dp_matches_brute_force_on_random_instances():
    generator = random.Random(3)
    for _ in range(300):
        matrix = random_matrix(generator)
        expected = solve_bruteforce(matrix)
        got = CoverageDP(matrix).solve()
        assert got.cost == pytest.approx(expected.cost)


def test_dp_assignment_reproduces_its_own_cost():
    generator = random.Random(5)
    for _ in range(200):
        matrix = random_matrix(generator)
        solution = CoverageDP(matrix).solve()
        assert matrix.cost(solution.haplotype, solution.assignment) == pytest.approx(
            solution.cost
        )


def test_dp_haplotype_is_an_optimal_haplotype():
    generator = random.Random(6)
    for _ in range(200):
        matrix = random_matrix(generator)
        solution = CoverageDP(matrix).solve()
        assert matrix.cost(solution.haplotype) == pytest.approx(solution.cost)


def test_block_wise_margins_match_constrained_re_solving():
    generator = random.Random(8)
    for _ in range(200):
        matrix = random_matrix(generator)
        got = tuple(margin.margin for margin in CoverageDP(matrix).margins())
        assert got == pytest.approx(bruteforce_margins(matrix))


def test_margin_relation_agrees_with_an_optimal_haplotype():
    generator = random.Random(9)
    for _ in range(120):
        matrix = random_matrix(generator)
        dp = CoverageDP(matrix)
        solution = dp.solve()
        for margin in dp.margins():
            if margin.margin > 1e-9:
                relation = (
                    solution.haplotype[margin.boundary]
                    ^ solution.haplotype[margin.boundary + 1]
                )
                assert relation == margin.relation


def test_boundary_costs_bracket_the_optimum():
    generator = random.Random(10)
    for _ in range(80):
        matrix = random_matrix(generator)
        dp = CoverageDP(matrix)
        optimum = dp.solve().cost
        for boundary in range(matrix.n_variants - 1):
            zero, one = dp.boundary_costs(boundary)
            assert min(zero, one) == pytest.approx(optimum)


def test_unspanned_boundaries_have_zero_margin():
    matrix = FragmentMatrix.from_rows(["01--", "10--", "--01", "--10"])
    margins = CoverageDP(matrix).margins()
    assert margins[1].margin == 0.0
    assert margins[1].spanned is False


def test_gapped_reads_still_link_the_boundary_they_span():
    matrix = FragmentMatrix.from_rows(["0-1", "0-1", "011", "100"])
    margins = CoverageDP(matrix).margins()
    assert all(margin.spanned for margin in margins)
    assert all(margin.margin > 0.0 for margin in margins)


def test_a_spanned_boundary_can_still_have_zero_margin():
    # Variant 1 is called by one read on each side, so its allele is free and
    # both phase relations cost the same. Zero margin means the decision is
    # unsupported, which is not the same as the boundary being unspanned.
    matrix = FragmentMatrix.from_rows(["0-1", "0-1", "010", "101"])
    margins = CoverageDP(matrix).margins()
    assert all(margin.spanned for margin in margins)
    assert all(margin.margin == 0.0 for margin in margins)


def test_the_dp_refuses_more_coverage_than_it_can_hold():
    matrix = simulate(4, coverage=40.0, read_span=4, error_rate=0.01, seed=1).matrix
    assert matrix.max_coverage() > 8
    with pytest.raises(CoverageLimitExceeded):
        CoverageDP(matrix, max_active_reads=8)


def test_empty_matrix_solves_to_zero():
    solution = CoverageDP(FragmentMatrix([])).solve()
    assert solution.cost == 0.0
    assert solution.haplotype == ()


def test_solving_an_error_free_simulation_recovers_the_haplotype():
    instance = simulate(10, coverage=6.0, read_span=5, error_rate=1e-9, seed=4)
    solution = CoverageDP(instance.matrix).solve()
    assert solution.cost == pytest.approx(0.0, abs=1e-6)
    truth = instance.haplotype
    flipped = instance.complement()
    assert solution.haplotype in (truth, flipped)
