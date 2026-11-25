"""The general-model dynamic program answers to the bipartition oracle."""

import random

import pytest

from phasebudget.dp import CoverageDP
from phasebudget.exact import solve_bipartition_bruteforce
from phasebudget.fragments import FragmentMatrix
from phasebudget.general import GeneralDP, assignment_cost, solve_general
from phasebudget._validation import CoverageLimitExceeded
from phasebudget.simulate import simulate, simulate_with_genotype_errors


def _small_instances():
    rng = random.Random(4242)
    for case in range(80):
        if case % 2:
            yield simulate_with_genotype_errors(
                rng.randint(2, 10),
                coverage=rng.choice([2.0, 3.0, 4.0]),
                read_span=rng.randint(2, 5),
                error_rate=rng.choice([0.01, 0.1, 0.2]),
                homozygous_rate=rng.choice([0.0, 0.2, 0.4]),
                seed=case,
            ).matrix
        else:
            yield simulate(
                rng.randint(2, 10),
                coverage=rng.choice([2.0, 3.0, 4.0]),
                read_span=rng.randint(2, 5),
                error_rate=rng.choice([0.01, 0.1, 0.3]),
                gap_rate=rng.choice([0.0, 0.2]),
                seed=case,
            ).matrix


def test_costs_equal_the_bipartition_oracle():
    checked = 0
    for matrix in _small_instances():
        if len(matrix) > 16:
            continue
        expected, _ = solve_bipartition_bruteforce(matrix)
        cost, _ = solve_general(matrix)
        assert cost == pytest.approx(expected, abs=1e-9)
        checked += 1
    assert checked >= 40


def test_the_returned_bipartition_achieves_the_returned_cost():
    for matrix in _small_instances():
        cost, assignment = solve_general(matrix)
        assert len(assignment) == len(matrix)
        assert assignment_cost(matrix, assignment) == pytest.approx(cost, abs=1e-9)


def test_the_general_model_never_costs_more_than_the_complementary_one():
    for matrix in _small_instances():
        assert solve_general(matrix)[0] <= CoverageDP(matrix).solve().cost + 1e-9


def test_a_homozygous_looking_column_costs_nothing():
    matrix = FragmentMatrix.from_rows(["011", "111", "110"])
    cost, _ = solve_general(matrix)
    assert cost == pytest.approx(solve_bipartition_bruteforce(matrix)[0])
    assert cost < CoverageDP(matrix).solve().cost


def test_empty_matrix_and_the_coverage_limit():
    assert solve_general(FragmentMatrix([])) == (0.0, ())
    deep = FragmentMatrix.from_rows(["01"] * 9)
    with pytest.raises(CoverageLimitExceeded):
        GeneralDP(deep, max_active_reads=8)


def test_it_scales_past_the_oracle():
    matrix = simulate_with_genotype_errors(
        60, coverage=8.0, read_span=5, error_rate=0.05, homozygous_rate=0.1, seed=3
    ).matrix
    assert len(matrix) > 60
    cost, assignment = solve_general(matrix)
    assert assignment_cost(matrix, assignment) == pytest.approx(cost)
