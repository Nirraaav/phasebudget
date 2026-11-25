import random

import pytest

from phasebudget import FragmentMatrix
from phasebudget.selection import (
    is_feasible,
    is_phase_informative,
    read_weights,
    select_bruteforce,
    select_greedy,
    select_optimal,
    select_random,
    selection_objective,
    span_score,
    variant_count_score,
)

from conftest import random_matrix


def test_optimal_selection_matches_exhaustive_enumeration():
    generator = random.Random(21)
    for _ in range(250):
        matrix = random_matrix(generator, max_variants=7, max_reads=7)
        weights = read_weights(matrix)
        for coverage in range(0, 4):
            optimal = select_optimal(matrix, coverage, weights)
            reference = select_bruteforce(matrix, coverage, weights)
            assert optimal.objective == pytest.approx(reference.objective)


def test_optimal_selection_is_feasible():
    generator = random.Random(22)
    for _ in range(250):
        matrix = random_matrix(generator)
        for coverage in range(0, 4):
            selected = select_optimal(matrix, coverage).selected
            assert is_feasible(matrix, selected, coverage)


def test_greedy_never_beats_the_exact_selector():
    generator = random.Random(23)
    for _ in range(250):
        matrix = random_matrix(generator)
        weights = read_weights(matrix)
        for coverage in range(1, 4):
            greedy = select_greedy(matrix, coverage, weights)
            optimal = select_optimal(matrix, coverage, weights)
            assert greedy.objective <= optimal.objective + 1e-9
            assert is_feasible(matrix, greedy.selected, coverage)


def test_greedy_is_strictly_suboptimal_on_a_blocking_read():
    # The long read scores highest and blocks both short reads, which together
    # are worth more. Every coverage-constrained greedy rule has an instance
    # of this shape; this one is the smallest.
    matrix = FragmentMatrix.from_rows(["0110", "01--", "--10"])
    weights = {0: 3.0, 1: 2.0, 2: 2.0}
    assert select_greedy(matrix, 1, weights).selected == (0,)
    assert select_optimal(matrix, 1, weights).selected == (1, 2)
    assert select_optimal(matrix, 1, weights).objective == 4.0


def test_a_cap_at_or_above_the_coverage_selects_every_useful_read():
    generator = random.Random(24)
    for _ in range(120):
        matrix = random_matrix(generator)
        coverage = matrix.max_coverage()
        selected = select_optimal(matrix, coverage).selected
        useful = tuple(
            fragment.read_id
            for fragment in matrix
            if is_phase_informative(fragment)
        )
        assert selected == useful


def test_a_read_calling_one_allele_is_never_selected():
    # It can always be assigned to the haplotype that matches it, so it costs
    # nothing under wMEC and constrains no phase relation, but it would still
    # consume coverage.
    matrix = FragmentMatrix.from_rows(["0110", "1---"])
    assert select_optimal(matrix, 4).selected == (0,)
    assert select_greedy(matrix, 4).selected == (0,)
    assert not is_phase_informative(matrix[1])


def test_zero_coverage_selects_nothing():
    matrix = FragmentMatrix.from_rows(["01", "10"])
    result = select_optimal(matrix, 0)
    assert result.selected == ()
    assert result.objective == 0.0


def test_non_positive_weights_are_never_selected():
    matrix = FragmentMatrix.from_rows(["01", "10"])
    weights = {0: 5.0, 1: 0.0}
    assert select_optimal(matrix, 2, weights).selected == (0,)


def test_selection_is_deterministic():
    matrix = FragmentMatrix.from_rows(["0110", "01--", "--10", "1001"])
    first = select_optimal(matrix, 1)
    second = select_optimal(matrix, 1)
    assert first == second


def test_random_selection_is_reproducible_and_feasible():
    matrix = FragmentMatrix.from_rows(["0110", "01--", "--10", "1001", "0-01"])
    first = select_random(matrix, 2, seed=17)
    assert first.selected == select_random(matrix, 2, seed=17).selected
    assert is_feasible(matrix, first.selected, 2)


def test_alternative_scores_change_the_selection():
    matrix = FragmentMatrix.from_rows(["0110", "0--1"])
    by_count = read_weights(matrix, variant_count_score)
    by_span = read_weights(matrix, span_score)
    assert by_count == {0: 4.0, 1: 2.0}
    assert by_span == {0: 4.0, 1: 4.0}


def test_objective_is_the_sum_of_the_selected_weights():
    weights = {0: 1.5, 1: 2.5, 2: 4.0}
    assert selection_objective(weights, (0, 2)) == 5.5


def test_a_missing_weight_is_an_error():
    matrix = FragmentMatrix.from_rows(["01", "10"])
    with pytest.raises(ValueError, match="no weight"):
        select_optimal(matrix, 1, {0: 1.0})


def test_negative_coverage_is_an_error():
    matrix = FragmentMatrix.from_rows(["01"])
    with pytest.raises(ValueError, match="non-negative"):
        select_optimal(matrix, -1)
