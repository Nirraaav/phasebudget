"""The two published read selectors, against their own definitions.

Fischer and Marschall 2016 and Mäkinen, Staneva, Tomescu and Valenzuela 2015
optimize different objectives from each other and from `select_optimal`.
Checking each against its own definition is the only comparison that means
anything; `benchmarks/benchmark_selection.py` reports all three side by side.
"""

import random

import pytest

from phasebudget import FragmentMatrix
from phasebudget.selection import (
    component_count,
    distinct_variants_covered,
    fischer_marschall_scores,
    is_feasible,
    min_coverage,
    read_weights,
    select_fischer_marschall,
    select_max_min_coverage,
    select_max_min_coverage_bruteforce,
    select_optimal,
    select_random,
)
from phasebudget.simulate import simulate

from conftest import random_matrix


# --- Fischer and Marschall 2016 -------------------------------------------


def test_the_static_score_subtracts_holes():
    matrix = FragmentMatrix.from_rows(["0-1-0", "01---"])
    dynamic, static, _ = fischer_marschall_scores(matrix[0], set())
    assert static == 3 - 2
    assert dynamic == static
    assert fischer_marschall_scores(matrix[1], set())[1] == 2


def test_the_dynamic_score_discounts_variants_already_covered():
    matrix = FragmentMatrix.from_rows(["0110"])
    assert fischer_marschall_scores(matrix[0], set())[0] == 4
    assert fischer_marschall_scores(matrix[0], {0, 1})[0] == 2
    assert fischer_marschall_scores(matrix[0], {0, 1, 2, 3})[0] == 0


def test_the_quality_score_is_the_worst_weight():
    matrix = FragmentMatrix.from_rows(["011"], [[5.0, 2.0, 9.0]])
    assert fischer_marschall_scores(matrix[0], set())[2] == 2.0


def test_selection_respects_the_cap_at_every_variant():
    generator = random.Random(51)
    for _ in range(120):
        matrix = random_matrix(generator, max_variants=8, max_reads=8)
        for coverage in (1, 2, 3):
            selected = select_fischer_marschall(matrix, coverage).selected
            assert is_feasible(matrix, selected, coverage)


def test_bridging_is_not_monotone_in_the_selected_set():
    # Bridging reads spend coverage in the slice that adds them, so a later
    # slice can no longer afford reads it would otherwise have taken. The
    # selection with bridging is therefore not a superset of the one without,
    # and any comparison of the two has to be made on outcomes rather than on
    # containment.
    generator = random.Random(52)
    differed = 0
    for _ in range(80):
        matrix = random_matrix(generator, max_variants=8, max_reads=8)
        with_bridge = set(select_fischer_marschall(matrix, 2).selected)
        without = set(select_fischer_marschall(matrix, 2, bridging=False).selected)
        assert is_feasible(matrix, tuple(sorted(with_bridge)), 2)
        differed += not without <= with_bridge
    assert differed > 0


def test_bridging_joins_two_blocks_a_slice_left_apart():
    # The first slice covers every variant with two disjoint reads and stops.
    # Only the spanning read joins the two components, and it scores worse on
    # every component of the tuple, so nothing but bridging picks it up.
    matrix = FragmentMatrix.from_rows(["01--", "--10", "-11-"])
    without = select_fischer_marschall(matrix, 1, bridging=False)
    assert 2 not in without.selected
    assert component_count(matrix, without.selected) == 2
    with_bridge = select_fischer_marschall(matrix, 2)
    assert 2 in with_bridge.selected
    assert component_count(matrix, with_bridge.selected) == 1


def test_phase_uninformative_reads_are_dropped_first():
    matrix = FragmentMatrix.from_rows(["0110", "1---", "--0-"])
    assert select_fischer_marschall(matrix, 3).selected == (0,)


def test_the_selector_covers_more_distinct_variants_than_random():
    # Goal 1 of the paper. Averaged over instances, not guaranteed per instance.
    generator = random.Random(53)
    published, sampled = 0, 0
    for trial in range(60):
        matrix = random_matrix(generator, max_variants=10, max_reads=10)
        published += distinct_variants_covered(
            matrix, select_fischer_marschall(matrix, 2).selected
        )
        sampled += distinct_variants_covered(
            matrix, select_random(matrix, 2, seed=trial).selected
        )
    assert published >= sampled


def test_selection_is_deterministic():
    instance = simulate(20, coverage=8.0, read_span=5, seed=54)
    first = select_fischer_marschall(instance.matrix, 3)
    assert first.selected == select_fischer_marschall(instance.matrix, 3).selected


def test_the_published_selector_does_not_maximize_the_additive_objective():
    # This is OQ-2 as a test. The two selectors solve different problems, so
    # neither dominating the other on one objective is the expected outcome.
    instance = simulate(30, coverage=12.0, read_span=6, error_rate=0.05, seed=55)
    weights = read_weights(instance.matrix)
    published = select_fischer_marschall(instance.matrix, 3, weights)
    exact = select_optimal(instance.matrix, 3, weights)
    assert published.objective < exact.objective


# --- Mäkinen, Staneva, Tomescu and Valenzuela 2015 ------------------------


def test_max_min_coverage_matches_enumeration():
    generator = random.Random(56)
    for _ in range(150):
        matrix = random_matrix(generator, max_variants=6, max_reads=6)
        for coverage in range(0, 4):
            result = select_max_min_coverage(matrix, coverage)
            expected, _ = select_max_min_coverage_bruteforce(matrix, coverage)
            assert min_coverage(matrix, result.selected) == expected
            assert is_feasible(matrix, result.selected, coverage)


def test_max_min_coverage_matches_enumeration_on_the_unrestricted_problem():
    # Without the phase-informative filter this is Problem 1 of the paper
    # exactly, including reads that call a single allele.
    generator = random.Random(57)
    for _ in range(150):
        matrix = random_matrix(generator, max_variants=6, max_reads=6)
        for coverage in range(0, 4):
            result = select_max_min_coverage(
                matrix, coverage, informative_only=False
            )
            expected, _ = select_max_min_coverage_bruteforce(
                matrix, coverage, informative_only=False
            )
            assert min_coverage(matrix, result.selected) == expected


def test_the_floor_never_exceeds_the_cap():
    instance = simulate(20, coverage=10.0, read_span=6, seed=58)
    for coverage in (1, 2, 4):
        result = select_max_min_coverage(instance.matrix, coverage)
        assert min_coverage(instance.matrix, result.selected) <= coverage


def test_an_uncovered_variant_forces_the_floor_to_zero():
    matrix = FragmentMatrix.from_rows(["01---", "---01"])
    assert matrix.coverage()[2] == 0
    result = select_max_min_coverage(matrix, 2)
    assert min_coverage(matrix, result.selected) == 0


def test_weight_breaks_ties_without_changing_the_floor():
    instance = simulate(16, coverage=8.0, read_span=5, error_rate=0.03, seed=59)
    weights = read_weights(instance.matrix)
    tied = select_max_min_coverage(instance.matrix, 3, weights)
    plain = select_max_min_coverage(
        instance.matrix, 3, weights, break_ties_by_weight=False
    )
    assert min_coverage(instance.matrix, tied.selected) == min_coverage(
        instance.matrix, plain.selected
    )
    assert tied.objective >= plain.objective


def test_the_two_exact_selectors_disagree():
    # Maximum weight and maximum minimum coverage are different objectives, so
    # the exact answer to one is not the exact answer to the other.
    instance = simulate(30, coverage=12.0, read_span=6, error_rate=0.05, seed=60)
    weights = read_weights(instance.matrix)
    by_weight = select_optimal(instance.matrix, 3, weights)
    by_floor = select_max_min_coverage(instance.matrix, 3, weights)
    assert by_weight.selected != by_floor.selected
    assert by_weight.objective >= by_floor.objective
    assert min_coverage(instance.matrix, by_floor.selected) >= min_coverage(
        instance.matrix, by_weight.selected
    )


def test_zero_coverage_selects_nothing():
    matrix = FragmentMatrix.from_rows(["011", "110"])
    assert select_max_min_coverage(matrix, 0).selected == ()


def test_negative_coverage_is_an_error():
    matrix = FragmentMatrix.from_rows(["011"])
    with pytest.raises(ValueError, match="non-negative"):
        select_max_min_coverage(matrix, -1)
