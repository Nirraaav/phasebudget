"""End-to-end checks that the parts compose the way the proposal claims."""

import pytest

from phasebudget import CoverageDP
from phasebudget.adaptive import adaptive_selection, static_selection
from phasebudget.metrics import phasing_from_solution, switch_positions, summarize
from phasebudget.selection import (
    is_feasible,
    read_weights,
    select_fischer_marschall,
    select_greedy,
    select_max_min_coverage,
    select_optimal,
    select_random,
)
from phasebudget.simulate import simulate_connected


def instance():
    return simulate_connected(
        24, coverage=12.0, read_span=6, error_rate=0.05, seed=101
    )


def test_downsampling_to_the_cap_makes_the_dp_tractable():
    data = instance()
    assert data.matrix.max_coverage() > 6
    subset = data.matrix.subset(select_optimal(data.matrix, 4).selected)
    assert subset.max_coverage() <= 4
    CoverageDP(subset, max_active_reads=4).solve()


def test_exact_selection_dominates_the_heuristics_on_the_surrogate_objective():
    data = instance()
    weights = read_weights(data.matrix)
    exact = select_optimal(data.matrix, 4, weights).objective
    greedy = select_greedy(data.matrix, 4, weights).objective
    sampled = select_random(data.matrix, 4, seed=1, weights=weights).objective
    assert exact >= greedy >= 0.0
    assert exact >= sampled


def test_error_free_data_phases_without_a_switch_error():
    data = simulate_connected(20, coverage=8.0, read_span=6, error_rate=1e-9, seed=102)
    subset = data.matrix.subset(select_optimal(data.matrix, 4).selected)
    phasing = phasing_from_solution(subset, CoverageDP(subset).solve())
    assert switch_positions(data.haplotype, phasing) == ()
    assert summarize(data.haplotype, phasing)["phased_fraction"] == 1.0


@pytest.mark.parametrize("selector", [select_fischer_marschall, select_max_min_coverage])
def test_published_selectors_feed_a_feasible_subset_to_the_dp(selector):
    data = instance()
    selection = selector(data.matrix, 4, read_weights(data.matrix))
    subset = data.matrix.subset(selection.selected)
    assert is_feasible(data.matrix, selection.selected, 4)
    assert CoverageDP(subset).solve().cost >= 0.0


def test_the_adaptive_loop_runs_the_whole_pipeline():
    data = instance()
    result = adaptive_selection(data.matrix, 4, iterations=3)
    phasing = phasing_from_solution(
        data.matrix.subset(result.selected), result.solution
    )
    report = summarize(data.haplotype, phasing, data.positions)
    assert report["phased_fraction"] > 0.0
    assert report["block_n50"] > 0.0
    assert result.full_cost <= static_selection(data.matrix, 4).full_cost + 1e-9


def test_a_higher_cap_never_raises_the_full_matrix_cost_of_the_solution():
    data = instance()
    costs = []
    for cap in (2, 3, 4, 5):
        result = static_selection(data.matrix, cap)
        costs.append(result.full_cost)
    assert costs[-1] <= costs[0] + 1e-9


@pytest.mark.parametrize("cap", [2, 3, 4])
def test_every_cap_produces_a_feasible_and_scored_result(cap):
    data = instance()
    result = static_selection(data.matrix, cap)
    subset = data.matrix.subset(result.selected)
    assert subset.max_coverage() <= cap
    assert result.full_cost >= 0.0
