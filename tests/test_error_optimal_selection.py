"""Selection under convex boundary costs against exhaustive enumeration."""

import itertools
import random

import pytest

from phasebudget.selection import (
    is_feasible,
    read_weights,
    select_connected,
    select_error_optimal,
    select_optimal,
    selection_objective,
    switch_model_costs,
)

from conftest import random_matrix


def _objective(matrix, weights, selected, costs):
    positive = [f for f in matrix if weights[f.read_id] > 0]
    by_id = {f.read_id: f for f in matrix}
    total = selection_objective(weights, selected)
    for b in range(matrix.n_variants - 1):
        if not any(f.start <= b and f.end - 1 >= b + 1 for f in positive):
            continue
        spanned = sum(1 for r in selected if by_id[r].start <= b and by_id[r].end - 1 >= b + 1)
        total -= costs[spanned]
    return total


def _best(matrix, weights, coverage, costs):
    ids = [f.read_id for f in matrix]
    return max(_objective(matrix, weights, subset, costs)
               for size in range(len(ids) + 1)
               for subset in itertools.combinations(ids, size)
               if is_feasible(matrix, subset, coverage))


def _random_convex_costs(rng, coverage):
    """Non-increasing convex: the decrements shrink from one count to the next."""
    steps = sorted((rng.choice([0.0, rng.uniform(0, 6)]) for _ in range(coverage)), reverse=True)
    costs = [sum(steps)]
    for step in steps:
        costs.append(costs[-1] - step)
    return costs


def test_the_flow_attains_the_exhaustive_maximum():
    rng = random.Random(31)
    for _ in range(250):
        matrix = random_matrix(rng, max_variants=8, max_reads=8)
        weights = read_weights(matrix)
        for coverage in (1, 2, 3):
            for costs in (_random_convex_costs(rng, coverage),
                          switch_model_costs(coverage, 0.05, rng.uniform(0, 50))):
                result = select_error_optimal(matrix, coverage, weights, costs=costs)
                assert is_feasible(matrix, result.selected, coverage)
                assert _objective(matrix, weights, result.selected, costs) == \
                    pytest.approx(_best(matrix, weights, coverage, costs))


def test_the_step_cost_is_the_connected_selector():
    rng = random.Random(32)
    for _ in range(150):
        matrix = random_matrix(rng, max_variants=9, max_reads=9)
        weights = read_weights(matrix)
        for coverage, k, penalty in ((1, 1, 2.0), (2, 2, 0.7), (3, 1, 40.0)):
            costs = [penalty * max(k - s, 0) for s in range(coverage + 1)]
            ours = select_error_optimal(matrix, coverage, weights, costs=costs)
            theirs = select_connected(matrix, coverage, weights, k=k, penalty=penalty)
            assert _objective(matrix, weights, ours.selected, costs) == \
                pytest.approx(_objective(matrix, weights, theirs.selected, costs))


def test_flat_costs_give_the_additive_optimum():
    rng = random.Random(33)
    for _ in range(100):
        matrix = random_matrix(rng, max_variants=9, max_reads=9)
        weights = read_weights(matrix)
        for coverage in (1, 2, 4):
            costs = [3.0] * (coverage + 1)
            assert select_error_optimal(matrix, coverage, weights, costs=costs).objective == \
                pytest.approx(select_optimal(matrix, coverage, weights).objective)


def test_a_heavy_cost_trades_weight_for_a_second_spanning_read():
    from phasebudget import FragmentMatrix

    # Cap 2. The heavy read "00-" spans only boundary 0, so the additive
    # optimum keeps it with the better long read and leaves boundary 1 with
    # one spanning read; charging that keeps both long reads instead.
    matrix = FragmentMatrix.from_rows(
        ["000", "000", "00-"],
        weights=[[1, 1, 1], [1, 1, 2], [10, 10, 0]],
    )
    weights = read_weights(matrix)
    additive = select_optimal(matrix, 2, weights)
    model = select_error_optimal(matrix, 2, weights, costs=[1000.0, 100.0, 0.0])
    assert set(additive.selected) == {1, 2}
    assert set(model.selected) == {0, 1}


def test_costs_that_break_exactness_are_rejected():
    from phasebudget import FragmentMatrix

    matrix = FragmentMatrix.from_rows(["00", "11"])
    with pytest.raises(ValueError, match="convex"):
        select_error_optimal(matrix, 2, costs=[2.0, 1.9, 0.0])
    with pytest.raises(ValueError, match="non-increasing"):
        select_error_optimal(matrix, 2, costs=[0.0, 1.0, 2.0])
    with pytest.raises(ValueError, match="one entry"):
        select_error_optimal(matrix, 2, costs=[1.0, 0.0])


def test_the_model_costs_are_the_bhattacharyya_bound():
    costs = switch_model_costs(3, 0.1, 2.0)
    q = 2 * 0.1 * 0.9
    rho = 2 * (q * (1 - q)) ** 0.5
    assert costs == pytest.approx([2.0, 2.0 * rho, 2.0 * rho ** 2, 2.0 * rho ** 3])
