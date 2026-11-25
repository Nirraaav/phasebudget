"""Connectivity-aware selection against exhaustive enumeration."""

import itertools
import random

import pytest

from phasebudget.selection import (
    is_feasible,
    read_weights,
    select_connected,
    select_optimal,
    selection_objective,
)

from conftest import random_matrix


def _objective(matrix, weights, selected, coverage, k, penalty):
    by_id = {f.read_id: f for f in matrix}
    positive = [f for f in matrix if weights[f.read_id] > 0]
    total = selection_objective(weights, selected)
    for b in range(matrix.n_variants - 1):
        if not any(f.start <= b and f.end - 1 >= b + 1 for f in positive):
            continue
        spanned = sum(1 for r in selected if by_id[r].start <= b and by_id[r].end - 1 >= b + 1)
        total -= penalty * max(0, min(k, coverage) - spanned)
    return total


def _best(matrix, weights, coverage, k, penalty):
    ids = [f.read_id for f in matrix]
    best = -float("inf")
    for size in range(len(ids) + 1):
        for subset in itertools.combinations(ids, size):
            if is_feasible(matrix, subset, coverage):
                best = max(best, _objective(matrix, weights, subset, coverage, k, penalty))
    return best


def test_the_flow_attains_the_exhaustive_maximum():
    rng = random.Random(21)
    for _ in range(250):
        matrix = random_matrix(rng, max_variants=8, max_reads=8)
        weights = read_weights(matrix)
        for coverage, k, penalty in ((1, 1, 0.5), (2, 1, 3.0), (2, 2, 1.0), (3, 1, 100.0)):
            result = select_connected(matrix, coverage, weights, k=k, penalty=penalty)
            assert is_feasible(matrix, result.selected, coverage)
            value = _objective(matrix, weights, result.selected, coverage, k, penalty)
            assert value == pytest.approx(_best(matrix, weights, coverage, k, penalty))


def test_zero_penalty_is_the_additive_optimum():
    rng = random.Random(22)
    for _ in range(200):
        matrix = random_matrix(rng, max_variants=9, max_reads=9)
        weights = read_weights(matrix)
        for coverage in (1, 2, 3):
            assert select_connected(matrix, coverage, weights, penalty=0.0).objective == \
                pytest.approx(select_optimal(matrix, coverage, weights).objective)


def test_a_hard_penalty_spans_every_boundary_the_cap_allows():
    from phasebudget import FragmentMatrix

    # Cap 1: the two heavy reads "0000--" and "--0000" overlap, so the additive
    # optimum keeps one of them plus a short read and leaves a break; the
    # connected selector keeps the long read that spans everything.
    matrix = FragmentMatrix.from_rows(
        ["0000--", "--0000", "000000", "----00"],
        weights=[[5, 5, 5, 5, 0, 0], [0, 0, 5, 5, 5, 5], [1, 1, 1, 1, 1, 1], [0, 0, 0, 0, 5, 5]],
    )
    weights = read_weights(matrix)
    additive = select_optimal(matrix, 1, weights)
    connected = select_connected(matrix, 1, weights, penalty=1000.0)
    assert connected.selected == (2,)
    assert additive.objective > connected.objective
