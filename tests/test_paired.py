"""The paired boundary-level comparison and its exact sign test."""

import math

import pytest

from phasebudget import Phasing
from phasebudget.paired import common_pairs, paired_comparison, sign_test, window_discordance


def test_sign_test_known_values():
    assert sign_test(0, 0) == 1.0
    assert sign_test(0, 10) == pytest.approx(2 / 1024)
    assert sign_test(10, 10) == pytest.approx(2 / 1024)
    assert sign_test(5, 10) == 1.0
    assert sign_test(2, 10) == pytest.approx(2 * (1 + 10 + 45) / 1024)


def test_sign_test_is_symmetric_and_bounded():
    for n in range(1, 30):
        for k in range(n + 1):
            p = sign_test(k, n)
            assert 0 < p <= 1
            assert p == pytest.approx(sign_test(n - k, n))


def test_shared_errors_drop_out():
    truth = (0, 1, 0, 1, 0)
    # Both phasings wrong at the second pair; only the second phasing wrong at the fourth.
    first = Phasing(alleles=(0, 1, 1, 0, 1), blocks=(0,) * 5)
    second = Phasing(alleles=(0, 1, 1, 0, 0), blocks=(0,) * 5)
    result = paired_comparison(truth, first, second)
    assert result.pairs == 4
    assert result.both_wrong == 1
    assert (result.only_first_wrong, result.only_second_wrong) == (0, 1)
    assert result.p_value == 1.0
    assert result.difference == pytest.approx(-0.25)


def test_only_pairs_inside_blocks_of_both_phasings_count():
    truth = (0, 1, 0, None, 1)
    first = Phasing(alleles=(0, 1, 0, 0, 1), blocks=(0, 0, 1, 1, 1))
    second = Phasing(alleles=(0, 1, 0, 1, 1), blocks=(0, 0, 0, 0, 0))
    rows = common_pairs(truth, first, second)
    # Pairs (0,1) and (2,4): the truth skips variant 3; (1,2) crosses a block of `first`.
    assert [q for q, _, _ in rows] == [1, 4]


def test_window_counts_sum_to_the_comparison():
    truth = tuple(i % 2 for i in range(20))
    first = Phasing(alleles=tuple((i % 2) ^ (i >= 10) for i in range(20)), blocks=(0,) * 20)
    second = Phasing(alleles=truth, blocks=(0,) * 20)
    rows = common_pairs(truth, first, second)
    windows = window_discordance(rows, tuple(range(0, 200, 10)), 50)
    assert sum(w[0] for w in windows) == len(rows)
    assert sum(w[1] for w in windows) == paired_comparison(truth, first, second).only_first_wrong
    assert sum(w[2] for w in windows) == 0
