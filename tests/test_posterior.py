"""Posterior margins answer to the enumeration in `exact`."""

import math
import random

import pytest

from phasebudget.dp import CoverageDP
from phasebudget.exact import bruteforce_posterior_margins
from phasebudget.fragments import FragmentMatrix
from phasebudget.posterior import boundary_log_sums, logsumexp, posterior_margins
from phasebudget.simulate import simulate


def _instances():
    rng = random.Random(99)
    for case in range(60):
        yield simulate(
            rng.randint(2, 11),
            coverage=rng.choice([2.0, 3.0, 5.0]),
            read_span=rng.randint(2, 5),
            error_rate=rng.choice([0.01, 0.1, 0.25]),
            gap_rate=rng.choice([0.0, 0.2]),
            quality_spread=rng.choice([0.0, 2.0]),
            seed=900 + case,
        ).matrix


def test_log_sums_equal_the_enumeration():
    for matrix in _instances():
        got = boundary_log_sums(matrix)
        expected = bruteforce_posterior_margins(matrix)
        assert len(got) == len(expected)
        for (zero, one), (want_zero, want_one) in zip(got, expected):
            assert zero == pytest.approx(want_zero, rel=1e-9, abs=1e-9)
            assert one == pytest.approx(want_one, rel=1e-9, abs=1e-9)


def test_both_relations_sum_to_the_same_total_at_every_boundary():
    for matrix in _instances():
        totals = [logsumexp(pair) for pair in boundary_log_sums(matrix)]
        for total in totals:
            assert total == pytest.approx(totals[0], rel=1e-9, abs=1e-9)


def test_an_unspanned_boundary_has_zero_posterior_margin():
    matrix = FragmentMatrix.from_rows(["01--", "10--", "--01", "--11"])
    margins = posterior_margins(matrix)
    assert margins[1].margin == 0.0
    assert not margins[1].spanned


def test_the_posterior_margin_never_exceeds_the_mec_margin_by_more_than_log_count():
    # Z_r is at least the best term and at most that term times the number of
    # solutions, so the two margins differ by at most the log of that count.
    for matrix in list(_instances())[:30]:
        hard = CoverageDP(matrix).margins()
        soft = posterior_margins(matrix)
        bound = matrix.n_variants * math.log(2) + len(matrix) * math.log(2)
        for h, s in zip(hard, soft):
            assert abs(h.margin - s.margin) <= bound + 1e-9


def test_with_a_dominant_optimum_the_posterior_margin_is_bounded_by_the_mec_margin():
    # Error-free reads with large weights: the optimum carries almost all the
    # mass of its relation, so log Z_map is its cost. The reversed relation
    # has several equally cheap alternatives, since a spanning read can take
    # its forced mismatch on either side, and the posterior counts them all.
    # So the posterior margin is at most the MEC margin, and below it by at
    # most the log of the number of alternatives.
    rows = ["0101", "1010", "0101", "1010", "0101"]
    weights = [[30.0] * 4 for _ in rows]
    matrix = FragmentMatrix.from_rows(rows, weights=weights)
    hard = CoverageDP(matrix).margins()
    soft = posterior_margins(matrix)
    for h, s in zip(hard, soft):
        assert s.relation == h.relation
        assert s.margin <= h.margin + 1e-6
        assert s.margin >= h.margin - len(rows) * math.log(2) - 1e-6
