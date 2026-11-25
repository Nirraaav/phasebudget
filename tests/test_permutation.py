"""The window sign-flip test in `benchmarks/permutation.py`."""

import importlib.util
import itertools
import random
import sys
from pathlib import Path

import pytest

from phasebudget.paired import sign_test


def _module():
    path = Path(__file__).parents[1] / "benchmarks" / "permutation.py"
    spec = importlib.util.spec_from_file_location("permutation", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["permutation"] = module  # dataclasses look the module up
    spec.loader.exec_module(module)
    return module


perm = _module()


def _enumerated(windows):
    """The oracle: every sign pattern of the non-zero windows, listed."""
    d = [a - b for _, a, b in windows if a != b]
    target = abs(sum(d))
    hits = sum(abs(sum(s * x for s, x in zip(signs, d))) >= target
               for signs in itertools.product((1, -1), repeat=len(d)))
    return hits / 2 ** len(d)


def test_fixed_small_case_by_hand():
    # d = (3, 1, 0): sums over signs are 4, 2, -2, -4; |.| >= 4 in 2 of 4.
    result = perm.window_sign_flip_test([(10, 3, 0), (5, 2, 1), (7, 1, 1)])
    assert result.exact and result.windows == 2 and result.difference == 4
    assert result.p_value == 0.5


def test_no_discordance_gives_one():
    result = perm.window_sign_flip_test([(4, 0, 0), (3, 2, 2)])
    assert result.p_value == 1.0 and result.windows == 0 and result.exact


def test_one_heavy_window_is_not_significant():
    # 30 discordant pairs in one window: the sign test says p ~ 2e-9, the
    # window test cannot go below one half.
    windows = [(100, 30, 0), (100, 0, 0)]
    assert sign_test(30, 30) < 1e-8
    assert perm.window_sign_flip_test(windows).p_value == 1.0


def test_exact_matches_enumeration():
    rng = random.Random(607)
    for _ in range(200):
        windows = [(50, rng.randint(0, 4), rng.randint(0, 4))
                   for _ in range(rng.randint(0, 12))]
        result = perm.window_sign_flip_test(windows)
        assert result.exact
        assert result.p_value == pytest.approx(_enumerated(windows), abs=1e-12)


def test_agrees_with_the_sign_test_at_one_pair_per_window():
    rng = random.Random(1)
    for _ in range(200):
        windows = [(9, *rng.choice([(1, 0), (0, 1), (0, 0)])) for _ in range(rng.randint(0, 20))]
        first = sum(w[1] for w in windows)
        second = sum(w[2] for w in windows)
        result = perm.window_sign_flip_test(windows)
        assert result.exact
        assert result.p_value == pytest.approx(sign_test(first, first + second), abs=1e-12)


def test_monte_carlo_near_the_sign_test_and_deterministic():
    windows = [(9, 1, 0)] * 25 + [(9, 0, 1)] * 15
    a = perm.window_sign_flip_test(windows, resamples=20000, seed=3)
    b = perm.window_sign_flip_test(windows, resamples=20000, seed=3)
    assert not a.exact and a == b
    assert a.p_value == pytest.approx(sign_test(25, 40), abs=0.01)
    assert perm.window_sign_flip_test(windows, resamples=500, seed=4) != a


def test_rejects_no_resamples():
    with pytest.raises(ValueError):
        perm.window_sign_flip_test([(1, 1, 0)], resamples=0)
