"""L4 planner (`docs/design/l4-calibration-planner.md`) against brute force over k."""

import math
import random
from fractions import Fraction

from phasebudget.confidence import (Units, certification_probability, critical_count,
                                    learn_then_test, units_needed)


def exact_cdf(k, n, p):
    p = Fraction(p)
    return sum(math.comb(n, i) * p**i * (1 - p)**(n - i) for i in range(k + 1))


def oracle(n, r, beta, delta):
    """Brute force: enumerate every wrong-unit count and run the real test on it."""
    total = Fraction(0)
    for k in range(n + 1):
        units = Units(tuple([1] * n), tuple([1] * k))  # one threshold, k wrong units
        if learn_then_test(units, [1], beta, delta) is not None:
            total += math.comb(n, k) * Fraction(r)**k * (1 - Fraction(r))**(n - k)
    return float(total)


def test_fixed_examples():
    # r = 0 floor: (1 - beta)^n <= delta first holds at n = 2995 for beta 1e-3, delta 0.05.
    assert critical_count(2994, 1e-3, 0.05) == -1
    assert critical_count(2995, 1e-3, 0.05) == 0
    assert certification_probability(2995, 0.0, 1e-3, 0.05) == 1.0
    assert certification_probability(100, 0.0, 1e-3, 0.05) == 0.0
    assert units_needed(0.0, 1e-3, 0.05, 0.95) == 2995
    assert units_needed(2e-3, 1e-3, 0.05, 0.8) is None
    # Closer to beta needs more units.
    assert units_needed(6e-4, 1e-3, 0.05, 0.8) > units_needed(3e-4, 1e-3, 0.05, 0.8)


def test_critical_count_matches_exact_cdf():
    for n, beta, delta in [(50, 0.1, 0.05), (200, 0.05, 0.1), (30, 0.3, 0.01)]:
        k = critical_count(n, beta, delta)
        assert k == max([j for j in range(n + 1) if exact_cdf(j, n, beta) <= delta], default=-1)


def test_random_against_oracle():
    rng = random.Random(4)
    for _ in range(150):
        n = rng.randint(1, 120)
        beta = rng.choice([0.02, 0.05, 0.1, 0.2])
        delta = rng.choice([0.01, 0.05, 0.1])
        r = rng.random() * 2 * beta
        assert abs(certification_probability(n, r, beta, delta) - oracle(n, r, beta, delta)) < 1e-9


def test_units_needed_is_first_window_pass():
    r, beta, delta, power = 0.02, 0.05, 0.05, 0.8
    n = units_needed(r, beta, delta, power)
    assert certification_probability(n, r, beta, delta) >= power
    assert all(certification_probability(n + n * i // 50, r, beta, delta) >= power
               for i in range(51))
