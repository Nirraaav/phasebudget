"""The latent class fit of `benchmarks/r1_truth_audit.py` (R1) against simulated raters."""

import importlib.util
import random
import sys
from collections import Counter
from pathlib import Path

BENCH = Path(__file__).resolve().parent.parent / "benchmarks"


def _module():
    sys.path.insert(0, str(BENCH))
    spec = importlib.util.spec_from_file_location("r1_truth_audit", BENCH / "r1_truth_audit.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _simulate(n, pi, errors, seed):
    rng = random.Random(seed)
    counts = Counter()
    for _ in range(n):
        z = 1 if rng.random() < pi else 0
        counts[tuple(1 - z if rng.random() < e else z for e in errors)] += 1
    return counts


def test_fit_recovers_known_error_rates():
    m = _module()
    errors = (2e-3, 1e-2, 3e-2, 1e-1)
    pi, fitted, _ = m.fit(_simulate(400_000, 0.55, errors, 1))
    assert abs(pi - 0.55) < 0.005
    for e, f in zip(errors, fitted):
        assert abs(f - e) < 0.2 * e + 2e-4, (e, f)


def test_fit_with_one_noisy_rater_puts_the_noise_on_it():
    m = _module()
    pi, fitted, _ = m.fit(_simulate(200_000, 0.5, (0.0, 0.0, 0.0, 0.05), 2))
    assert fitted[3] > 0.04 and max(fitted[:3]) < 1e-4


def test_posterior_follows_the_majority_of_good_raters():
    m = _module()
    err = (1e-3, 1e-3, 1e-3, 0.2)
    assert m.posterior((1, 1, 1, 0), 0.5, err) > 0.999
    assert m.posterior((0, 0, 1, 1), 0.5, err) < 0.5


def test_relation_needs_one_block_for_a_phaser():
    m = _module()
    assert m.relation(0, 0) == 1 and m.relation(0, 1) == 0
    assert m.relation([0, 7], [1, 7]) == 0 and m.relation([1, 7], [1, 7]) == 1
    assert m.relation([0, 7], [0, 8]) is None and m.relation(None, 1) is None


def test_g2_is_small_when_the_model_holds():
    m = _module()
    counts = _simulate(200_000, 0.5, (1e-2, 2e-2, 3e-2, 4e-2), 3)
    pi, err, _ = m.fit(counts)
    stat, df = m.g2(counts, pi, err)
    assert df == 10 and stat < 40
