"""The window bootstrap registered for the real-data run."""

import random

import pytest

from phasebudget import Phasing
from phasebudget.bootstrap import rate_interval, ratio_interval, window_counts
from phasebudget.dp import CoverageDP
from phasebudget.metrics import comparable_pairs, phasing_from_solution, switch_positions
from phasebudget.selection import read_weights, select_optimal, select_whatshap
from phasebudget.simulate import simulate_genomic


def _phased(instance, selector, cap):
    matrix = instance.matrix
    chosen = selector(matrix, cap, read_weights(matrix)).selected
    subset = matrix.subset(chosen)
    return phasing_from_solution(subset, CoverageDP(subset).solve())


@pytest.fixture(scope="module")
def instance():
    return simulate_genomic(300, coverage=20.0, heterozygosity=0.001, quality_spread=2.0,
                            read_length=10_000, error_rate=0.1, seed=3)


def test_window_counts_sum_to_the_metrics(instance):
    phasing = _phased(instance, select_optimal, 3)
    for window in (1, 50_000, 1_000_000, 10_000_000):
        counts = window_counts(instance.haplotype, phasing, instance.positions, window)
        assert sum(s for s, _ in counts) == len(switch_positions(instance.haplotype, phasing))
        assert sum(p for _, p in counts) == comparable_pairs(phasing, instance.haplotype)


def test_window_counts_respect_a_partial_truth(instance):
    phasing = _phased(instance, select_optimal, 3)
    rng = random.Random(1)
    truth = [None if rng.random() < 0.2 else a for a in instance.haplotype]
    counts = window_counts(truth, phasing, instance.positions, 100_000)
    assert sum(s for s, _ in counts) == len(switch_positions(truth, phasing))
    assert sum(p for _, p in counts) == comparable_pairs(phasing, truth)


def test_a_hand_example():
    # Variants at 0, 10, 20, 30 in one block; the relation between the last
    # two is wrong, so the one switch is recorded at position 30.
    phasing = Phasing(alleles=(0, 1, 0, 0), blocks=(0, 0, 0, 0))
    counts = window_counts((0, 1, 0, 1), phasing, (0, 10, 20, 30), 20)
    assert counts == ((0, 1), (1, 2))


def test_the_rate_interval_contains_the_estimate_and_is_reproducible(instance):
    phasing = _phased(instance, select_optimal, 3)
    counts = window_counts(instance.haplotype, phasing, instance.positions, 50_000)
    first = rate_interval(counts, resamples=500, seed=7)
    assert first == rate_interval(counts, resamples=500, seed=7)
    assert first.low <= first.estimate <= first.high
    assert first.estimate == pytest.approx(
        len(switch_positions(instance.haplotype, phasing))
        / comparable_pairs(phasing, instance.haplotype))


def test_a_pipeline_against_itself_has_ratio_one(instance):
    counts = window_counts(instance.haplotype, _phased(instance, select_optimal, 3),
                           instance.positions, 50_000)
    interval = ratio_interval(counts, counts, resamples=300)
    assert (interval.estimate, interval.low, interval.high) == (1.0, 1.0, 1.0)


def test_the_ratio_is_computed_on_jointly_resampled_windows(instance):
    shipped = window_counts(instance.haplotype,
                            _phased(instance, lambda m, c, w: select_whatshap(m, c, w), 3),
                            instance.positions, 50_000)
    optimal = window_counts(instance.haplotype, _phased(instance, select_optimal, 3),
                            instance.positions, 50_000)
    interval = ratio_interval(shipped, optimal, resamples=500)
    total = lambda counts, i: sum(c[i] for c in counts)
    assert interval.estimate == pytest.approx(
        (total(shipped, 0) / total(shipped, 1)) / (total(optimal, 0) / total(optimal, 1)))
    assert interval.low <= interval.estimate <= interval.high


def test_undefined_statistics_are_rejected():
    with pytest.raises(ValueError):
        rate_interval(((0, 0), (0, 0)))
    with pytest.raises(ValueError):
        window_counts((0, 1), Phasing(alleles=(0, 1), blocks=(0, 0)), (0, 10), 0)


def test_a_custom_statistic_reproduces_the_ratio_interval(instance):
    from phasebudget.bootstrap import statistic_interval

    shipped = window_counts(instance.haplotype,
                            _phased(instance, lambda m, c, w: select_whatshap(m, c, w), 3),
                            instance.positions, 50_000)
    optimal = window_counts(instance.haplotype, _phased(instance, select_optimal, 3),
                            instance.positions, 50_000)

    def ratio(picks):
        a = sum(shipped[k][1] for k in picks)
        b = sum(optimal[k][1] for k in picks)
        top = sum(shipped[k][0] for k in picks)
        bottom = sum(optimal[k][0] for k in picks)
        if not a or not b or not bottom:
            return None
        return (top / a) / (bottom / b)

    assert statistic_interval(ratio, len(shipped), resamples=300) == \
        ratio_interval(shipped, optimal, resamples=300)
