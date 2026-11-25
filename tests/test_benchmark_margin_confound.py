"""Contract checks for the margin-source benchmark, and its rank statistics."""

import importlib.util
import math
from pathlib import Path


def _benchmark_module():
    path = Path(__file__).parents[1] / "benchmarks" / "benchmark_margin_confound.py"
    spec = importlib.util.spec_from_file_location("margin_confound_benchmark", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_spearman_is_one_on_a_monotone_pair():
    benchmark = _benchmark_module()
    assert benchmark.spearman([1.0, 2.0, 3.0, 4.0], [5.0, 6.0, 7.0, 8.0]) == 1.0
    assert benchmark.spearman([1.0, 2.0, 3.0, 4.0], [8.0, 7.0, 6.0, 5.0]) == -1.0


def test_spearman_matches_the_hand_computed_value():
    # d = (-1, 1, -1, 1, 0), so rho = 1 - 6*4/(5*24) = 0.8.
    benchmark = _benchmark_module()
    value = benchmark.spearman([1.0, 2.0, 3.0, 4.0, 5.0], [2.0, 1.0, 4.0, 3.0, 5.0])
    assert value == 0.8


def test_spearman_is_undefined_on_a_constant_vector():
    benchmark = _benchmark_module()
    assert math.isnan(benchmark.spearman([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]))
    assert math.isnan(benchmark.spearman([1.0, 2.0], [1.0, 2.0]))


def test_tied_values_share_an_average_rank():
    benchmark = _benchmark_module()
    assert benchmark._ranks([5.0, 1.0, 1.0, 3.0]) == [4.0, 1.5, 1.5, 3.0]


def test_rank_probability_is_one_half_when_the_signals_coincide():
    benchmark = _benchmark_module()
    assert benchmark.rank_probability([1.0], [1.0]) == 0.5
    assert benchmark.rank_probability([0.0], [1.0]) == 1.0
    assert math.isnan(benchmark.rank_probability([], [1.0]))


def test_the_confound_benchmark_emits_one_row_per_cell():
    benchmark = _benchmark_module()
    rows = benchmark.run(
        error_rates=[0.05, 0.2],
        caps=[2, 3],
        replicates=2,
        n_variants=12,
        coverage=6.0,
        span=4,
        iterations=2,
    )
    assert len(rows) == 4
    assert {(row["error_rate"], row["cap"]) for row in rows} == {
        (0.05, 2),
        (0.05, 3),
        (0.2, 2),
        (0.2, 3),
    }
    for row in rows:
        assert row["instances"] + row["skipped"] == 2
        assert 0.0 <= row["hidden"] <= 1.0 or math.isnan(row["hidden"])
        assert 0.0 <= row["switch_subset"] <= 1.0
        assert 0.0 <= row["switch_full"] <= 1.0
