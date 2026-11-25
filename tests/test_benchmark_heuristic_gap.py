"""Contract checks for the headline optimality-gap measurement.

The claims the table makes are that the ratio column is a certified optimality
gap, that the four levels are measured separately, and that the two wMEC
columns are different quantities. Each is checkable on a tiny grid.
"""

import importlib.util
from pathlib import Path

from phasebudget.simulate import simulate_connected


def _benchmark_module():
    path = Path(__file__).parents[1] / "benchmarks" / "benchmark_heuristic_gap.py"
    spec = importlib.util.spec_from_file_location("heuristic_gap_benchmark", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_grid_emits_one_row_per_selector_per_cell():
    benchmark = _benchmark_module()
    rows = benchmark.run(
        caps=[2, 3],
        error_rates=[0.05],
        replicates=2,
        n_variants=14,
        coverage=8.0,
        span=4,
    )
    assert len(rows) == 2 * len(benchmark.SELECTORS)
    assert {row["cap"] for row in rows} == {2, 3}


def test_no_selector_beats_the_certified_optimum():
    benchmark = _benchmark_module()
    rows = benchmark.run(
        caps=[2, 4],
        error_rates=[0.05, 0.2],
        replicates=3,
        n_variants=16,
        coverage=10.0,
        span=5,
    )
    for row in rows:
        assert row["ratio"] <= 1.0 + 1e-9
        assert row["ratio_worst"] <= row["ratio"] + 1e-9
        if row["selector"] == "optimal":
            assert row["ratio"] == 1.0
            assert row["ratio_worst"] == 1.0


def test_the_worst_case_is_a_worst_case():
    benchmark = _benchmark_module()
    rows = benchmark.run(
        caps=[3],
        error_rates=[0.1],
        replicates=4,
        n_variants=16,
        coverage=10.0,
        span=5,
    )
    for row in rows:
        assert row["ratio_worst"] <= row["ratio"] + 1e-9
        assert row["switch_worst"] >= row["switch"] - 1e-9


def test_the_two_mec_columns_are_different_quantities():
    """Cost on the retained reads is not cost on all of them."""
    benchmark = _benchmark_module()
    instance = simulate_connected(
        20, coverage=12.0, read_span=5, error_rate=0.1, seed=7000
    )
    record = benchmark.measure(instance, 3, "random", 0)
    assert record["mec_kept"] <= record["mec_full"]
    assert record["reads"] < len(instance.matrix)


def test_selection_and_solve_time_are_reported_apart():
    benchmark = _benchmark_module()
    instance = simulate_connected(
        20, coverage=12.0, read_span=5, error_rate=0.05, seed=7013
    )
    record = benchmark.measure(instance, 3, "fischer_marschall", 0)
    assert record["select_ms"] > 0.0
    assert record["solve_ms"] > 0.0
    assert "select_ms" in record and "solve_ms" in record
