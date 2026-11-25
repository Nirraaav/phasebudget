"""Contract checks for the adaptive sensitivity sweep."""

import importlib.util
from pathlib import Path


def _benchmark_module():
    path = Path(__file__).parents[1] / "benchmarks" / "benchmark_sensitivity.py"
    spec = importlib.util.spec_from_file_location("sensitivity_benchmark", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_sweep_emits_one_row_per_parameter_setting():
    benchmark = _benchmark_module()
    rows = benchmark.run(
        boosts=[1.5, 2.0],
        quantiles=[0.25],
        iteration_caps=[2, 4],
        error_rates=[0.1],
        caps=[3],
        replicates=2,
        n_variants=16,
        coverage=8.0,
        span=5,
        margin_source="subset",
    )
    assert len(rows) == 4
    assert {(row["boost"], row["iterations"]) for row in rows} == {
        (1.5, 2),
        (1.5, 4),
        (2.0, 2),
        (2.0, 4),
    }


def test_every_setting_shares_one_static_baseline():
    """A difference between two cells must not be a difference in instances."""
    benchmark = _benchmark_module()
    rows = benchmark.run(
        boosts=[1.5, 2.0, 4.0],
        quantiles=[0.25, 0.5],
        iteration_caps=[3],
        error_rates=[0.1],
        caps=[3],
        replicates=2,
        n_variants=16,
        coverage=8.0,
        span=5,
        margin_source="subset",
    )
    assert len({row["static"] for row in rows}) == 1
    for row in rows:
        assert row["delta"] == row["switch"] - row["static"]


def test_the_summary_counts_the_whole_grid():
    benchmark = _benchmark_module()
    rows = benchmark.run(
        boosts=[1.5, 2.0],
        quantiles=[0.25, 0.5],
        iteration_caps=[2],
        error_rates=[0.1, 0.2],
        caps=[3],
        replicates=2,
        n_variants=16,
        coverage=8.0,
        span=5,
        margin_source="subset",
    )
    summary = benchmark.summarize(rows)
    assert len(summary) == 2
    for entry in summary:
        assert entry["cells"] == 4
        assert entry["better"] + entry["worse"] <= entry["cells"]
        assert entry["best"] <= entry["median"] <= entry["worst"]
