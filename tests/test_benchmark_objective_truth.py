"""Small contract checks for the objective-truth benchmark grid."""

import importlib.util
from pathlib import Path


def _benchmark_module():
    path = Path(__file__).parents[1] / "benchmarks" / "benchmark_objective_truth.py"
    spec = importlib.util.spec_from_file_location("objective_truth_benchmark", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_objective_truth_benchmark_emits_one_row_per_grid_cell():
    benchmark = _benchmark_module()
    rows = benchmark.run(
        error_rates=[0.05, 0.2],
        replicates=2,
        n_variants=12,
        coverages=[4.0, 6.0],
        spans=[4],
    )
    assert len(rows) == 4
    assert {(row["coverage"], row["span"]) for row in rows} == {
        (4.0, 4),
        (6.0, 4),
    }
    assert all(0.0 <= row["truth_optimal"] <= 1.0 for row in rows)
