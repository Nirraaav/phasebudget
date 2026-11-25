"""Contract checks for the heterozygosity sweep."""

import importlib.util
import math
from pathlib import Path


def _benchmark_module():
    path = Path(__file__).parents[1] / "benchmarks" / "benchmark_heterozygosity.py"
    spec = importlib.util.spec_from_file_location("heterozygosity_benchmark", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_sweep_emits_one_row_per_cell():
    benchmark = _benchmark_module()
    rows = benchmark.run(
        heterozygosities=[5e-4, 1e-3],
        error_rates=[0.05, 0.2],
        replicates=2,
        n_variants=20,
        coverage=5.0,
        read_length=5000,
        max_active_reads=20,
    )
    assert len(rows) == 4
    assert {(row["heterozygosity"], row["error_rate"]) for row in rows} == {
        (5e-4, 0.05),
        (5e-4, 0.2),
        (1e-3, 0.05),
        (1e-3, 0.2),
    }
    for row in rows:
        assert row["instances"] + row["skipped"] == 2
        assert 0.0 <= row["phased"] <= 1.0
        assert 0.0 <= row["switch"] <= 1.0
        assert math.isnan(row["truth_optimal"]) or 0.0 <= row["truth_optimal"] <= 1.0


def test_denser_variants_put_more_of_them_on_each_read():
    """The axis the sweep exists to vary has to actually move."""
    benchmark = _benchmark_module()
    rows = benchmark.run(
        heterozygosities=[5e-4, 2e-3],
        error_rates=[0.05],
        replicates=4,
        n_variants=60,
        coverage=5.0,
        read_length=5000,
        max_active_reads=20,
    )
    sparse, dense = rows[0], rows[1]
    assert dense["per_read"] > sparse["per_read"]
    assert dense["blocks"] <= sparse["blocks"]
