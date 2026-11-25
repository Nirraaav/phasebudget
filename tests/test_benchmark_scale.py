"""Small contract checks for the milestone-1 to milestone-4 benchmark scripts."""

import importlib.util
import math
from pathlib import Path


def _load(name):
    path = Path(__file__).parents[1] / "benchmarks" / name
    spec = importlib.util.spec_from_file_location(name[:-3], path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


NONE = {"systematic_rate": 0.0, "systematic_strength": 0.5, "burst_rate": 0.0,
        "burst_length": 5, "burst_error": 0.35}


def test_selection_at_scale_scores_every_selector_and_certifies_the_optimum():
    bench = _load("benchmark_selection_at_scale.py")
    profile, cap, replicate, reads, measured = bench.unit(("ont", 5, 0, 120, 20.0, 0.001, 2.0, NONE))
    assert set(measured) == set(bench.ALL)
    best = measured["optimal"]["objective"]
    for name in bench.SELECTORS:
        assert measured[name]["objective"] <= best + 1e-9
    assert math.isnan(measured["merge_optimal"]["objective"])
    assert 0.0 <= measured["merge_optimal"]["mixed"] <= 1.0


def test_adaptive_at_scale_pairs_the_arms_on_one_instance():
    bench = _load("benchmark_adaptive_at_scale.py")
    arguments = ("noisy", 3, 0, 120, 20.0, 0.001, 2.0, "posterior", 0.05, "stitched", 20,
                 True, 3, 2.0, False, 8000, NONE, {"p": 0.15, "threshold": 1e6})
    profile, cap, exact, adaptive, control, blocks = bench.unit(arguments)
    assert set(blocks) == {"exact", "adaptive", "control"}
    assert (profile, cap) == ("noisy", 3)
    for rate in (exact, adaptive, control):
        assert 0.0 <= rate <= 1.0


def test_phase_diagram_cell_reports_every_column():
    bench = _load("benchmark_phase_diagram.py")
    row = bench.cell((6, 0.1, 0.001, 2.0, 3, 20))
    assert row["replicates"] == 3
    assert 0.0 <= row["truth_optimal"] <= 1.0
    assert row["gap"] >= -1e-9
    assert {"margin_info", "posterior_info", "skipped"} <= set(row)


def test_the_score_rule_diagnostic_still_finds_the_line_it_replaces():
    bench = _load("diagnose_whatshap_score_rule.py")
    paper_rule = bench._paper_rule_selector()
    profile, cap, out = bench.unit(("ont", 3, 0, 150))
    assert set(out) == {"shipped", "paper rule"}
    for record in out.values():
        assert 0.0 < record["ratio"] <= 1.0 + 1e-9
        assert record["blocks"] >= 1
    assert callable(paper_rule)


def test_resources_skips_what_would_not_fit_and_measures_a_step(tmp_path):
    import os

    import pytest

    bench = _load("benchmark_resources.py")
    assert not bench.feasible("native_full_solve", 8000, 20, 6.0, 10, 2000)
    assert bench.feasible("native_sqrt_solve", 8000, 20, 6.0, 10, 2000)
    assert not bench.feasible("python_solve", 500, 15, 6.0, 10, 2000)
    if not os.path.exists("/proc/self/clear_refs"):
        pytest.skip("peak-memory reset needs Linux")
    path = str(tmp_path / "instance.pickle")
    bench.prepare(80, 20.0, 4, 1, path)
    record = bench.measure("python_solve", path, 4)
    assert record["active"] <= 4
    assert record["total_s"] >= 0.0 and record["peak_mb"] >= 0.0


def test_the_simulated_trio_is_mendelian_and_the_unit_runs():
    import pytest

    bench = _load("benchmark_pedigree.py")
    positions, genomes, reads = bench.simulate_trio(
        150, depth=10, read_length=10_000, error_rate=0.05, spread=2.0, seed=3)
    assert positions == sorted(positions)
    father, mother, child = genomes["father"], genomes["mother"], genomes["child"]
    for j in range(len(positions)):
        assert child[0][j] in (father[0][j], father[1][j])
        assert child[1][j] in (mother[0][j], mother[1][j])
        assert any(genomes[m][0][j] != genomes[m][1][j] for m in bench.MEMBERS)
    for member in bench.MEMBERS:
        assert reads[member].n_variants == len(positions)
    pytest.importorskip("whatshap")
    profile, cap, out = bench.unit(("ont", 3, 0, 150, 10.0, 2.0, 30000))
    assert set(out) == set(bench.SELECTORS)
    for record in out.values():
        for switches, pairs in record.values():
            assert 0 <= switches <= pairs
