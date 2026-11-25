"""Contract checks for the genotype-error comparison and its two metrics."""

import importlib.util
import math
from pathlib import Path

from phasebudget import FragmentMatrix


def _benchmark_module():
    path = Path(__file__).parents[1] / "benchmarks" / "benchmark_genotype_errors.py"
    spec = importlib.util.spec_from_file_location("genotype_errors_benchmark", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_bipartition_error_ignores_the_global_flip():
    benchmark = _benchmark_module()
    assert benchmark.bipartition_error((0, 0, 1, 1), (0, 0, 1, 1)) == 0.0
    assert benchmark.bipartition_error((1, 1, 0, 0), (0, 0, 1, 1)) == 0.0
    assert benchmark.bipartition_error((0, 1, 1, 1), (0, 0, 1, 1)) == 0.25
    assert math.isnan(benchmark.bipartition_error((), ()))


def test_general_phase_marks_a_collapsed_column_undecided():
    """Both parts read allele 0 at variant 1, so the model declines to phase it."""
    benchmark = _benchmark_module()
    matrix = FragmentMatrix.from_rows(["00", "10"])
    haplotype, decided = benchmark.general_phase(matrix, (0, 1))
    assert decided[0] is True
    assert decided[1] is False
    assert haplotype[0] == 0


def test_general_phase_marks_an_uncovered_column_undecided():
    benchmark = _benchmark_module()
    matrix = FragmentMatrix.from_rows(["01-", "10-"])
    _, decided = benchmark.general_phase(matrix, (0, 1))
    assert decided[2] is False


def test_switch_error_scores_every_pair_without_a_decided_mask():
    benchmark = _benchmark_module()
    truth = (0, 1, 0, 1)
    predicted = (0, 1, 1, 0)
    error, fraction = benchmark.switch_error_on_sites(truth, predicted, (0, 1, 2, 3))
    assert fraction == 1.0
    assert error == 1 / 3


def test_a_decided_mask_drops_the_pairs_the_model_declined():
    benchmark = _benchmark_module()
    truth = (0, 1, 0, 1)
    predicted = (0, 1, 1, 0)
    decided = (True, True, False, True)
    error, fraction = benchmark.switch_error_on_sites(
        truth, predicted, (0, 1, 2, 3), decided
    )
    # Only the pair (0, 1) survives, and it agrees with the truth.
    assert fraction == 1 / 3
    assert error == 0.0


def test_a_model_that_phased_nothing_reports_a_zero_fraction():
    benchmark = _benchmark_module()
    error, fraction = benchmark.switch_error_on_sites(
        (0, 1, 0), (0, 1, 0), (0, 1, 2), (False, False, False)
    )
    assert math.isnan(error)
    assert fraction == 0.0


def test_the_comparison_emits_one_row_per_cell():
    benchmark = _benchmark_module()
    rows = benchmark.run(
        homozygous_rates=[0.0, 0.25],
        error_rates=[0.05],
        replicates=3,
        n_variants=8,
        coverage=4.0,
        span=4,
    )
    assert len(rows) == 2
    for row in rows:
        # The general model is a relaxation, so it can never cost more.
        assert row["general_cost"] <= row["complementary_cost"] + 1e-9
        assert 0.0 <= row["general_phased"] <= 1.0
