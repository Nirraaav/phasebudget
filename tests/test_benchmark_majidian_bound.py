"""The counterexample to the step between Theorems 1 and 2 of Majidian et al. 2020."""

import importlib.util
from pathlib import Path

from phasebudget import dp
from phasebudget.exact import solve_bruteforce
from phasebudget.fragments import FragmentMatrix


def _benchmark_module():
    path = Path(__file__).parents[1] / "benchmarks" / "benchmark_majidian_bound.py"
    spec = importlib.util.spec_from_file_location("majidian_bound_benchmark", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_no_column_is_half_wrong_and_the_truth_still_loses():
    benchmark = _benchmark_module()
    example = benchmark.COUNTEREXAMPLE
    matrix = FragmentMatrix.from_rows(example["rows"])
    errors, called = benchmark.column_errors(
        matrix, example["truth"], example["origin"]
    )
    assert errors == [0, 1, 1, 1, 0]
    assert called == [1, 3, 4, 3, 2]
    assert all(2 * e < c for e, c in zip(errors, called) if c)
    assert not benchmark.has_bad_column(errors, called)
    assert matrix.cost(example["truth"]) == 3
    assert matrix.cost(example["better"]) == 2
    assert solve_bruteforce(matrix).cost <= 2
    assert dp.solve(matrix).cost == solve_bruteforce(matrix).cost


def test_every_read_in_the_counterexample_covers_two_variants():
    # Majidian et al. drop reads covering fewer than two SNPs before building
    # the fragment matrix, so the instance has to survive that filter.
    matrix = FragmentMatrix.from_rows(_benchmark_module().COUNTEREXAMPLE["rows"])
    assert all(
        sum(allele is not None for allele in fragment.alleles) >= 2
        for fragment in matrix.fragments
    )


def test_theorem2_counts_exactly_half_wrong_as_recoverable():
    benchmark = _benchmark_module()
    # One column of coverage 2 at error rate 0.5: n = 0 and n = 1 both count.
    assert benchmark.theorem2_probability([2], 0.5) == 0.75
    assert benchmark.theorem2_probability([0, 0], 0.3) == 1.0


def test_the_bound_benchmark_emits_one_row_per_cell():
    benchmark = _benchmark_module()
    rows = benchmark.run(
        n_variants=10,
        span=4,
        coverages=[4.0],
        error_rates=[0.05, 0.2],
        replicates=3,
        seed_base=9000,
    )
    assert len(rows) == 2
    for row in rows:
        assert row["gap_cases"] <= row["no_bad_column"]
        assert row["truth_optimal"] <= 1.0
