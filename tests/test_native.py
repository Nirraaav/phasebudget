"""The native dynamic program answers to the Python one, exactly.

Exact means equal floats, not close ones: the native pass repeats the Python
arithmetic in the same order, so any difference is a bug, including a
different choice among tied optima.
"""

import random

import pytest

from phasebudget import _native
from phasebudget.dp import CoverageDP
from phasebudget.exact import solve_bruteforce
from phasebudget.fragments import FragmentMatrix
from phasebudget.simulate import simulate

pytestmark = pytest.mark.skipif(
    not _native.available(), reason="native library not built; run `make -C cpp`"
)


def _instances():
    rng = random.Random(20260925)
    for case in range(120):
        yield simulate(
            rng.randint(2, 40),
            coverage=rng.choice([2.0, 4.0, 6.0, 8.0]),
            read_span=rng.randint(2, 7),
            error_rate=rng.choice([0.01, 0.05, 0.15, 0.3]),
            gap_rate=rng.choice([0.0, 0.0, 0.2]),
            quality_spread=rng.choice([0.0, 2.0]),
            weight_model=rng.choice(["log_odds", "uniform"]),
            seed=case,
        ).matrix


def _both(matrix):
    return CoverageDP(matrix, backend="python"), CoverageDP(matrix, backend="native")


def test_native_solutions_equal_python_solutions_exactly():
    for matrix in _instances():
        python, native = _both(matrix)
        assert native.solve() == python.solve()


def test_native_boundary_costs_equal_python_boundary_costs_exactly():
    for matrix in _instances():
        python, native = _both(matrix)
        for boundary in range(matrix.n_variants - 1):
            assert native.boundary_costs(boundary) == python.boundary_costs(boundary)
        assert native.margins() == python.margins()


def test_native_solution_after_margins_is_the_same_solution():
    for matrix in list(_instances())[:30]:
        dp = CoverageDP(matrix, backend="native")
        margins = dp.margins()
        first = dp.solve()
        again = CoverageDP(matrix, backend="native").solve()
        assert first == again
        assert len(margins) == max(matrix.n_variants - 1, 0)


def test_every_checkpoint_spacing_gives_the_same_answer():
    for matrix in list(_instances())[:40]:
        dp = CoverageDP(matrix, backend="python")
        reference = (dp.solve(), [dp.boundary_costs(b) for b in range(matrix.n_variants - 1)])
        offsets, active, c0, c1 = dp.native_inputs()
        for interval in {1, 2, 3, 5, matrix.n_variants, matrix.n_variants + 4}:
            cost, haplotype, assignment, zero, one = _native.run(
                matrix.n_variants, offsets, active, c0, c1, len(matrix),
                want_margins=True, interval=interval,
            )
            assert (cost, haplotype, assignment) == (
                reference[0].cost, reference[0].haplotype, reference[0].assignment
            )
            assert list(zip(zero, one)) == reference[1]


def test_native_costs_match_the_exhaustive_oracle():
    rng = random.Random(7)
    for case in range(60):
        matrix = simulate(
            rng.randint(2, 12),
            coverage=rng.choice([3.0, 5.0]),
            read_span=rng.randint(2, 5),
            error_rate=0.1,
            seed=500 + case,
        ).matrix
        native = CoverageDP(matrix, backend="native").solve()
        assert native.cost == pytest.approx(solve_bruteforce(matrix).cost)
        assert matrix.cost(native.haplotype, native.assignment) == pytest.approx(native.cost)


def test_native_handles_empty_and_single_column_matrices():
    empty = FragmentMatrix([])
    assert CoverageDP(empty, backend="native").solve() == CoverageDP(empty, backend="python").solve()
    single = FragmentMatrix.from_rows(["0", "1", "1"])
    assert CoverageDP(single, backend="native").solve() == CoverageDP(single, backend="python").solve()
    assert CoverageDP(single, backend="native").margins() == ()


def test_disconnected_blocks_have_zero_margin_on_both_backends():
    matrix = FragmentMatrix.from_rows(["01--", "10--", "--01", "--11"])
    python, native = _both(matrix)
    assert native.margins() == python.margins()
    assert native.margins()[1].margin == 0.0


def test_column_sum_matches_the_builtin_sum_on_this_interpreter():
    import sys

    from phasebudget.dp import column_sum

    rng = random.Random(3)
    for _ in range(2000):
        values = [rng.choice([0.0, rng.uniform(0, 30), 1e-12, 1e12]) for _ in range(rng.randint(0, 9))]
        if sys.version_info >= (3, 12):
            assert column_sum(values) == sum(values)
        assert column_sum(values) == pytest.approx(sum(values))


def test_native_posterior_agrees_with_python_and_the_oracle():
    from phasebudget.exact import bruteforce_posterior_margins
    from phasebudget.posterior import boundary_log_sums

    for index, matrix in enumerate(_instances()):
        native = boundary_log_sums(matrix, backend="native")
        python = boundary_log_sums(matrix, backend="python")
        assert len(native) == len(python)
        for got, want in zip(native, python):
            assert got == pytest.approx(want, rel=1e-9, abs=1e-9)
        if matrix.n_variants <= 12 and index % 3 == 0:
            for got, want in zip(native, bruteforce_posterior_margins(matrix)):
                assert got == pytest.approx(want, rel=1e-9, abs=1e-9)


def test_native_posterior_is_the_same_at_every_checkpoint_spacing():
    for matrix in list(_instances())[:30]:
        dp = CoverageDP(matrix, backend="python")
        offsets, active, c0, c1 = dp.native_inputs()
        reference = _native.posterior(matrix.n_variants, offsets, active, c0, c1, len(matrix))
        for interval in (1, 2, 5, matrix.n_variants + 1):
            got = _native.posterior(
                matrix.n_variants, offsets, active, c0, c1, len(matrix), interval=interval
            )
            for a, b in zip(got, reference):
                assert a == pytest.approx(b, rel=1e-12, abs=1e-12)


def test_the_backend_is_validated():
    matrix = FragmentMatrix.from_rows(["01"])
    with pytest.raises(ValueError):
        CoverageDP(matrix, backend="gpu")
    assert CoverageDP(matrix).backend == "native"
    assert CoverageDP(matrix, backend="python").backend == "python"
