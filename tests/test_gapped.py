"""The called-allele coverage convention, and why it is not the one used.

`docs/design/coverage-convention.md` makes four claims that are checkable:
called coverage never exceeds span coverage, the relaxation is fractional on
the triangle, the reduction from independent set is exact, and the span
optimum is a lower bound on the called optimum. Each has a test here.
"""

from __future__ import annotations

import itertools
import random

import pytest

from phasebudget import FragmentMatrix
from phasebudget.gapped import (
    called_coverage,
    called_lp_bound,
    called_runs,
    independent_set_instance,
    is_feasible_called,
    max_called_coverage,
    max_independent_set_bruteforce,
    select_called_bruteforce,
    select_called_ilp,
    two_union_instance,
    unit_weights,
)
from phasebudget.selection import is_feasible, read_weights, select_optimal

from conftest import random_matrix

TRIANGLE = FragmentMatrix.from_rows(["11-", "1-1", "-11"])


def graphs_up_to(n_vertices: int):
    """Every simple graph on `n_vertices` labelled vertices, no isolated ones."""
    pairs = list(itertools.combinations(range(n_vertices), 2))
    for mask in range(1 << len(pairs)):
        edges = [pair for index, pair in enumerate(pairs) if mask >> index & 1]
        touched = {vertex for edge in edges for vertex in edge}
        if len(touched) != n_vertices:
            continue
        yield edges


def test_called_runs_counts_maximal_runs():
    matrix = FragmentMatrix.from_rows(["111--", "1-1--", "1--11", "-1---"])
    assert [called_runs(fragment) for fragment in matrix] == [1, 2, 2, 1]


def test_called_coverage_never_exceeds_span_coverage():
    generator = random.Random(31)
    for _ in range(300):
        matrix = random_matrix(generator)
        called = called_coverage(matrix)
        span = matrix.coverage()
        assert len(called) == len(span)
        assert all(a <= b for a, b in zip(called, span))


def test_called_coverage_equals_span_coverage_without_holes():
    matrix = FragmentMatrix.from_rows(["11--", "-111", "--11"])
    assert called_coverage(matrix) == matrix.coverage()


def test_the_hole_is_what_separates_the_two_conventions():
    assert TRIANGLE.coverage() == (2, 3, 2)
    assert called_coverage(TRIANGLE) == (2, 2, 2)


def test_span_feasible_selections_are_called_feasible():
    generator = random.Random(32)
    for _ in range(300):
        matrix = random_matrix(generator)
        for coverage in range(0, 4):
            selected = select_optimal(matrix, coverage).selected
            assert is_feasible(matrix, selected, coverage)
            assert is_feasible_called(matrix, selected, coverage)


def test_called_feasibility_can_admit_what_the_span_cap_refuses():
    selected = (0, 1, 2)
    assert max_called_coverage(TRIANGLE, selected) == 2
    assert TRIANGLE.subset(selected).max_coverage() == 3
    assert is_feasible_called(TRIANGLE, selected, 2)
    assert not is_feasible(TRIANGLE, selected, 2)


def test_span_optimum_is_a_lower_bound_on_the_called_optimum():
    generator = random.Random(33)
    for _ in range(200):
        matrix = random_matrix(generator, max_variants=6, max_reads=6)
        weights = read_weights(matrix)
        for coverage in range(1, 4):
            span = select_optimal(matrix, coverage, weights)
            called = select_called_bruteforce(matrix, coverage, weights)
            assert span.objective <= called.objective + 1e-9


def test_the_span_convention_is_sometimes_strictly_worse():
    weights = read_weights(TRIANGLE)
    span = select_optimal(TRIANGLE, 2, weights)
    called = select_called_bruteforce(TRIANGLE, 2, weights)
    assert span.objective == pytest.approx(4.0)
    assert called.objective == pytest.approx(6.0)


def test_brute_force_selections_respect_the_called_cap():
    generator = random.Random(34)
    for _ in range(200):
        matrix = random_matrix(generator, max_variants=6, max_reads=6)
        for coverage in range(0, 3):
            result = select_called_bruteforce(matrix, coverage)
            assert is_feasible_called(matrix, result.selected, coverage)


def test_brute_force_refuses_an_instance_it_cannot_enumerate():
    matrix = random_matrix(random.Random(35), max_variants=6, max_reads=6)
    with pytest.raises(ValueError, match="brute-force limit"):
        select_called_bruteforce(matrix, 2, max_fragments=0)


def test_the_reduction_is_exact_on_every_small_graph():
    for n_vertices in (2, 3, 4, 5):
        for edges in graphs_up_to(n_vertices):
            matrix = independent_set_instance(n_vertices, edges)
            result = select_called_bruteforce(matrix, 1, unit_weights(matrix))
            independent = max_independent_set_bruteforce(n_vertices, edges)
            assert result.objective == pytest.approx(float(len(independent)))
            assert len(result.selected) == len(independent)
            chosen = set(result.selected)
            assert not any(
                left in chosen and right in chosen for left, right in edges
            )


def test_the_reduction_maps_the_triangle_to_the_triangle():
    matrix = independent_set_instance(3, [(0, 1), (0, 2), (1, 2)])
    assert matrix.to_rows() == ("00-", "0-0", "-00")


def test_the_reduction_reaches_three_runs_at_degree_three():
    # Vertex 0 is incident to edges 0, 2, and 4, which are not consecutive, so
    # its read has two holes. Three runs is where section 5 bites.
    matrix = independent_set_instance(4, [(0, 1), (2, 3), (0, 2), (1, 3), (0, 3)])
    runs = [called_runs(fragment) for fragment in matrix]
    assert runs[0] == 3
    assert max(runs) == 3


def test_the_reduction_rejects_isolated_vertices_and_loops():
    with pytest.raises(ValueError, match="isolated"):
        independent_set_instance(3, [(0, 1)])
    with pytest.raises(ValueError, match="loop"):
        independent_set_instance(2, [(0, 0), (0, 1)])


def test_the_triangle_relaxation_is_strictly_fractional():
    pytest.importorskip("pulp")
    weights = unit_weights(TRIANGLE)
    integral = select_called_bruteforce(TRIANGLE, 1, weights)
    assert integral.objective == pytest.approx(1.0)
    assert called_lp_bound(TRIANGLE, 1, weights) == pytest.approx(1.5)


def test_the_span_relaxation_is_integral_on_the_same_instance():
    pytest.importorskip("pulp")
    weights = unit_weights(TRIANGLE)
    assert select_optimal(TRIANGLE, 1, weights).objective == pytest.approx(1.0)


def test_the_called_program_agrees_with_enumeration():
    pytest.importorskip("pulp")
    generator = random.Random(36)
    for _ in range(25):
        matrix = random_matrix(generator, max_variants=6, max_reads=6)
        weights = read_weights(matrix)
        for coverage in (1, 2):
            solved = select_called_ilp(matrix, coverage, weights)
            reference = select_called_bruteforce(matrix, coverage, weights)
            assert solved.objective == pytest.approx(reference.objective)
            assert is_feasible_called(matrix, solved.selected, coverage)


def test_the_called_program_is_empty_at_cap_zero():
    pytest.importorskip("pulp")
    assert select_called_ilp(TRIANGLE, 0).selected == ()
    assert called_lp_bound(TRIANGLE, 0) == 0.0


def _segments_meet(left, right):
    return max(left[0], right[0]) < min(left[1], right[1])


def test_the_two_union_embedding_has_one_hole_run_per_read():
    matrix = two_union_instance([((0, 2), (3, 5)), ((1, 3), (0, 2))])
    assert [called_runs(fragment) for fragment in matrix] == [2, 2]


def test_the_two_union_embedding_preserves_every_conflict():
    """Reads share a called variant exactly when the vertices conflict."""
    generator = random.Random(37)
    for _ in range(300):
        segments = []
        for _ in range(generator.randint(1, 7)):
            a = generator.randint(0, 6)
            c = generator.randint(0, 6)
            segments.append(((a, a + 2), (c, c + 2)))
        matrix = two_union_instance(segments)
        calls = [set(p for p, _, _ in fragment.observed()) for fragment in matrix]
        for u in range(len(segments)):
            for v in range(u + 1, len(segments)):
                conflict = _segments_meet(segments[u][0], segments[v][0]) or (
                    _segments_meet(segments[u][1], segments[v][1])
                )
                assert bool(calls[u] & calls[v]) == conflict


def test_a_cap_one_selection_on_the_embedding_is_a_maximum_independent_set():
    generator = random.Random(38)
    for _ in range(60):
        segments = []
        for _ in range(generator.randint(1, 7)):
            a = generator.randint(0, 5)
            c = generator.randint(0, 5)
            segments.append(((a, a + 2), (c, c + 2)))
        edges = [
            (u, v)
            for u in range(len(segments))
            for v in range(u + 1, len(segments))
            if _segments_meet(segments[u][0], segments[v][0])
            or _segments_meet(segments[u][1], segments[v][1])
        ]
        matrix = two_union_instance(segments)
        result = select_called_bruteforce(matrix, 1, unit_weights(matrix))
        independent = max_independent_set_bruteforce(len(segments), edges)
        assert result.objective == pytest.approx(float(len(independent)))
