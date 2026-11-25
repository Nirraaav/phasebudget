"""Read selection under the called-allele coverage convention.

The repository counts coverage over a read's whole span, holes included. That
is the convention the dynamic program pays for and the one that keeps the
selection problem in P. This module implements the other convention, where a
read consumes coverage only where it calls an allele, because the comparison
between the two is a result and because the called convention is where the
open theory is. `docs/design/coverage-convention.md` states both, proves the
called convention NP-hard, and fixes which one the pipeline uses.

Nothing in `phasebudget` outside this module calls anything in it.
"""

from __future__ import annotations

from typing import Iterable, Optional, Sequence

from . import _pulp
from ._types import Fragment, SelectionResult
from ._validation import check_coverage, check_weights
from .fragments import FragmentMatrix
from .selection import read_weights, selection_objective


def called_positions(fragment: Fragment) -> tuple[int, ...]:
    """Variants where the read calls an allele, in order."""
    return tuple(position for position, _, _ in fragment.observed())


def called_runs(fragment: Fragment) -> int:
    """Number of maximal runs of consecutive called variants.

    One run is a hole-free read. The reduction in
    `docs/design/coverage-convention.md` section 5 is NP-hard already at three
    runs, and two runs is the open case.
    """
    positions = called_positions(fragment)
    if not positions:
        return 0
    return 1 + sum(
        1 for before, after in zip(positions, positions[1:]) if after != before + 1
    )


def unit_weights(matrix: FragmentMatrix) -> dict[int, float]:
    """One unit of weight per read, whatever it calls.

    `phasebudget.selection.read_weights` scores a read calling a single allele
    at zero, because such a read constrains no phase relation. That is a
    phasing policy and it is wrong for the reduction in
    `independent_set_instance`, where a degree-one vertex is worth as much as
    any other.
    """
    return {fragment.read_id: 1.0 for fragment in matrix}


def called_coverage(
    matrix: FragmentMatrix, selected: Optional[Iterable[int]] = None
) -> tuple[int, ...]:
    """Per-variant coverage counting called alleles only.

    This is at most `FragmentMatrix.coverage` entry by entry, which is the
    ordering in section 3 of the convention document.
    """
    subset = matrix if selected is None else matrix.subset(selected)
    counts = [0] * matrix.n_variants
    for fragment in subset:
        for position, _, _ in fragment.observed():
            counts[position] += 1
    return tuple(counts)


def max_called_coverage(
    matrix: FragmentMatrix, selected: Optional[Iterable[int]] = None
) -> int:
    return max(called_coverage(matrix, selected), default=0)


def is_feasible_called(
    matrix: FragmentMatrix, selected: Sequence[int], coverage: int
) -> bool:
    """Whether a selection respects the cap under the called convention."""
    check_coverage(coverage)
    return max_called_coverage(matrix, selected) <= coverage


def select_called_bruteforce(
    matrix: FragmentMatrix,
    coverage: int,
    weights: Optional[dict[int, float]] = None,
    *,
    max_fragments: int = 18,
) -> SelectionResult:
    """Enumerate every subset. The oracle for the called convention.

    Reads calling a single allele are kept. The phase-informative filter that
    `phasebudget.selection` applies is a policy of the phasing pipeline, and
    this function is used to check a combinatorial reduction in which a
    degree-one vertex is a legitimate read.
    """
    check_coverage(coverage)
    if weights is None:
        weights = read_weights(matrix)
    check_weights(weights, matrix.read_ids)
    if len(matrix) > max_fragments:
        raise ValueError(
            f"{len(matrix)} fragments exceeds the brute-force limit of {max_fragments}"
        )
    ids = matrix.read_ids
    calls = {
        fragment.read_id: called_positions(fragment) for fragment in matrix
    }
    best_value = 0.0
    best: tuple[int, ...] = ()
    for code in range(1 << len(ids)):
        candidate = tuple(ids[i] for i in range(len(ids)) if code >> i & 1)
        counts = [0] * matrix.n_variants
        feasible = True
        for read_id in candidate:
            for position in calls[read_id]:
                counts[position] += 1
                if counts[position] > coverage:
                    feasible = False
                    break
            if not feasible:
                break
        if not feasible:
            continue
        value = selection_objective(weights, candidate)
        if value > best_value:
            best_value, best = value, candidate
    return SelectionResult(
        selected=tuple(sorted(best)),
        objective=best_value,
        method="called_bruteforce",
        coverage=coverage,
    )


def _called_program(matrix: FragmentMatrix, coverage: int, weights, relaxed: bool):
    """Build the called-convention program in PuLP. Returns (problem, vars)."""
    import pulp

    sense = "Continuous" if relaxed else "Binary"
    problem = pulp.LpProblem("called_selection", pulp.LpMaximize)
    variables = {
        fragment.read_id: _pulp.variable(
            problem, f"x_{fragment.read_id}", low=0, up=1, cat=sense
        )
        for fragment in matrix
    }
    problem += pulp.lpSum(
        weights[read_id] * variable for read_id, variable in variables.items()
    )
    for position in range(matrix.n_variants):
        touching = [
            variables[fragment.read_id]
            for fragment in matrix
            if fragment.allele_at(position) is not None
        ]
        if touching:
            problem += pulp.lpSum(touching) <= coverage, f"cap_{position}"
    return problem, variables


def _solve(problem, solver=None) -> None:
    status = _pulp.solve(problem, solver)
    if status != "Optimal":
        raise RuntimeError(f"solver returned status {status}")


def select_called_ilp(
    matrix: FragmentMatrix,
    coverage: int,
    weights: Optional[dict[int, float]] = None,
    *,
    solver=None,
) -> SelectionResult:
    """Maximum-weight selection under the called convention, by 0-1 program.

    There is no flow formulation to fall back on: section 4 of
    `docs/design/coverage-convention.md` gives a three-read instance whose
    linear relaxation is strictly fractional. Raises `ImportError` without
    PuLP.
    """
    check_coverage(coverage)
    if weights is None:
        weights = read_weights(matrix)
    check_weights(weights, matrix.read_ids)
    if coverage == 0 or len(matrix) == 0 or matrix.n_variants == 0:
        return SelectionResult((), 0.0, "called_ilp", coverage)
    try:
        import pulp  # noqa: F401
    except ImportError as error:  # pragma: no cover - depends on the environment
        raise ImportError(
            "the called-convention program needs PuLP; install with "
            "pip install 'phasebudget[ilp]'"
        ) from error
    problem, variables = _called_program(matrix, coverage, weights, relaxed=False)
    _solve(problem, solver)
    selected = tuple(
        read_id
        for read_id, variable in sorted(variables.items())
        if round(variable.value() or 0.0) == 1
    )
    return SelectionResult(
        selected=selected,
        objective=selection_objective(weights, selected),
        method="called_ilp",
        coverage=coverage,
    )


def called_lp_bound(
    matrix: FragmentMatrix,
    coverage: int,
    weights: Optional[dict[int, float]] = None,
    *,
    solver=None,
) -> float:
    """Optimal value of the linear relaxation under the called convention.

    Strictly above the integer optimum on the triangle instance, which is what
    rules out the min-cost-flow argument. Raises `ImportError` without PuLP.
    """
    check_coverage(coverage)
    if weights is None:
        weights = read_weights(matrix)
    check_weights(weights, matrix.read_ids)
    if coverage == 0 or len(matrix) == 0 or matrix.n_variants == 0:
        return 0.0
    try:
        import pulp
    except ImportError as error:  # pragma: no cover - depends on the environment
        raise ImportError(
            "the called-convention relaxation needs PuLP; install with "
            "pip install 'phasebudget[ilp]'"
        ) from error
    problem, _ = _called_program(matrix, coverage, weights, relaxed=True)
    _solve(problem, solver)
    return float(pulp.value(problem.objective))


def independent_set_instance(
    n_vertices: int, edges: Sequence[tuple[int, int]]
) -> FragmentMatrix:
    """The reduction of maximum independent set to called-convention selection.

    One variant per edge, one read per vertex, and a read calls the variants of
    the edges incident to its vertex. At cap 1 the constraint at the variant of
    edge `{u, v}` is exactly that `u` and `v` are not both selected, so a
    feasible selection under unit weights is an independent set of the same
    size. The argument is section 5 of `docs/design/coverage-convention.md`.

    Isolated vertices are rejected rather than handled. They belong to every
    maximum independent set and can be removed before the reduction, and a read
    calling nothing is not a fragment.
    """
    if n_vertices < 0:
        raise ValueError("n_vertices must be non-negative")
    incident: list[list[int]] = [[] for _ in range(n_vertices)]
    for index, (left, right) in enumerate(edges):
        if not 0 <= left < n_vertices or not 0 <= right < n_vertices:
            raise ValueError(f"edge {index} names a vertex outside the graph")
        if left == right:
            raise ValueError(f"edge {index} is a loop")
        incident[left].append(index)
        incident[right].append(index)
    fragments = []
    for vertex in range(n_vertices):
        positions = incident[vertex]
        if not positions:
            raise ValueError(f"vertex {vertex} is isolated")
        start, stop = positions[0], positions[-1] + 1
        called = set(positions)
        alleles = [0 if j in called else None for j in range(start, stop)]
        weights = [1.0 if j in called else 0.0 for j in range(start, stop)]
        fragments.append(
            Fragment(
                read_id=vertex,
                start=start,
                alleles=tuple(alleles),
                weights=tuple(weights),
            )
        )
    return FragmentMatrix(fragments, len(edges))


def max_independent_set_bruteforce(
    n_vertices: int, edges: Sequence[tuple[int, int]]
) -> tuple[int, ...]:
    """A maximum independent set by enumeration, for checking the reduction."""
    best: tuple[int, ...] = ()
    for code in range(1 << n_vertices):
        chosen = tuple(v for v in range(n_vertices) if code >> v & 1)
        member = set(chosen)
        if any(left in member and right in member for left, right in edges):
            continue
        if len(chosen) > len(best):
            best = chosen
    return best


def two_union_instance(
    segments: Sequence[tuple[tuple[int, int], tuple[int, int]]]
) -> FragmentMatrix:
    """Embed a 2-union instance as reads with exactly one hole each.

    Vertex `v` of a 2-union graph is a pair of integer segments, `[a, b)` on a
    first line and `[c, d)` on a second, and two vertices conflict when their
    segments meet on either line. Laying the second line to the right of the
    first turns each vertex into one read that calls `[a, b)`, skips to the
    offset copy of `[c, d)`, and calls that. Two reads then share a called
    variant exactly when the two vertices conflict, so a cap-1 selection is an
    independent set.

    Bar-Yehuda, Halldórsson, Naor, Shachnai and Shapira, "Scheduling split
    intervals", SIAM J. Comput. 36(1), 2006, Corollary 2.3, prove maximum
    independent set APX-hard on unweighted 2-union graphs whose segments all
    have length 2. Through this embedding the called convention is APX-hard
    with one hole per read, which settles OQ-11.
    """
    if not segments:
        return FragmentMatrix((), 0)
    for index, ((a, b), (c, d)) in enumerate(segments):
        if not (0 <= a < b and 0 <= c < d):
            raise ValueError(f"vertex {index} has an empty or negative segment")
    offset = max(b for (_, b), _ in segments) + 1
    fragments = []
    for vertex, ((a, b), (c, d)) in enumerate(segments):
        called = set(range(a, b)) | set(range(offset + c, offset + d))
        start, stop = a, offset + d
        fragments.append(
            Fragment(
                read_id=vertex,
                start=start,
                alleles=tuple(0 if j in called else None for j in range(start, stop)),
                weights=tuple(1.0 if j in called else 0.0 for j in range(start, stop)),
            )
        )
    width = offset + max(d for _, (_, d) in segments)
    return FragmentMatrix(fragments, width)
