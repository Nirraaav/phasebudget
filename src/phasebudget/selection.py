"""Coverage-constrained read selection.

Each read `r` occupies an interval `[l_r, u_r]` of variant indices. Selecting
`x_r in {0,1}` subject to

    sum over r with j in [l_r, u_r] of x_r <= c    for every variant j

and maximizing `sum_r w_r x_r` is an integer program whose constraint matrix
has the consecutive-ones property, so its linear relaxation is integral. The
same structure gives a min-cost-flow formulation, which is what
`select_optimal` solves: a path of `c`-capacity arcs along the variants, one
unit-capacity arc of cost `-w_r` per read, and `c` units pushed from the first
node to the last. The argument is written out in
`docs/design/selection-integrality.md`.
"""

from __future__ import annotations

import heapq
import random
from typing import Callable, Iterable, Optional, Sequence

from ._types import Fragment, SelectionResult
from ._validation import check_coverage, check_weights
from .fragments import FragmentMatrix

INF = float("inf")


def is_phase_informative(fragment: Fragment) -> bool:
    """Whether a read constrains any phase relation at all.

    A read calling one allele can always be assigned to whichever haplotype
    matches it, so it costs nothing under wMEC and says nothing about phase.
    Fischer and Marschall 2016 discard these reads in the first line of their
    procedure, and the reason is the same here: such a read consumes coverage
    and buys nothing.
    """
    return fragment.n_observed >= 2


def total_weight_score(fragment: Fragment) -> float:
    """Sum of the confidence weights of a read's called alleles."""
    return fragment.total_weight if is_phase_informative(fragment) else 0.0


def variant_count_score(fragment: Fragment) -> float:
    """Number of called alleles, which is the unweighted information content."""
    return float(fragment.n_observed) if is_phase_informative(fragment) else 0.0


def span_score(fragment: Fragment) -> float:
    """Length of the covered interval, gaps included."""
    return float(len(fragment)) if is_phase_informative(fragment) else 0.0


def holes(fragment: Fragment) -> int:
    """Uncalled variants strictly inside a read's span."""
    return len(fragment) - fragment.n_observed


def read_weights(
    matrix: FragmentMatrix,
    score: Callable[[Fragment], float] = total_weight_score,
) -> dict[int, float]:
    """Score every read in a matrix, keyed by read id."""
    return {fragment.read_id: float(score(fragment)) for fragment in matrix}


def selection_objective(weights: dict[int, float], selected: Iterable[int]) -> float:
    return sum(weights[read_id] for read_id in selected)


def is_feasible(matrix: FragmentMatrix, selected: Sequence[int], coverage: int) -> bool:
    """Whether a selection respects the coverage cap at every variant."""
    return matrix.subset(selected).max_coverage() <= coverage


class _MinCostFlow:
    """Successive shortest paths with potentials, on a directed acyclic graph.

    Costs may be negative. The initial potentials come from a single pass in
    topological order, which is why every arc must run from a lower node index
    to a higher one.
    """

    def __init__(self, n_nodes: int) -> None:
        self.n_nodes = n_nodes
        self.to: list[int] = []
        self.cap: list[int] = []
        self.cost: list[float] = []
        self.out: list[list[int]] = [[] for _ in range(n_nodes)]

    def add_edge(self, tail: int, head: int, cap: int, cost: float) -> int:
        if tail >= head:
            raise ValueError("edges must run forward for the topological pass")
        index = len(self.to)
        self.to.append(head)
        self.cap.append(cap)
        self.cost.append(cost)
        self.out[tail].append(index)
        self.to.append(tail)
        self.cap.append(0)
        self.cost.append(-cost)
        self.out[head].append(index + 1)
        return index

    def _initial_potentials(self, source: int) -> list[float]:
        distance = [INF] * self.n_nodes
        distance[source] = 0.0
        for node in range(source, self.n_nodes):
            if distance[node] == INF:
                continue
            for edge in self.out[node]:
                if self.cap[edge] <= 0:
                    continue
                head = self.to[edge]
                relaxed = distance[node] + self.cost[edge]
                if relaxed < distance[head]:
                    distance[head] = relaxed
        return distance

    def _dijkstra(
        self, source: int, potential: list[float]
    ) -> tuple[list[float], list[int]]:
        distance = [INF] * self.n_nodes
        arrived = [-1] * self.n_nodes
        distance[source] = 0.0
        queue: list[tuple[float, int]] = [(0.0, source)]
        settled = [False] * self.n_nodes
        while queue:
            here, node = heapq.heappop(queue)
            if settled[node]:
                continue
            settled[node] = True
            for edge in self.out[node]:
                if self.cap[edge] <= 0:
                    continue
                head = self.to[edge]
                if potential[node] == INF or potential[head] == INF:
                    continue
                reduced = self.cost[edge] + potential[node] - potential[head]
                if reduced < 0.0:
                    # Only rounding can make a reduced cost negative here, and
                    # a negative step would break the Dijkstra invariant.
                    reduced = 0.0
                if here + reduced < distance[head]:
                    distance[head] = here + reduced
                    arrived[head] = edge
                    heapq.heappush(queue, (distance[head], head))
        return distance, arrived

    def run(self, source: int, sink: int, units: int) -> int:
        """Push up to `units` at minimum cost; return the units actually sent."""
        potential = self._initial_potentials(source)
        sent = 0
        while sent < units:
            distance, arrived = self._dijkstra(source, potential)
            if distance[sink] == INF:
                break
            for node in range(self.n_nodes):
                if distance[node] < INF and potential[node] < INF:
                    potential[node] += distance[node]
            bottleneck = units - sent
            node = sink
            while node != source:
                edge = arrived[node]
                bottleneck = min(bottleneck, self.cap[edge])
                node = self.to[edge ^ 1]
            node = sink
            while node != source:
                edge = arrived[node]
                self.cap[edge] -= bottleneck
                self.cap[edge ^ 1] += bottleneck
                node = self.to[edge ^ 1]
            sent += bottleneck
        return sent


def select_optimal(
    matrix: FragmentMatrix,
    coverage: int,
    weights: Optional[dict[int, float]] = None,
) -> SelectionResult:
    """Maximum-weight read subset under the coverage cap, solved exactly.

    Reads with non-positive weight never enter the flow network. Keeping them
    would leave the objective alone, since the slack can travel on the free
    path arcs instead, but a zero-weight read admitted on a tie would still
    consume coverage that a later stage could have used.
    """
    check_coverage(coverage)
    if weights is None:
        weights = read_weights(matrix)
    check_weights(weights, matrix.read_ids)
    if coverage == 0 or len(matrix) == 0 or matrix.n_variants == 0:
        return SelectionResult((), 0.0, "optimal", coverage)

    nodes = matrix.n_variants + 1
    flow = _MinCostFlow(nodes)
    for j in range(matrix.n_variants):
        flow.add_edge(j, j + 1, coverage, 0.0)
    arc_of: dict[int, int] = {}
    for fragment in matrix:
        weight = weights[fragment.read_id]
        if weight <= 0.0:
            continue
        arc_of[fragment.read_id] = flow.add_edge(
            fragment.start, fragment.end, 1, -weight
        )
    flow.run(0, matrix.n_variants, coverage)
    selected = tuple(
        read_id for read_id, arc in sorted(arc_of.items()) if flow.cap[arc] == 0
    )
    return SelectionResult(
        selected=selected,
        objective=selection_objective(weights, selected),
        method="optimal",
        coverage=coverage,
    )


def select_connected(
    matrix: FragmentMatrix,
    coverage: int,
    weights: Optional[dict[int, float]] = None,
    *,
    k: int = 1,
    penalty: float = 0.0,
) -> SelectionResult:
    """Maximum weight under the cap, minus `penalty` per missing spanning read.

    `docs/design/connectivity-selection.md` states the program and the flow
    network. Every boundary that some positively weighted read spans asks for
    `k` selected reads spanning it; each read short of that costs `penalty`.
    `penalty = 0` gives the additive optimum. The result's `objective` is the
    additive objective of the selection, without the penalty, so it compares
    with every other selector's.
    """
    check_coverage(coverage)
    if k < 0:
        raise ValueError("k must be non-negative")
    if penalty < 0:
        raise ValueError("penalty must be non-negative")
    k = min(k, coverage)
    costs = [penalty * max(k - s, 0) for s in range(coverage + 1)]
    return _select_with_boundary_costs(matrix, coverage, weights, costs, "connected")


def select_error_optimal(
    matrix: FragmentMatrix,
    coverage: int,
    weights: Optional[dict[int, float]] = None,
    *,
    costs: Sequence[float],
) -> SelectionResult:
    """Maximum weight under the cap, minus `costs[s]` at every boundary spanned `s` times.

    `docs/design/error-optimal-selection.md` states the program. `costs` has
    one entry for each count from 0 to `coverage`, and must be non-increasing
    and convex, which is what keeps the flow exact; `switch_model_costs` gives
    the model's choice. Only boundaries that some positively weighted read
    spans are charged. The result's `objective` is the additive objective.
    """
    check_coverage(coverage)
    if len(costs) != coverage + 1:
        raise ValueError("costs needs one entry for each count from 0 to coverage")
    steps = [costs[s - 1] - costs[s] for s in range(1, coverage + 1)]
    tolerance = 1e-12 * max(1.0, max(abs(x) for x in costs))
    if any(step < -tolerance for step in steps):
        raise ValueError("costs must be non-increasing")
    if any(later > earlier + tolerance for earlier, later in zip(steps, steps[1:])):
        raise ValueError("costs must be convex")
    return _select_with_boundary_costs(matrix, coverage, weights, list(costs),
                                       "error_optimal")


def switch_model_costs(coverage: int, allele_error: float, scale: float) -> list[float]:
    """`scale * rho**s` for `s` in `0..coverage`, the Bhattacharyya bound on a vote of `s` reads.

    A read calls a phase relation wrong when exactly one of its two alleles is
    wrong, with probability `q = 2e(1 - e)`, and `rho = 2 sqrt(q(1 - q))`.
    """
    if not 0.0 < allele_error < 0.5:
        raise ValueError("allele_error must lie strictly between 0 and 0.5")
    if scale < 0:
        raise ValueError("scale must be non-negative")
    q = 2 * allele_error * (1 - allele_error)
    rho = 2 * (q * (1 - q)) ** 0.5
    return [scale * rho ** s for s in range(coverage + 1)]


def _select_with_boundary_costs(
    matrix: FragmentMatrix,
    coverage: int,
    weights: Optional[dict[int, float]],
    costs: list[float],
    method: str,
) -> SelectionResult:
    """The split network of `connectivity-selection.md` with convex boundary costs.

    The flow on boundary `b`'s arcs is `coverage - s_b`, so its `u`-th unit
    costs `costs[coverage - u] - costs[coverage - u + 1]`; runs of equal unit
    costs share one arc.
    """
    if weights is None:
        weights = read_weights(matrix)
    check_weights(weights, matrix.read_ids)
    n = matrix.n_variants
    if coverage == 0 or len(matrix) == 0 or n == 0:
        return SelectionResult((), 0.0, method, coverage)
    reads = [f for f in matrix if weights[f.read_id] > 0.0]
    spanning = [0] * max(n - 1, 0)
    for fragment in reads:
        for b in range(fragment.start, fragment.end - 1):
            spanning[b] += 1
    units = [costs[coverage - u] - costs[coverage - u + 1] for u in range(1, coverage + 1)]
    arcs: list[list[float]] = []
    for unit in units:
        if arcs and unit == arcs[-1][1]:
            arcs[-1][0] += 1
        else:
            arcs.append([1, unit])
    # Node 2j is before variant j, 2j + 1 after it; boundary j joins 2j + 1 to 2j + 2.
    flow = _MinCostFlow(2 * n)
    for j in range(n):
        flow.add_edge(2 * j, 2 * j + 1, coverage, 0.0)
    for b in range(n - 1):
        if spanning[b]:
            for capacity, cost in arcs:
                flow.add_edge(2 * b + 1, 2 * b + 2, int(capacity), cost)
        else:
            flow.add_edge(2 * b + 1, 2 * b + 2, coverage, 0.0)
    arc_of = {
        fragment.read_id: flow.add_edge(2 * fragment.start, 2 * fragment.end - 1, 1,
                                        -weights[fragment.read_id])
        for fragment in reads
    }
    flow.run(0, 2 * n - 1, coverage)
    selected = tuple(read_id for read_id, arc in sorted(arc_of.items()) if flow.cap[arc] == 0)
    return SelectionResult(
        selected=selected,
        objective=selection_objective(weights, selected),
        method=method,
        coverage=coverage,
    )


def select_greedy(
    matrix: FragmentMatrix,
    coverage: int,
    weights: Optional[dict[int, float]] = None,
) -> SelectionResult:
    """Take reads in decreasing weight, skipping any that would exceed the cap.

    This is the greedy skeleton shared by deployed read selectors, including
    the one in WhatsHap. It is a baseline here, and `docs/literature/`
    records which parts of the published Fischer and Marschall procedure it
    does and does not reproduce. Ties are broken by read id so the result is
    reproducible.
    """
    check_coverage(coverage)
    if weights is None:
        weights = read_weights(matrix)
    check_weights(weights, matrix.read_ids)
    counts = [0] * matrix.n_variants
    selected: list[int] = []
    order = sorted(matrix, key=lambda f: (-weights[f.read_id], f.read_id))
    for fragment in order:
        if weights[fragment.read_id] <= 0.0:
            continue
        if any(counts[j] >= coverage for j in range(fragment.start, fragment.end)):
            continue
        for j in range(fragment.start, fragment.end):
            counts[j] += 1
        selected.append(fragment.read_id)
    selected.sort()
    return SelectionResult(
        selected=tuple(selected),
        objective=selection_objective(weights, selected),
        method="greedy",
        coverage=coverage,
    )


def select_random(
    matrix: FragmentMatrix,
    coverage: int,
    seed: int = 0,
    weights: Optional[dict[int, float]] = None,
) -> SelectionResult:
    """Downsample uniformly at random under the coverage cap.

    This is the baseline WhatsHap used before Fischer and Marschall 2016,
    described there as enumerating the reads in random order and retaining
    each one "if it is phase informative and adding it does not violate the
    coverage constraint". The phase-informative filter is part of the
    baseline, so it is applied here too. Without it the comparison would
    credit random selection for covering variants with reads that constrain
    nothing.
    """
    check_coverage(coverage)
    if weights is None:
        weights = read_weights(matrix)
    rng = random.Random(seed)
    order = [f for f in matrix if is_phase_informative(f)]
    rng.shuffle(order)
    counts = [0] * matrix.n_variants
    selected: list[int] = []
    for fragment in order:
        if any(counts[j] >= coverage for j in range(fragment.start, fragment.end)):
            continue
        for j in range(fragment.start, fragment.end):
            counts[j] += 1
        selected.append(fragment.read_id)
    selected.sort()
    return SelectionResult(
        selected=tuple(selected),
        objective=selection_objective(weights, selected),
        method="random",
        coverage=coverage,
    )


def select_bruteforce(
    matrix: FragmentMatrix,
    coverage: int,
    weights: Optional[dict[int, float]] = None,
    *,
    max_fragments: int = 18,
) -> SelectionResult:
    """Enumerate every subset. The oracle `select_optimal` is checked against."""
    check_coverage(coverage)
    if weights is None:
        weights = read_weights(matrix)
    check_weights(weights, matrix.read_ids)
    if len(matrix) > max_fragments:
        raise ValueError(
            f"{len(matrix)} fragments exceeds the brute-force limit of {max_fragments}"
        )
    ids = matrix.read_ids
    best_value = 0.0
    best: tuple[int, ...] = ()
    for code in range(1 << len(ids)):
        candidate = tuple(ids[i] for i in range(len(ids)) if code >> i & 1)
        counts = [0] * matrix.n_variants
        feasible = True
        for read_id in candidate:
            fragment = matrix[matrix.index_of(read_id)]
            for j in range(fragment.start, fragment.end):
                counts[j] += 1
                if counts[j] > coverage:
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
        method="bruteforce",
        coverage=coverage,
    )


class _UnionFind:
    """Disjoint sets over variant indices, for the bridging step."""

    def __init__(self, size: int) -> None:
        self._parent = list(range(size))

    def find(self, item: int) -> int:
        root = item
        while self._parent[root] != root:
            root = self._parent[root]
        while self._parent[item] != root:
            self._parent[item], item = root, self._parent[item]
        return root

    def union(self, left: int, right: int) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            self._parent[max(left_root, right_root)] = min(left_root, right_root)

    def components_of(self, items: Iterable[int]) -> set[int]:
        return {self.find(item) for item in items}


class _CoverageMonitor:
    """Running per-variant coverage with a cap, shared across slices."""

    def __init__(self, n_variants: int, cap: int) -> None:
        self._counts = [0] * n_variants
        self._cap = cap

    def admits(self, fragment: Fragment) -> bool:
        return all(
            self._counts[j] < self._cap for j in range(fragment.start, fragment.end)
        )

    def add(self, fragment: Fragment) -> None:
        for j in range(fragment.start, fragment.end):
            self._counts[j] += 1


def fischer_marschall_scores(
    fragment: Fragment, covered: set[int]
) -> tuple[int, int, float]:
    """The `(dynamic, static, quality)` tuple of Fischer and Marschall 2016.

    `scoredyn` counts the variants this read would newly cover in the current
    slice, `scorestatic` is the same count without the slice context, and both
    subtract the read's holes. `scorequal` is the worst base quality the read
    calls, which here is the smallest confidence weight. Reads are ordered by
    the tuple, larger first.
    """
    called = {position for position, _, _ in fragment.observed()}
    gaps = holes(fragment)
    dynamic = len(called) - len(called & covered) - gaps
    static = len(called) - gaps
    quality = min(
        (weight for _, _, weight in fragment.observed()), default=0.0
    )
    return dynamic, static, quality


def select_fischer_marschall(
    matrix: FragmentMatrix,
    coverage: int,
    weights: Optional[dict[int, float]] = None,
    *,
    bridging: bool = True,
) -> SelectionResult:
    """Read selection as published by Fischer and Marschall 2016.

    This follows the paper's pseudocode. WhatsHap's code departs from it and
    chooses different reads on most instances; `select_whatshap` is the port
    of the deployed code, and it is the one to use as the WhatsHap baseline.

    The procedure builds disjoint slices. Each slice greedily takes the
    best-scoring read that covers a variant the slice has not covered yet, and
    then, optionally, adds reads that join two connected components of the
    slice. Reads that would breach the coverage cap are discarded permanently.

    Two things in the published pseudocode are underdetermined and are resolved
    here, both in the only way that keeps the coverage constraint true:

    * the cap is checked against every read selected so far, not against the
      current slice alone;
    * bridging reads are removed from the undecided set, so a later slice
      cannot select them a second time.

    The objective is not the additive one that `select_optimal` maximizes.
    `scoredyn` pays for a variant the first time it is covered and not
    afterwards, which makes the slice objective a coverage function rather than
    a sum over reads. Comparing the two selectors on either objective alone
    therefore compares two different problems, which is what OQ-2 was about.
    """
    check_coverage(coverage)
    if weights is None:
        weights = read_weights(matrix)
    check_weights(weights, matrix.read_ids)

    monitor = _CoverageMonitor(matrix.n_variants, coverage)
    undecided = [f for f in matrix if is_phase_informative(f)]
    selected: list[int] = []

    while undecided:
        slice_reads, violating = _select_slice(undecided, monitor)
        if not slice_reads and not violating:
            break
        selected.extend(fragment.read_id for fragment in slice_reads)
        spent = {fragment.read_id for fragment in slice_reads}
        spent.update(fragment.read_id for fragment in violating)
        undecided = [f for f in undecided if f.read_id not in spent]
        if bridging:
            bridge = _bridge(slice_reads, undecided, monitor, matrix.n_variants)
            selected.extend(fragment.read_id for fragment in bridge)
            bridged = {fragment.read_id for fragment in bridge}
            undecided = [f for f in undecided if f.read_id not in bridged]

    selected.sort()
    return SelectionResult(
        selected=tuple(selected),
        objective=selection_objective(weights, selected),
        method="fischer_marschall" if bridging else "fischer_marschall_no_bridge",
        coverage=coverage,
    )


def _ordered_queue(fragments: Sequence[Fragment], covered: set[int]):
    """A max-heap of `(score, read_id)` with lazy rescoring."""
    heap = [
        (_negate(fischer_marschall_scores(fragment, covered)), fragment.read_id)
        for fragment in fragments
    ]
    heapq.heapify(heap)
    return heap


def _negate(score: tuple[int, int, float]) -> tuple[int, int, float]:
    return (-score[0], -score[1], -score[2])


def _select_slice(
    undecided: Sequence[Fragment], monitor: _CoverageMonitor
) -> tuple[list[Fragment], list[Fragment]]:
    """One slice: cover every variant once, best-scoring read first."""
    by_id = {fragment.read_id: fragment for fragment in undecided}
    covered: set[int] = set()
    chosen: list[Fragment] = []
    violating: list[Fragment] = []
    heap = _ordered_queue(undecided, covered)
    while heap:
        stored, read_id = heapq.heappop(heap)
        fragment = by_id[read_id]
        current = _negate(fischer_marschall_scores(fragment, covered))
        if current != stored:
            # A read admitted since this entry was pushed changed its dynamic
            # score. Re-push rather than act on the stale value.
            heapq.heappush(heap, (current, read_id))
            continue
        if not monitor.admits(fragment):
            violating.append(fragment)
            continue
        called = {position for position, _, _ in fragment.observed()}
        if not called - covered:
            continue
        monitor.add(fragment)
        chosen.append(fragment)
        covered |= called
    return chosen, violating


def _bridge(
    slice_reads: Sequence[Fragment],
    undecided: Sequence[Fragment],
    monitor: _CoverageMonitor,
    n_variants: int,
) -> list[Fragment]:
    """Add reads that join two connected components of the current slice."""
    components = _UnionFind(n_variants)
    for fragment in slice_reads:
        positions = [position for position, _, _ in fragment.observed()]
        for other in positions[1:]:
            components.union(positions[0], other)
    covered: set[int] = set()
    by_id = {fragment.read_id: fragment for fragment in undecided}
    heap = _ordered_queue(undecided, covered)
    added: list[Fragment] = []
    while heap:
        _, read_id = heapq.heappop(heap)
        fragment = by_id[read_id]
        positions = [position for position, _, _ in fragment.observed()]
        if len(components.components_of(positions)) < 2:
            continue
        if not monitor.admits(fragment):
            continue
        monitor.add(fragment)
        added.append(fragment)
        for other in positions[1:]:
            components.union(positions[0], other)
    return added


class _StdIntSet:
    """The iteration order of libstdc++'s `std::unordered_set<int>`.

    WhatsHap walks the variants a read newly covers in this order, which fixes
    the order of later score changes and so the heap's tie-breaking. The
    bucket counts are those of GCC 13's prime rehash policy, measured, and
    elements are kept in libstdc++'s node order: a new element goes to the
    front of its bucket's run, or to the front of the list if its bucket is
    empty. `clear` keeps the bucket count, as the C++ container does.
    """

    _BUCKETS = (13, 29, 59, 127, 257, 541, 1109, 2357, 5087, 10273, 20753, 42043)

    def __init__(self) -> None:
        self.buckets = 1
        self.order: list[int] = []

    def clear(self) -> None:
        self.order = []

    def insert(self, value: int) -> None:
        if value in self.order:
            return
        if len(self.order) + 1 > self.buckets:
            self.buckets = next(b for b in self._BUCKETS if b > self.buckets)
            old, self.order = self.order, []
            for item in old:
                self._place(item)
        self._place(value)

    def _place(self, value: int) -> None:
        bucket = value % self.buckets
        for index, item in enumerate(self.order):
            if item % self.buckets == bucket:
                self.order.insert(index, value)
                return
        self.order.insert(0, value)


class _WhatsHapHeap:
    """A port of `whatshap/priorityqueue.pyx`: a binary max-heap with an index.

    Ties between equal scores are broken by the heap's layout, which depends on
    the order of pushes and score changes. Reproducing WhatsHap's selection
    exactly therefore needs this heap, not `heapq` with a tie-breaking key.
    """

    def __init__(self) -> None:
        self.heap: list[list] = []
        self.positions: dict[int, int] = {}

    def __bool__(self) -> bool:
        return bool(self.heap)

    def push(self, score: tuple, item: int) -> None:
        self.heap.append([score, item])
        self.positions[item] = len(self.heap) - 1
        self._sift_up(len(self.heap) - 1)

    def pop(self) -> tuple[tuple, int]:
        first = self.heap[0]
        last = self.heap.pop()
        del self.positions[first[1]]
        if self.heap:
            self.heap[0] = last
            self.positions[last[1]] = 0
            self._sift_down(0)
        return first[0], first[1]

    def score(self, item: int) -> Optional[tuple]:
        index = self.positions.get(item)
        return None if index is None else self.heap[index][0]

    def change_score(self, item: int, score: tuple) -> None:
        index = self.positions[item]
        old = self.heap[index][0]
        self.heap[index][0] = score
        if old < score:
            self._sift_up(index)
        else:
            self._sift_down(index)

    def _lower(self, first: int, second: int) -> bool:
        return self.heap[first][0] < self.heap[second][0]

    def _swap(self, first: int, second: int) -> None:
        heap = self.heap
        heap[first], heap[second] = heap[second], heap[first]
        self.positions[heap[first][1]] = first
        self.positions[heap[second][1]] = second

    def _sift_up(self, index: int) -> None:
        while index > 0:
            parent = (index - 1) // 2
            if not self._lower(parent, index):
                return
            self._swap(parent, index)
            index = parent

    def _sift_down(self, index: int) -> None:
        size = len(self.heap)
        while True:
            left, right = 2 * index + 1, 2 * index + 2
            if right < size:
                child = right if self._lower(left, right) else left
            elif left < size:
                child = left
            else:
                return
            if not self._lower(index, child):
                return
            self._swap(child, index)
            index = child


def whatshap_score_drop(variants: Sequence[int], newly_covered: set[int]) -> int:
    """How much WhatsHap 2.8 lowers a queued read's first score after an admission.

    `_update_score_for_reads` in `readselect.pyx` subtracts one for every
    variant of the queued read that is absent from the set the admitted read
    newly covered. The docstring of `_compute_score_for_read` describes the
    score as new variants minus gaps, which would subtract the variants that
    are present instead. This function reproduces the code, not the docstring.
    """
    return sum(1 for variant in variants if variant not in newly_covered)


def select_whatshap(
    matrix: FragmentMatrix,
    coverage: int,
    weights: Optional[dict[int, float]] = None,
    *,
    bridging: bool = True,
    positions: Optional[Sequence[int]] = None,
) -> SelectionResult:
    """Read selection as WhatsHap 2.8 implements it in `readselect.pyx`.

    This is a line-by-line port of `readselection` with one read source, run
    on the phase-informative reads in matrix order (WhatsHap's read index is
    the position in its `ReadSet`). It departs from the published pseudocode
    that `select_fischer_marschall` follows in four places, and
    `tests/test_whatshap_port.py` checks the port against WhatsHap itself:

    * When a read is admitted, every overlapping read still in the queue has
      its first score lowered by the number of its variants that the admitted
      read did not newly cover. The paper's `scoredyn` would lower it by the
      number that were newly covered.
    * Scores change in place inside the heap and ties follow the heap layout.
    * During bridging the coverage cap is checked before the component test,
      and a read that would breach it is discarded for good.
    * Holes and coverage are counted over the variant positions that some
      phase-informative read calls, not over all `n_variants`.

    The quality score is the smallest weight a read calls, as in
    `select_fischer_marschall`; WhatsHap's integer base qualities give the
    same order when the weights are integers. `positions` gives the genomic
    position of each variant index. It changes nothing but the order in which
    tied reads leave the heap, because WhatsHap hashes positions; without it
    the variant index is used.
    """
    check_coverage(coverage)
    if weights is None:
        weights = read_weights(matrix)
    check_weights(weights, matrix.read_ids)

    if positions is None:
        positions = range(matrix.n_variants)
    genomic = list(positions)
    if len(genomic) != matrix.n_variants:
        raise ValueError("positions must give one position per variant")
    reads = [f for f in matrix if is_phase_informative(f)]
    called = [[genomic[j] for j, _, _ in f.observed()] for f in reads]
    quality = [min(weight for _, _, weight in f.observed()) for f in reads]
    present = sorted({position for row in called for position in row})
    index_of = {position: index for index, position in enumerate(present)}
    variant_to_reads: dict[int, list[int]] = {}
    for read, row in enumerate(called):
        for position in row:
            variant_to_reads.setdefault(index_of[position], []).append(read)
    counts = [0] * len(present)

    def score_of(read: int) -> tuple:
        row = called[read]
        good = len(row)
        bad = index_of[row[-1]] - index_of[row[0]] + 1 - good
        return (good - bad, good - bad, quality[read])

    def queue_of(indices: Iterable[int]) -> _WhatsHapHeap:
        queue = _WhatsHapHeap()
        for read in indices:
            queue.push(score_of(read), read)
        return queue

    def interval(read: int) -> tuple[int, int]:
        return index_of[called[read][0]], index_of[called[read][-1]] + 1

    def full(begin: int, end: int) -> bool:
        return max(counts[begin:end]) >= coverage

    def admit(begin: int, end: int) -> None:
        for j in range(begin, end):
            counts[j] += 1

    selected: set[int] = set()
    undecided = set(range(len(reads)))
    while undecided:
        queue = queue_of(undecided)
        covered: set[int] = set()
        in_slice: set[int] = set()
        violating: set[int] = set()
        new = _StdIntSet()
        while queue:
            _, read = queue.pop()
            new.clear()
            for position in called[read]:
                if position not in covered:
                    new.insert(position)
            begin, end = interval(read)
            if full(begin, end):
                violating.add(read)
            elif new.order:
                admit(begin, end)
                in_slice.add(read)
                affected: set[int] = set()
                for position in new.order:
                    covered.add(position)
                    affected.update(variant_to_reads[index_of[position]])
                fresh = set(new.order)
                for other in affected.difference(in_slice):
                    old = queue.score(other)
                    if old is not None:
                        drop = whatshap_score_drop(called[other], fresh)
                        queue.change_score(other, (old[0] - drop, old[1], old[2]))
        selected.update(in_slice)
        undecided -= in_slice
        undecided -= violating

        components = _UnionFind(len(present))
        for read in in_slice:
            for position in called[read][1:]:
                components.union(index_of[called[read][0]], index_of[position])
        if bridging:
            queue = queue_of(undecided)
            while queue:
                _, read = queue.pop()
                blocks = components.components_of(index_of[p] for p in called[read])
                begin, end = interval(read)
                if full(begin, end):
                    undecided.remove(read)
                    continue
                if len(blocks) < 2:
                    continue
                selected.add(read)
                admit(begin, end)
                undecided.remove(read)
                for position in called[read][1:]:
                    components.union(index_of[called[read][0]], index_of[position])

    chosen = sorted(reads[read].read_id for read in selected)
    return SelectionResult(
        selected=tuple(chosen),
        objective=selection_objective(weights, chosen),
        method="whatshap" if bridging else "whatshap_no_bridge",
        coverage=coverage,
    )


def min_coverage(matrix: FragmentMatrix, selected: Iterable[int]) -> int:
    """Smallest per-variant coverage of a selection, over every variant."""
    subset = matrix.subset(selected)
    return min(subset.coverage(), default=0)


def distinct_variants_covered(matrix: FragmentMatrix, selected: Iterable[int]) -> int:
    """Variants with at least one called allele among the selected reads."""
    subset = matrix.subset(selected)
    positions: set[int] = set()
    for fragment in subset:
        positions.update(position for position, _, _ in fragment.observed())
    return len(positions)


def component_count(matrix: FragmentMatrix, selected: Iterable[int]) -> int:
    """Number of phase blocks a selection leaves, counting covered variants."""
    blocks = matrix.subset(selected).blocks()
    return len({block for block in blocks if block is not None})


def _coverage_feasible(
    matrix: FragmentMatrix,
    coverage: int,
    floor: int,
    weights: Optional[dict[int, float]] = None,
    *,
    informative_only: bool = True,
) -> Optional[tuple[int, ...]]:
    """A selection with `maxcov <= coverage` and `mincov >= floor`, or `None`.

    The network is the one in Mäkinen, Staneva, Tomescu and Valenzuela 2015:
    backbone arcs of capacity `coverage - floor` between consecutive variant
    boundaries, capacity `coverage` at the two ends, unit arcs for the reads,
    and `coverage` units pushed end to end. Their Theorem 1 says the instance
    is feasible exactly when that flow saturates.

    Passing `weights` costs the read arcs by `-w_r`, which picks the
    maximum-weight selection among the feasible ones. That is a min-cost flow
    on the same network and needs no separate argument.
    """
    m = matrix.n_variants
    source, sink = 0, m + 2
    flow = _MinCostFlow(m + 3)
    flow.add_edge(source, 1, coverage, 0.0)
    backbone = max(coverage - floor, 0)
    for j in range(m):
        flow.add_edge(j + 1, j + 2, backbone, 0.0)
    flow.add_edge(m + 1, sink, coverage, 0.0)
    arc_of: dict[int, int] = {}
    for fragment in matrix:
        if informative_only and not is_phase_informative(fragment):
            continue
        cost = 0.0 if weights is None else -weights[fragment.read_id]
        arc_of[fragment.read_id] = flow.add_edge(
            fragment.start + 1, fragment.end + 1, 1, cost
        )
    sent = flow.run(source, sink, coverage)
    if sent < coverage:
        return None
    return tuple(
        read_id for read_id, arc in sorted(arc_of.items()) if flow.cap[arc] == 0
    )


def select_max_min_coverage(
    matrix: FragmentMatrix,
    coverage: int,
    weights: Optional[dict[int, float]] = None,
    *,
    break_ties_by_weight: bool = True,
    informative_only: bool = True,
) -> SelectionResult:
    """Maximize the smallest per-variant coverage, subject to the cap.

    This is the objective of Mäkinen, Staneva, Tomescu and Valenzuela 2015,
    solved exactly by binary search over the coverage floor with one max-flow
    per step. Feasibility is monotone in the floor, so the search is valid: a
    selection with `mincov >= t` also has `mincov >= t-1`.

    With `break_ties_by_weight`, the selection returned at the optimal floor is
    the maximum-weight one among those achieving it. The floor is decided
    first and the weight only breaks ties, so the reported floor is the
    optimum either way.

    A variant that no read covers forces the optimum to zero, and the answer is
    then any feasible selection. `min_coverage` reports the floor so a caller
    can see when that has happened.

    `informative_only` is a policy choice, not part of the published problem.
    The paper selects over whatever intervals it is given, and a read calling
    one allele raises the coverage floor while constraining no phase relation.
    Set it to `False` to solve the problem exactly as stated there.
    """
    check_coverage(coverage)
    if weights is None:
        weights = read_weights(matrix)
    check_weights(weights, matrix.read_ids)
    if coverage == 0 or matrix.n_variants == 0:
        return SelectionResult((), 0.0, "max_min_coverage", coverage)

    best = _coverage_feasible(
        matrix, coverage, 0, None, informative_only=informative_only
    )
    if best is None:
        return SelectionResult((), 0.0, "max_min_coverage", coverage)
    low, high = 0, coverage
    while low < high:
        middle = (low + high + 1) // 2
        candidate = _coverage_feasible(
            matrix, coverage, middle, None, informative_only=informative_only
        )
        if candidate is None:
            high = middle - 1
        else:
            low, best = middle, candidate
    if break_ties_by_weight:
        refined = _coverage_feasible(
            matrix, coverage, low, weights, informative_only=informative_only
        )
        if refined is not None:
            best = refined
    return SelectionResult(
        selected=best,
        objective=selection_objective(weights, best),
        method="max_min_coverage",
        coverage=coverage,
    )


def select_max_min_coverage_bruteforce(
    matrix: FragmentMatrix,
    coverage: int,
    *,
    max_fragments: int = 16,
    informative_only: bool = True,
) -> tuple[int, tuple[int, ...]]:
    """Enumerate subsets and return `(best floor, a selection achieving it)`."""
    check_coverage(coverage)
    if len(matrix) > max_fragments:
        raise ValueError(f"{len(matrix)} fragments exceeds the brute-force limit")
    ids = tuple(
        fragment.read_id
        for fragment in matrix
        if is_phase_informative(fragment) or not informative_only
    )
    best_floor = -1
    best: tuple[int, ...] = ()
    for code in range(1 << len(ids)):
        candidate = tuple(ids[i] for i in range(len(ids)) if code >> i & 1)
        subset = matrix.subset(candidate)
        counts = subset.coverage()
        if max(counts, default=0) > coverage:
            continue
        floor = min(counts) if counts else 0
        if floor > best_floor:
            best_floor, best = floor, candidate
    return best_floor, best


SELECTORS: dict[str, Callable[..., SelectionResult]] = {
    "optimal": select_optimal,
    "fischer_marschall": select_fischer_marschall,
    "whatshap": select_whatshap,
    "connected": select_connected,
    "max_min_coverage": select_max_min_coverage,
    "greedy": select_greedy,
    "random": select_random,
    "bruteforce": select_bruteforce,
}
