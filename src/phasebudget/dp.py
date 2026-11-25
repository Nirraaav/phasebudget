"""Coverage-parameterized dynamic program for wMEC, and its phase margins.

The state at variant `j` is a labelled bipartition of the reads active at `j`
together with the allele of the first haplotype at `j`. The number of states is
`2^(coverage+1)`, so the running time is linear in the number of variants and
exponential in the active coverage. That is what makes coverage-constrained
read selection an optimization problem rather than a formality.

The margin of the boundary between variants `j` and `j+1` is the increase in
optimal cost forced by reversing `h[j] ^ h[j+1]`. Both values are read off the
forward and backward tables at that boundary, so all `n-1` margins cost one
forward pass and one backward pass rather than `n-1` re-solves.

Two backends compute the same thing. The pure-Python pass is the reference.
The native pass in `cpp/dp.cpp`, reached through `phasebudget._native`, repeats
its arithmetic and tie-breaking exactly and keeps memory at `O(sqrt(n) 2^k)`;
`backend="auto"` uses it when it is built.
"""

from __future__ import annotations

import math
from typing import Optional, Sequence

from . import _native
from ._types import Margin, WMECSolution
from ._validation import CoverageLimitExceeded
from .fragments import FragmentMatrix

#: Refuse rather than allocate 2^25 states per column.
MAX_ACTIVE_READS = 22

INF = float("inf")


def column_sum(values: Sequence[float]) -> float:
    """`sum(values)` as CPython 3.12 computes it for floats, in any version.

    CPython 3.12 changed the built-in `sum` to Neumaier's compensated
    summation. Spelling it out keeps the column totals, and so every cost in
    the dynamic program, identical across Python versions and identical to
    the native pass, which repeats this loop.
    """
    total = 0.0
    compensation = 0.0
    for value in values:
        value = float(value)
        partial = total + value
        if abs(total) >= abs(value):
            compensation += (total - partial) + value
        else:
            compensation += (value - partial) + total
        total = partial
    if compensation and math.isfinite(compensation):
        total += compensation
    return total


class CoverageDP:
    """Exact wMEC and exact phase margins under the complementary model.

    `solve` runs the forward pass and a traceback; `margins` additionally runs
    the backward pass. Both are cached. `backend` is "python", "native", or
    "auto", which is native when the library is built.
    """

    def __init__(
        self,
        matrix: FragmentMatrix,
        *,
        max_active_reads: int = MAX_ACTIVE_READS,
        backend: str = "auto",
    ) -> None:
        if backend not in ("auto", "python", "native"):
            raise ValueError(f"unknown backend {backend!r}")
        if backend == "native" and not _native.available():
            raise RuntimeError("the native dynamic program is not built; run `make -C cpp`")
        self.matrix = matrix
        self.n_variants = matrix.n_variants
        self._active = [matrix.active(j) for j in range(self.n_variants)]
        widest = max((len(a) for a in self._active), default=0)
        if widest > max_active_reads:
            raise CoverageLimitExceeded(
                f"active coverage {widest} exceeds the state limit of "
                f"{max_active_reads}; select reads first"
            )
        self.backend = (
            "native" if backend == "native" or (backend == "auto" and _native.available())
            else "python"
        )
        self._delta_tables: Optional[list[tuple[float, float, list[float]]]] = None
        self._boundaries: Optional[tuple[list[float], list[float]]] = None
        self._forward: Optional[list[list[float]]] = None
        self._parent: Optional[list[list[int]]] = None
        self._backward: Optional[list[list[float]]] = None
        self._solution: Optional[WMECSolution] = None

    # ---- column costs -------------------------------------------------

    @property
    def _delta(self) -> list[tuple[float, float, list[float]]]:
        """Per-column cost tables for the Python pass, built on first use."""
        if self._delta_tables is None:
            self._delta_tables = [self._column_tables(j) for j in range(self.n_variants)]
        return self._delta_tables

    def _column_weights(self, position: int) -> tuple[list[float], list[float]]:
        """`c0[i]` and `c1[i]` for the read in bit position `i` at `position`."""
        c0: list[float] = []
        c1: list[float] = []
        for index in self._active[position]:
            fragment = self.matrix[index]
            allele = fragment.allele_at(position)
            weight = fragment.weight_at(position)
            if allele is None:
                c0.append(0.0)
                c1.append(0.0)
            elif allele == 0:
                c0.append(0.0)
                c1.append(weight)
            else:
                c0.append(weight)
                c1.append(0.0)
        return c0, c1

    def _column_tables(self, position: int) -> tuple[float, float, list[float]]:
        """`(S0, S1, D)` with `cost(mask, a=0) = S0 + D[mask]`.

        `c0[i]` is the weighted cost of giving read `i` the value 0 at this
        variant and `c1[i]` the cost of giving it 1. A read in part 0 takes the
        value `a` and a read in part 1 takes `1-a`, so `cost(mask, 0)` is
        `sum(c0) + sum_(i in mask) (c1[i] - c0[i])` and, by the same argument
        with the roles exchanged, `cost(mask, 1) = sum(c1) - sum_(i in mask)
        (c1[i] - c0[i])`.
        """
        active = self._active[position]
        c0, c1 = self._column_weights(position)
        size = 1 << len(active)
        delta = [0.0] * size
        for i in range(len(active)):
            step = c1[i] - c0[i]
            bit = 1 << i
            for mask in range(bit, size):
                if mask & bit:
                    delta[mask] = delta[mask ^ bit] + step
        return column_sum(c0), column_sum(c1), delta

    def column_cost(self, position: int, mask: int, allele: int) -> float:
        """Weighted mismatches at one variant for one state."""
        s0, s1, delta = self._delta[position]
        return s0 + delta[mask] if allele == 0 else s1 - delta[mask]

    def _state_costs(self, position: int) -> list[float]:
        """Flat table indexed by `mask * 2 + allele`."""
        s0, s1, delta = self._delta[position]
        table = [0.0] * (len(delta) * 2)
        for mask, shift in enumerate(delta):
            table[mask * 2] = s0 + shift
            table[mask * 2 + 1] = s1 - shift
        return table

    # ---- projections between adjacent columns -------------------------

    def _shared(self, left: int, right: int) -> tuple[tuple[int, int], ...]:
        """Bit positions of reads active at both columns, as `(left, right)`."""
        right_at = {read: i for i, read in enumerate(self._active[right])}
        return tuple(
            (i, right_at[read])
            for i, read in enumerate(self._active[left])
            if read in right_at
        )

    @staticmethod
    def _project(mask: int, bits: Sequence[int]) -> int:
        code = 0
        for slot, bit in enumerate(bits):
            if mask >> bit & 1:
                code |= 1 << slot
        return code

    # ---- forward and backward passes ----------------------------------

    def _run_forward(self) -> None:
        if self._forward is not None:
            return
        forward: list[list[float]] = []
        parent: list[list[int]] = []
        if self.n_variants == 0:
            self._forward, self._parent = forward, parent
            return
        forward.append(self._state_costs(0))
        parent.append([-1] * len(forward[0]))
        for j in range(1, self.n_variants):
            shared = self._shared(j - 1, j)
            left_bits = [pair[0] for pair in shared]
            right_bits = [pair[1] for pair in shared]
            groups: dict[int, list[float]] = {}
            winners: dict[int, list[int]] = {}
            previous = forward[j - 1]
            for state, cost in enumerate(previous):
                if cost == INF:
                    continue
                code = self._project(state >> 1, left_bits)
                allele = state & 1
                best = groups.get(code)
                if best is None:
                    best = [INF, INF]
                    groups[code] = best
                    winners[code] = [-1, -1]
                if cost < best[allele]:
                    best[allele] = cost
                    winners[code][allele] = state
            costs = self._state_costs(j)
            table = [INF] * len(costs)
            back = [-1] * len(costs)
            for state, own in enumerate(costs):
                code = self._project(state >> 1, right_bits)
                best = groups.get(code)
                if best is None:
                    continue
                allele = 0 if best[0] <= best[1] else 1
                if best[allele] == INF:
                    continue
                table[state] = own + best[allele]
                back[state] = winners[code][allele]
            forward.append(table)
            parent.append(back)
        self._forward, self._parent = forward, parent

    def _run_backward(self) -> None:
        if self._backward is not None:
            return
        if self.n_variants == 0:
            self._backward = []
            return
        last = self.n_variants - 1
        backward: list[Optional[list[float]]] = [None] * self.n_variants
        backward[last] = self._state_costs(last)
        for j in range(last - 1, -1, -1):
            shared = self._shared(j, j + 1)
            left_bits = [pair[0] for pair in shared]
            right_bits = [pair[1] for pair in shared]
            groups: dict[int, float] = {}
            following = backward[j + 1]
            assert following is not None
            for state, cost in enumerate(following):
                if cost == INF:
                    continue
                code = self._project(state >> 1, right_bits)
                if cost < groups.get(code, INF):
                    groups[code] = cost
            costs = self._state_costs(j)
            table = [INF] * len(costs)
            for state, own in enumerate(costs):
                code = self._project(state >> 1, left_bits)
                tail = groups.get(code, INF)
                if tail < INF:
                    table[state] = own + tail
            backward[j] = table
        self._backward = [table for table in backward if table is not None]

    # ---- the native pass ---------------------------------------------

    def native_inputs(self) -> tuple[list[int], list[int], list[float], list[float]]:
        """`(offsets, active, c0, c1)`: the flat column layout the native pass reads."""
        offsets = [0]
        active: list[int] = []
        c0: list[float] = []
        c1: list[float] = []
        for j in range(self.n_variants):
            weights0, weights1 = self._column_weights(j)
            active.extend(self._active[j])
            c0.extend(weights0)
            c1.extend(weights1)
            offsets.append(len(active))
        return offsets, active, c0, c1

    def _run_native(self, want_margins: bool) -> None:
        offsets, active, c0, c1 = self.native_inputs()
        cost, haplotype, assignment, zero, one = _native.run(
            self.n_variants,
            offsets,
            active,
            c0,
            c1,
            len(self.matrix),
            want_margins=want_margins,
        )
        self._solution = WMECSolution(cost=cost, haplotype=haplotype, assignment=assignment)
        if zero is not None and one is not None:
            self._boundaries = (zero, one)

    # ---- results ------------------------------------------------------

    def solve(self) -> WMECSolution:
        """The optimal wMEC solution, with a traceback for the assignment."""
        if self._solution is not None:
            return self._solution
        if self.n_variants == 0:
            self._solution = WMECSolution(0.0, (), tuple(0 for _ in self.matrix))
            return self._solution
        if self.backend == "native":
            self._run_native(want_margins=False)
            assert self._solution is not None
            return self._solution
        self._run_forward()
        assert self._forward is not None and self._parent is not None
        last = self.n_variants - 1
        final = self._forward[last]
        best_state = min(range(len(final)), key=lambda s: final[s])
        if final[best_state] == INF:
            raise RuntimeError("no feasible state; the fragment matrix is malformed")
        states = [0] * self.n_variants
        state = best_state
        for j in range(last, -1, -1):
            states[j] = state
            if j:
                state = self._parent[j][state]
                if state < 0:
                    raise RuntimeError(f"traceback broke at variant {j}")
        haplotype = tuple(states[j] & 1 for j in range(self.n_variants))
        assignment = [0] * len(self.matrix)
        for j in range(self.n_variants):
            mask = states[j] >> 1
            for slot, index in enumerate(self._active[j]):
                assignment[index] = (mask >> slot) & 1
        self._solution = WMECSolution(
            cost=final[best_state],
            haplotype=haplotype,
            assignment=tuple(assignment),
        )
        return self._solution

    def boundary_costs(self, boundary: int) -> tuple[float, float]:
        """`(best with relation 0, best with relation 1)` at one boundary.

        The forward table covers variants `0..boundary` and the backward table
        covers `boundary+1..n-1`, and the transition between them adds no cost
        of its own, so the two tables add directly once the shared reads agree.
        """
        if not 0 <= boundary < self.n_variants - 1:
            raise ValueError(f"boundary {boundary} is outside the variant range")
        if self.backend == "native":
            if self._boundaries is None:
                self._run_native(want_margins=True)
            assert self._boundaries is not None
            return self._boundaries[0][boundary], self._boundaries[1][boundary]
        self._run_forward()
        self._run_backward()
        assert self._forward is not None and self._backward is not None
        shared = self._shared(boundary, boundary + 1)
        left_bits = [pair[0] for pair in shared]
        right_bits = [pair[1] for pair in shared]
        left: dict[int, list[float]] = {}
        for state, cost in enumerate(self._forward[boundary]):
            if cost == INF:
                continue
            code = self._project(state >> 1, left_bits)
            slot = left.setdefault(code, [INF, INF])
            allele = state & 1
            if cost < slot[allele]:
                slot[allele] = cost
        right: dict[int, list[float]] = {}
        for state, cost in enumerate(self._backward[boundary + 1]):
            if cost == INF:
                continue
            code = self._project(state >> 1, right_bits)
            slot = right.setdefault(code, [INF, INF])
            allele = state & 1
            if cost < slot[allele]:
                slot[allele] = cost
        best = [INF, INF]
        for code, head in left.items():
            tail = right.get(code)
            if tail is None:
                continue
            for allele in (0, 1):
                if head[allele] == INF:
                    continue
                for relation in (0, 1):
                    other = tail[allele ^ relation]
                    if other == INF:
                        continue
                    total = head[allele] + other
                    if total < best[relation]:
                        best[relation] = total
        return best[0], best[1]

    def margins(self) -> tuple[Margin, ...]:
        """The MEC margin of every adjacent phase decision.

        A boundary no read spans has equal cost either way, because the two
        sides can be relabelled independently, so its margin is zero. That is
        the block boundary condition, arrived at by the arithmetic rather than
        by a separate rule.
        """
        spanned = self.matrix.spanned_boundaries()
        result = []
        for boundary in range(self.n_variants - 1):
            zero, one = self.boundary_costs(boundary)
            relation = 0 if zero <= one else 1
            gap = abs(zero - one)
            result.append(
                Margin(
                    boundary=boundary,
                    relation=relation,
                    margin=gap,
                    spanned=spanned[boundary],
                )
            )
        return tuple(result)


def solve(matrix: FragmentMatrix, **kwargs) -> WMECSolution:
    """Convenience wrapper around `CoverageDP.solve`."""
    return CoverageDP(matrix, **kwargs).solve()


def margins(matrix: FragmentMatrix, **kwargs) -> tuple[Margin, ...]:
    """Convenience wrapper around `CoverageDP.margins`."""
    return CoverageDP(matrix, **kwargs).margins()
