"""Exact MEC under the general model, where the two haplotypes are free.

`phasebudget.dp` fixes the second haplotype to the complement of the first,
which cannot represent a site that was called heterozygous and is not. This is
the same coverage-parameterized recurrence without that restriction: each side
of the read bipartition is corrected to its own weighted consensus at every
variant. `docs/design/genotype-errors.md` states the recurrence, and
`exact.solve_bipartition_bruteforce` is the oracle it answers to.
"""

from __future__ import annotations

from typing import Sequence

from ._validation import CoverageLimitExceeded
from .fragments import FragmentMatrix

#: Refuse rather than allocate 2^25 states per column.
MAX_ACTIVE_READS = 22

INF = float("inf")


def _project(mask: int, bits: Sequence[int]) -> int:
    code = 0
    for slot, bit in enumerate(bits):
        if mask >> bit & 1:
            code |= 1 << slot
    return code


class GeneralDP:
    """Optimal general-model MEC and one optimal read bipartition."""

    def __init__(
        self, matrix: FragmentMatrix, *, max_active_reads: int = MAX_ACTIVE_READS
    ) -> None:
        self.matrix = matrix
        self.n_variants = matrix.n_variants
        self._active = [matrix.active(j) for j in range(self.n_variants)]
        widest = max((len(a) for a in self._active), default=0)
        if widest > max_active_reads:
            raise CoverageLimitExceeded(
                f"active coverage {widest} exceeds the state limit of "
                f"{max_active_reads}; select reads first"
            )

    def column_costs(self, position: int) -> list[float]:
        """`cost_j(S)` for every mask `S` of the reads active at `position`."""
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
        size = 1 << len(c0)
        a0 = [0.0] * size
        a1 = [0.0] * size
        for i in range(len(c0)):
            bit = 1 << i
            for mask in range(bit, size):
                if mask & bit:
                    a0[mask] = a0[mask ^ bit] + c0[i]
                    a1[mask] = a1[mask ^ bit] + c1[i]
        total0, total1 = a0[size - 1], a1[size - 1]
        return [
            min(total0 - a0[mask], total1 - a1[mask]) + min(a0[mask], a1[mask])
            for mask in range(size)
        ]

    def _shared(self, left: int, right: int) -> tuple[list[int], list[int]]:
        right_at = {read: i for i, read in enumerate(self._active[right])}
        pairs = [
            (i, right_at[read])
            for i, read in enumerate(self._active[left])
            if read in right_at
        ]
        return [p[0] for p in pairs], [p[1] for p in pairs]

    def solve(self) -> tuple[float, tuple[int, ...]]:
        """`(cost, assignment)`, the shape `solve_bipartition_bruteforce` returns."""
        n_reads = len(self.matrix)
        if self.n_variants == 0:
            return 0.0, tuple(0 for _ in range(n_reads))
        forward = [self.column_costs(0)]
        parents: list[list[int]] = [[-1] * len(forward[0])]
        for j in range(1, self.n_variants):
            left_bits, right_bits = self._shared(j - 1, j)
            best: dict[int, tuple[float, int]] = {}
            for mask, cost in enumerate(forward[j - 1]):
                code = _project(mask, left_bits)
                if code not in best or cost < best[code][0]:
                    best[code] = (cost, mask)
            costs = self.column_costs(j)
            table = [INF] * len(costs)
            back = [-1] * len(costs)
            for mask, own in enumerate(costs):
                entry = best.get(_project(mask, right_bits))
                if entry is None:
                    continue
                table[mask] = own + entry[0]
                back[mask] = entry[1]
            forward.append(table)
            parents.append(back)
        last = forward[-1]
        state = min(range(len(last)), key=lambda mask: last[mask])
        cost = last[state]
        assignment = [0] * n_reads
        for j in range(self.n_variants - 1, -1, -1):
            for slot, index in enumerate(self._active[j]):
                assignment[index] = (state >> slot) & 1
            if j:
                state = parents[j][state]
        return cost, tuple(assignment)


def assignment_cost(matrix: FragmentMatrix, assignment: Sequence[int]) -> float:
    """General-model cost of a fixed bipartition: each side takes its consensus."""
    cost = 0.0
    for position in range(matrix.n_variants):
        for side in (0, 1):
            to_zero = to_one = 0.0
            for index in matrix.active(position):
                if assignment[index] != side:
                    continue
                allele = matrix[index].allele_at(position)
                if allele is None:
                    continue
                if allele == 0:
                    to_one += matrix[index].weight_at(position)
                else:
                    to_zero += matrix[index].weight_at(position)
            cost += min(to_zero, to_one)
    return cost


def solve_general(matrix: FragmentMatrix, **kwargs) -> tuple[float, tuple[int, ...]]:
    """Convenience wrapper around `GeneralDP.solve`."""
    return GeneralDP(matrix, **kwargs).solve()
