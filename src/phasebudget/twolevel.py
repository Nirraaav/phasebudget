"""Two-level exact phasing: orient the subset's blocks using every read.

`docs/design/two-level-phasing.md` states the objective: keep the exact
haplotype of the selected subset inside each of its phase blocks, and choose
one orientation per block to minimize the full-matrix wMEC. Reads that call
alleles in several blocks are factors over a window of consecutive blocks,
so a dynamic program along the chain of blocks, whose state is the
orientation of the last blocks of the window, finds the exact minimum.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

from ._types import Phasing, WMECSolution
from .fragments import FragmentMatrix

INF = float("inf")


@dataclass(frozen=True)
class TwoLevelResult:
    """The oriented haplotype, its blocks, and what the orientation cost."""

    haplotype: tuple[int, ...]
    phasing: Phasing
    orientations: tuple[int, ...]
    cost: float
    subset_blocks: int
    window: int
    dropped_reads: int


@dataclass(frozen=True)
class _Factor:
    blocks: tuple[int, ...]
    agree_if_kept: tuple[float, ...]     # cost of the read's alleles when the block is not flipped
    agree_if_flipped: tuple[float, ...]  # cost when it is flipped

    def cost(self, orientation_of) -> float:
        side0 = side1 = 0.0
        for block, a, b in zip(self.blocks, self.agree_if_kept, self.agree_if_flipped):
            if orientation_of(block):
                side0 += b
                side1 += a
            else:
                side0 += a
                side1 += b
        return min(side0, side1)


def _block_order(subset: FragmentMatrix) -> tuple[list[Optional[int]], int]:
    """Per variant, the index of its subset block in chromosome order, or None."""
    raw = subset.blocks()
    order: dict[int, int] = {}
    indexed: list[Optional[int]] = []
    for block in raw:
        if block is None:
            indexed.append(None)
            continue
        if block not in order:
            order[block] = len(order)
        indexed.append(order[block])
    return indexed, len(order)


def _factors(matrix: FragmentMatrix, haplotype: Sequence[int], block_of, max_window: int):
    factors: list[_Factor] = []
    dropped = 0
    for fragment in matrix:
        kept: dict[int, float] = {}
        flipped: dict[int, float] = {}
        for position, allele, weight in fragment.observed():
            block = block_of[position]
            if block is None:
                continue
            if allele != haplotype[position]:
                kept[block] = kept.get(block, 0.0) + weight
                flipped.setdefault(block, 0.0)
            else:
                flipped[block] = flipped.get(block, 0.0) + weight
                kept.setdefault(block, 0.0)
        if len(kept) < 2:
            continue
        blocks = tuple(sorted(kept))
        if blocks[-1] - blocks[0] + 1 > max_window:
            dropped += 1
            continue
        factors.append(_Factor(blocks, tuple(kept[b] for b in blocks),
                               tuple(flipped[b] for b in blocks)))
    return factors, dropped


def _solve_chain(n_blocks: int, factors: Sequence[_Factor]) -> tuple[list[int], int]:
    """Exact orientations minimizing the summed factor costs, with block 0 fixed at 0."""
    if n_blocks == 0:
        return [], 0
    window = max((f.blocks[-1] - f.blocks[0] + 1 for f in factors), default=1)
    mask = (1 << window) - 1
    ending: list[list[_Factor]] = [[] for _ in range(n_blocks)]
    for factor in factors:
        ending[factor.blocks[-1]].append(factor)
    # State bit t holds the orientation of block k - t, bit 0 the current block.
    costs = {0: 0.0}
    back: list[dict[int, int]] = [{}]
    for k in range(1, n_blocks):
        new_costs: dict[int, float] = {}
        new_back: dict[int, int] = {}
        for state, value in costs.items():
            for orientation in (0, 1):
                nxt = ((state << 1) | orientation) & mask

                def orientation_of(block, nxt=nxt, k=k):
                    return (nxt >> (k - block)) & 1

                total = value + sum(f.cost(orientation_of) for f in ending[k])
                if total < new_costs.get(nxt, INF):
                    new_costs[nxt] = total
                    new_back[nxt] = state
        costs = new_costs
        back.append(new_back)
    best = min(costs, key=lambda s: (costs[s], s))
    orientations = [0] * n_blocks
    state = best
    for k in range(n_blocks - 1, 0, -1):
        orientations[k] = state & 1
        state = back[k][state]
    return orientations, window


def phased_cost(matrix: FragmentMatrix, haplotype: Sequence[int], block_of) -> float:
    """Full-matrix wMEC counting only alleles at variants inside a block.

    Variants no selected read covers are reported unphased, so their alleles
    in the subset's solution are arbitrary and are not scored.
    """
    total = 0.0
    for fragment in matrix:
        side0 = side1 = 0.0
        for position, allele, weight in fragment.observed():
            if block_of[position] is None:
                continue
            if allele != haplotype[position]:
                side0 += weight
            else:
                side1 += weight
        total += min(side0, side1)
    return total


class _UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left: int, right: int) -> None:
        a, b = self.find(left), self.find(right)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def orient_blocks(
    matrix: FragmentMatrix,
    subset: FragmentMatrix,
    solution: WMECSolution,
    *,
    max_window: int = 12,
) -> TwoLevelResult:
    """Orient the phase blocks of `subset`'s solution by every read of `matrix`.

    `solution` is the exact solution of `subset`, over all of `matrix`'s
    variants. Reads calling alleles in blocks more than `max_window` blocks
    apart are left out of the orientation and counted in `dropped_reads`.
    """
    if subset.n_variants != matrix.n_variants:
        raise ValueError("subset and matrix must share the variants")
    block_of, n_blocks = _block_order(subset)
    return orient_segments(matrix, solution.haplotype, block_of, n_blocks,
                           max_window=max_window)


def orient_segments(
    matrix: FragmentMatrix,
    haplotype: Sequence[int],
    segment_of: Sequence[Optional[int]],
    n_segments: int,
    *,
    max_window: int = 12,
) -> TwoLevelResult:
    """Orient given segments of `haplotype` by every read of `matrix`, exactly.

    `segment_of[j]` is the index, in chromosome order, of the segment holding
    variant `j`, or None where the variant is unphased. `orient_blocks` calls
    this with the blocks of a selected subset; stitched acceptance calls it
    with its windows cut further at the support's block breaks.
    """
    if max_window < 1:
        raise ValueError("max_window must be positive")
    if len(segment_of) != matrix.n_variants or len(haplotype) != matrix.n_variants:
        raise ValueError("haplotype and segments must cover the matrix's variants")
    haplotype = tuple(haplotype)
    block_of, n_blocks = list(segment_of), n_segments
    factors, dropped = _factors(matrix, haplotype, block_of, max_window)
    orientations, window = _solve_chain(n_blocks, factors)
    oriented = tuple(
        allele ^ orientations[block_of[j]] if block_of[j] is not None else allele
        for j, allele in enumerate(haplotype)
    )
    joined = _UnionFind(n_blocks)
    for factor in factors:
        for block in factor.blocks[1:]:
            joined.union(factor.blocks[0], block)
    phasing = Phasing(
        alleles=tuple(oriented[j] if block_of[j] is not None else None
                      for j in range(len(oriented))),
        blocks=tuple(joined.find(block_of[j]) if block_of[j] is not None else None
                     for j in range(len(oriented))),
    )
    return TwoLevelResult(
        haplotype=oriented,
        phasing=phasing,
        orientations=tuple(orientations),
        cost=phased_cost(matrix, oriented, block_of),
        subset_blocks=n_blocks,
        window=window,
        dropped_reads=dropped,
    )
