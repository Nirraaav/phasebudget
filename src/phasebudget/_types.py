"""Immutable result objects and the fragment record."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, Optional

Allele = Optional[int]


@dataclass(frozen=True)
class Fragment:
    """One sequencing read reduced to its alleles at heterozygous variants.

    `start` is the index of the first covered variant. `alleles` is indexed by
    offset from `start`; `None` is an internal gap, a variant inside the read's
    span with no called allele. The first and last entries are never `None`, so
    the covered interval is exactly `[start, end)`.
    """

    read_id: int
    start: int
    alleles: tuple[Allele, ...]
    weights: tuple[float, ...]

    @property
    def end(self) -> int:
        """One past the last variant in the read's span."""
        return self.start + len(self.alleles)

    @property
    def last(self) -> int:
        """Index of the last variant in the read's span."""
        return self.end - 1

    def __len__(self) -> int:
        return len(self.alleles)

    def covers(self, position: int) -> bool:
        """True when `position` lies in the span, gap or not."""
        return self.start <= position < self.end

    def spans(self, boundary: int) -> bool:
        """True when the read covers both sides of the gap after `boundary`."""
        return self.start <= boundary and self.end > boundary + 1

    def allele_at(self, position: int) -> Allele:
        if not self.covers(position):
            return None
        return self.alleles[position - self.start]

    def weight_at(self, position: int) -> float:
        if not self.covers(position):
            return 0.0
        return self.weights[position - self.start]

    def observed(self) -> Iterator[tuple[int, int, float]]:
        """Yield `(position, allele, weight)` for called alleles only."""
        for offset, (allele, weight) in enumerate(zip(self.alleles, self.weights)):
            if allele is not None:
                yield self.start + offset, allele, weight

    @property
    def n_observed(self) -> int:
        return sum(1 for allele in self.alleles if allele is not None)

    @property
    def total_weight(self) -> float:
        return sum(w for a, w in zip(self.alleles, self.weights) if a is not None)


@dataclass(frozen=True)
class WMECSolution:
    """A weighted minimum error correction solution.

    `haplotype` is the first haplotype; the second is its bitwise complement,
    which is the diploid all-heterozygous model described in
    `docs/design/problem-statement.md`. `assignment[i]` is the haplotype index
    of the fragment at position `i` of the matrix that produced the solution.
    """

    cost: float
    haplotype: tuple[int, ...]
    assignment: tuple[int, ...]


@dataclass(frozen=True)
class Margin:
    """The MEC margin of one adjacent phase decision.

    `boundary` is `j`, the gap between variants `j` and `j+1`. `relation` is
    `h[j] ^ h[j+1]` in the optimal solution. `margin` is the increase in
    optimal cost forced by reversing that relation, and is zero exactly when no
    fragment spans the boundary.
    """

    boundary: int
    relation: int
    margin: float
    spanned: bool


@dataclass(frozen=True)
class SelectionResult:
    """The outcome of one coverage-constrained read selection."""

    selected: tuple[int, ...]
    objective: float
    method: str
    coverage: int


@dataclass(frozen=True)
class Phasing:
    """A per-variant phase call with block structure.

    `alleles[j]` is `None` where the variant is unphased. `blocks[j]` is the
    identifier of the phase block containing variant `j`, and `None` where the
    variant is unphased. Phase is comparable only within a block.
    """

    alleles: tuple[Allele, ...]
    blocks: tuple[Optional[int], ...]

    def __len__(self) -> int:
        return len(self.alleles)
