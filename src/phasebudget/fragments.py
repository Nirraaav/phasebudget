"""The fragment matrix and the quantities every other module reads off it."""

from __future__ import annotations

from typing import Iterable, Optional, Sequence

from ._types import Allele, Fragment
from ._validation import check_fragment, check_haplotype, check_unique_ids


class FragmentMatrix:
    """An ordered collection of fragments over `n_variants` ordered variants.

    Fragment order is the order given and is the order of `assignment` in a
    `WMECSolution`. `n_variants` is kept across `subset`, so variant indices
    mean the same thing in a selected subset as in the full matrix.
    """

    __slots__ = ("_fragments", "_n_variants", "_active", "_coverage")

    def __init__(
        self,
        fragments: Iterable[Fragment],
        n_variants: Optional[int] = None,
    ) -> None:
        ordered = tuple(fragments)
        check_unique_ids(ordered)
        implied = max((f.end for f in ordered), default=0)
        if n_variants is None:
            n_variants = implied
        elif n_variants < implied:
            raise ValueError(
                f"n_variants={n_variants} is smaller than the largest span end {implied}"
            )
        for fragment in ordered:
            check_fragment(fragment, n_variants)
        self._fragments = ordered
        self._n_variants = n_variants
        self._active: Optional[tuple[tuple[int, ...], ...]] = None
        self._coverage: Optional[tuple[int, ...]] = None

    @property
    def fragments(self) -> tuple[Fragment, ...]:
        return self._fragments

    @property
    def n_variants(self) -> int:
        return self._n_variants

    @property
    def read_ids(self) -> tuple[int, ...]:
        return tuple(f.read_id for f in self._fragments)

    def __len__(self) -> int:
        return len(self._fragments)

    def __iter__(self):
        return iter(self._fragments)

    def __getitem__(self, index: int) -> Fragment:
        return self._fragments[index]

    def index_of(self, read_id: int) -> int:
        for index, fragment in enumerate(self._fragments):
            if fragment.read_id == read_id:
                return index
        raise KeyError(f"no read with id {read_id}")

    def active(self, position: int) -> tuple[int, ...]:
        """Indices of the fragments whose span contains `position`."""
        if self._active is None:
            buckets: list[list[int]] = [[] for _ in range(self._n_variants)]
            for index, fragment in enumerate(self._fragments):
                for j in range(fragment.start, fragment.end):
                    buckets[j].append(index)
            self._active = tuple(tuple(bucket) for bucket in buckets)
        return self._active[position]

    def coverage(self) -> tuple[int, ...]:
        """Active coverage at every variant, counting gaps inside a span."""
        if self._coverage is None:
            counts = [0] * self._n_variants
            for fragment in self._fragments:
                for j in range(fragment.start, fragment.end):
                    counts[j] += 1
            self._coverage = tuple(counts)
        return self._coverage

    def max_coverage(self) -> int:
        return max(self.coverage(), default=0)

    def subset(self, read_ids: Iterable[int]) -> "FragmentMatrix":
        """The matrix over the given reads, keeping the variant coordinates."""
        wanted = set(read_ids)
        unknown = wanted - set(self.read_ids)
        if unknown:
            raise KeyError(f"unknown read ids: {sorted(unknown)}")
        kept = tuple(f for f in self._fragments if f.read_id in wanted)
        return FragmentMatrix(kept, self._n_variants)

    def spanned_boundaries(self) -> tuple[bool, ...]:
        """For each `j` in `[0, n_variants-1)`, whether a read spans the gap."""
        spanned = [False] * max(self._n_variants - 1, 0)
        for fragment in self._fragments:
            for j in range(fragment.start, min(fragment.last, len(spanned))):
                spanned[j] = True
        return tuple(spanned)

    def blocks(self) -> tuple[Optional[int], ...]:
        """Block identifier per variant, `None` where no fragment calls it.

        Two adjacent covered variants are in the same block when some fragment
        spans the boundary between them. Blocks are numbered from zero in
        variant order. A variant inside a fragment's span that no fragment
        calls (every fragment has a gap there) keeps its block's continuity
        but is `None` itself: nothing constrains its allele, so it is not
        phased. Coverage, which the cap counts, still includes it.
        """
        coverage = self.coverage()
        spanned = self.spanned_boundaries()
        called = [False] * self._n_variants
        for fragment in self._fragments:
            for j, _, _ in fragment.observed():
                called[j] = True
        labels: list[Optional[int]] = [None] * self._n_variants
        next_id = 0
        current: Optional[int] = None
        for j in range(self._n_variants):
            if coverage[j] == 0:
                current = None
                continue
            if current is None or j == 0 or not spanned[j - 1]:
                current = next_id
                next_id += 1
            if called[j]:
                labels[j] = current
        return tuple(labels)

    def fragment_cost(self, index: int, haplotype: Sequence[int], side: int) -> float:
        """Weighted mismatches of one fragment against one haplotype side."""
        fragment = self._fragments[index]
        total = 0.0
        for position, allele, weight in fragment.observed():
            expected = haplotype[position] ^ side
            if allele != expected:
                total += weight
        return total

    def best_assignment(self, haplotype: Sequence[int]) -> tuple[int, ...]:
        """Cheapest haplotype side per fragment; ties go to side 0."""
        check_haplotype(haplotype, self._n_variants)
        sides = []
        for index in range(len(self._fragments)):
            cost0 = self.fragment_cost(index, haplotype, 0)
            cost1 = self.fragment_cost(index, haplotype, 1)
            sides.append(0 if cost0 <= cost1 else 1)
        return tuple(sides)

    def cost(
        self,
        haplotype: Sequence[int],
        assignment: Optional[Sequence[int]] = None,
    ) -> float:
        """Weighted MEC cost of a haplotype, minimized over read assignment."""
        check_haplotype(haplotype, self._n_variants)
        if assignment is None:
            return sum(
                min(
                    self.fragment_cost(index, haplotype, 0),
                    self.fragment_cost(index, haplotype, 1),
                )
                for index in range(len(self._fragments))
            )
        if len(assignment) != len(self._fragments):
            raise ValueError(
                f"assignment has {len(assignment)} entries for "
                f"{len(self._fragments)} fragments"
            )
        return sum(
            self.fragment_cost(index, haplotype, side)
            for index, side in enumerate(assignment)
        )

    @classmethod
    def from_rows(
        cls,
        rows: Sequence[str],
        weights: Optional[Sequence[Sequence[float]]] = None,
        read_ids: Optional[Sequence[int]] = None,
    ) -> "FragmentMatrix":
        """Build a matrix from strings of `0`, `1`, and `-`.

        Every row has the same length, which becomes `n_variants`. Leading and
        trailing gaps are dropped from each span. A row with no called allele
        is rejected, since it constrains nothing and would still consume
        coverage.
        """
        if not rows:
            return cls((), 0)
        width = len(rows[0])
        if any(len(row) != width for row in rows):
            raise ValueError("all rows must have the same length")
        if weights is not None and len(weights) != len(rows):
            raise ValueError("one weight row is needed per allele row")
        if read_ids is not None and len(read_ids) != len(rows):
            raise ValueError("one read id is needed per allele row")
        fragments = []
        for index, row in enumerate(rows):
            called = [i for i, symbol in enumerate(row) if symbol != "-"]
            if not called:
                raise ValueError(f"row {index} has no called allele")
            start, stop = called[0], called[-1] + 1
            alleles: list[Allele] = []
            row_weights: list[float] = []
            for position in range(start, stop):
                symbol = row[position]
                if symbol == "-":
                    alleles.append(None)
                    row_weights.append(0.0)
                elif symbol in "01":
                    alleles.append(int(symbol))
                    row_weights.append(
                        1.0 if weights is None else float(weights[index][position])
                    )
                else:
                    raise ValueError(f"row {index}: unexpected symbol {symbol!r}")
            fragments.append(
                Fragment(
                    read_id=index if read_ids is None else int(read_ids[index]),
                    start=start,
                    alleles=tuple(alleles),
                    weights=tuple(row_weights),
                )
            )
        return cls(fragments, width)

    def to_rows(self) -> tuple[str, ...]:
        """The inverse of `from_rows`, dropping weights."""
        rows = []
        for fragment in self._fragments:
            row = ["-"] * self._n_variants
            for position, allele, _ in fragment.observed():
                row[position] = str(allele)
            rows.append("".join(row))
        return tuple(rows)
