"""Input checking. Every public entry point validates before it computes."""

from __future__ import annotations

from typing import Iterable, Optional, Sequence

from ._types import Fragment


class CoverageLimitExceeded(ValueError):
    """Raised when an exponential-state algorithm is given too much coverage."""


def check_fragment(fragment: Fragment, n_variants: Optional[int] = None) -> None:
    if not isinstance(fragment, Fragment):
        raise TypeError(f"expected Fragment, got {type(fragment).__name__}")
    if not isinstance(fragment.read_id, int) or isinstance(fragment.read_id, bool):
        raise TypeError("read_id must be an int")
    if fragment.read_id < 0:
        raise ValueError(f"read_id must be non-negative, got {fragment.read_id}")
    if fragment.start < 0:
        raise ValueError(f"read {fragment.read_id}: start must be non-negative")
    if len(fragment.alleles) != len(fragment.weights):
        raise ValueError(
            f"read {fragment.read_id}: {len(fragment.alleles)} alleles but "
            f"{len(fragment.weights)} weights"
        )
    if not fragment.alleles:
        raise ValueError(f"read {fragment.read_id}: empty fragment")
    if fragment.alleles[0] is None or fragment.alleles[-1] is None:
        raise ValueError(
            f"read {fragment.read_id}: span must be trimmed to called alleles"
        )
    for allele in fragment.alleles:
        if allele is not None and allele not in (0, 1):
            raise ValueError(
                f"read {fragment.read_id}: allele {allele!r} is not 0, 1, or None"
            )
    for allele, weight in zip(fragment.alleles, fragment.weights):
        if weight < 0:
            raise ValueError(f"read {fragment.read_id}: negative weight {weight}")
        if allele is None and weight != 0.0:
            raise ValueError(
                f"read {fragment.read_id}: a gap must carry weight 0.0, got {weight}"
            )
    if n_variants is not None and fragment.end > n_variants:
        raise ValueError(
            f"read {fragment.read_id}: span ends at {fragment.end}, past the "
            f"{n_variants} variants of the matrix"
        )


def check_unique_ids(fragments: Sequence[Fragment]) -> None:
    seen: set[int] = set()
    for fragment in fragments:
        if fragment.read_id in seen:
            raise ValueError(f"duplicate read_id {fragment.read_id}")
        seen.add(fragment.read_id)


def check_haplotype(haplotype: Sequence[int], n_variants: int) -> None:
    if len(haplotype) != n_variants:
        raise ValueError(
            f"haplotype has {len(haplotype)} entries, matrix has {n_variants} variants"
        )
    for allele in haplotype:
        if allele not in (0, 1):
            raise ValueError(f"haplotype entry {allele!r} is not 0 or 1")


def check_coverage(coverage: int) -> None:
    if not isinstance(coverage, int) or isinstance(coverage, bool):
        raise TypeError("coverage must be an int")
    if coverage < 0:
        raise ValueError(f"coverage must be non-negative, got {coverage}")


def check_weights(weights: dict[int, float], read_ids: Iterable[int]) -> None:
    for read_id in read_ids:
        if read_id not in weights:
            raise ValueError(f"no weight given for read {read_id}")
        weight = weights[read_id]
        if weight != weight or weight in (float("inf"), float("-inf")):
            raise ValueError(f"read {read_id}: weight must be finite, got {weight}")
