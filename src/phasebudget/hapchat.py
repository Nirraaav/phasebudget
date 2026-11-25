"""HapCHAT's read reduction: merge likely same-haplotype reads, then select greedily.

Beretta et al. 2018 (read as the 2017 bioRxiv preprint, section 3.1 and the
experimental setup) reduce coverage in two steps. Reads that very probably
come from one haplotype are merged: for two reads that share `m + x` called
sites, agreeing on `m` and disagreeing on `x`, the paper gives
`p_s = p^x (1 - p/3^m)` for the same haplotype and `p_d = p^m (1 - p/3^x)` for
different ones, with `p = 0.15`, and partitions the reads into clusters in
which every pair has `p_s / p_d >= 10^6`. Each merged site takes the Phred-
weighted majority. Then, if the cap is still exceeded, reads are selected
greedily by Phred score with ties broken at random.

What the paper leaves open is fixed here and stated in
`docs/design/published-selection-objectives.md`: clusters are formed greedily
in order of start, a merged call carries the weight of the majority minus the
minority, and a read's score is its total weight. Merged reads are new reads,
so this selector is compared on accuracy and on wMEC against the full matrix,
never on the additive objective of the original reads.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Optional

from ._types import Fragment, SelectionResult
from ._validation import check_coverage
from .fragments import FragmentMatrix

#: The per-base error probability the paper assumes for PacBio reads.
ERROR_PROBABILITY = 0.15

#: The likelihood ratio a pair must reach to share a cluster.
RATIO_THRESHOLD = 1e6


def agreement(first: Fragment, second: Fragment) -> tuple[int, int]:
    """`(m, x)`: shared called sites on which the two reads agree and disagree."""
    agree = disagree = 0
    for position in range(max(first.start, second.start), min(first.end, second.end)):
        a, b = first.allele_at(position), second.allele_at(position)
        if a is None or b is None:
            continue
        if a == b:
            agree += 1
        else:
            disagree += 1
    return agree, disagree


def _shrunk(p: float, sites: int) -> float:
    """`p / 3**sites`, as a float even when `3**sites` is not.

    Ultra-long ONT reads can share more than 646 sites, where `3**sites`
    overflows a float. Below 600 the value is computed exactly as before; above
    it the term is under 1e-286, and zero changes no log-ratio.
    """
    return p / 3 ** sites if sites < 600 else 0.0


def log_ratio(agree: int, disagree: int, p: float = ERROR_PROBABILITY) -> float:
    """`log(p_s / p_d)` with the paper's `p_s` and `p_d`."""
    same = disagree * math.log(p) + math.log1p(-_shrunk(p, agree))
    different = agree * math.log(p) + math.log1p(-_shrunk(p, disagree))
    return same - different


@dataclass(frozen=True)
class Merged:
    """The merged matrix and, per merged read, the original read ids it holds."""

    matrix: FragmentMatrix
    members: tuple[tuple[int, ...], ...]


def merge_reads(
    matrix: FragmentMatrix,
    *,
    p: float = ERROR_PROBABILITY,
    threshold: float = RATIO_THRESHOLD,
) -> Merged:
    """Cluster reads so every pair in a cluster passes the ratio, then merge each."""
    limit = math.log(threshold)
    order = sorted(range(len(matrix)), key=lambda i: (matrix[i].start, matrix[i].read_id))
    clusters: list[list[int]] = []
    cluster_end: list[int] = []
    open_clusters: list[int] = []
    for index in order:
        fragment = matrix[index]
        open_clusters = [c for c in open_clusters if cluster_end[c] > fragment.start]
        home: Optional[int] = None
        for c in open_clusters:
            if all(
                log_ratio(*agreement(matrix[other], fragment), p=p) >= limit
                for other in clusters[c]
            ):
                home = c
                break
        if home is None:
            clusters.append([index])
            cluster_end.append(fragment.end)
            open_clusters.append(len(clusters) - 1)
        else:
            clusters[home].append(index)
            cluster_end[home] = max(cluster_end[home], fragment.end)

    fragments: list[Fragment] = []
    members: list[tuple[int, ...]] = []
    for cluster in clusters:
        start = min(matrix[i].start for i in cluster)
        end = max(matrix[i].end for i in cluster)
        alleles: list[Optional[int]] = []
        weights: list[float] = []
        for position in range(start, end):
            support = [0.0, 0.0]
            for i in cluster:
                allele = matrix[i].allele_at(position)
                if allele is not None:
                    support[allele] += matrix[i].weight_at(position)
            if support[0] == support[1]:
                alleles.append(None)
                weights.append(0.0)
            else:
                winner = 0 if support[0] > support[1] else 1
                alleles.append(winner)
                weights.append(abs(support[0] - support[1]))
        called = [k for k, allele in enumerate(alleles) if allele is not None]
        if not called:
            continue
        low, high = called[0], called[-1] + 1
        fragments.append(
            Fragment(
                read_id=len(fragments),
                start=start + low,
                alleles=tuple(alleles[low:high]),
                weights=tuple(weights[low:high]),
            )
        )
        members.append(tuple(sorted(matrix[i].read_id for i in cluster)))
    return Merged(FragmentMatrix(fragments, matrix.n_variants), tuple(members))


def select_hapchat(
    matrix: FragmentMatrix,
    coverage: int,
    *,
    seed: int = 0,
    p: float = ERROR_PROBABILITY,
    threshold: float = RATIO_THRESHOLD,
) -> tuple[Merged, SelectionResult]:
    """Merge, then keep merged reads greedily by total weight under the cap.

    The returned selection names reads of the merged matrix, not of `matrix`.
    Ties in total weight are broken at random with `seed`, as the paper does.
    """
    check_coverage(coverage)
    merged = merge_reads(matrix, p=p, threshold=threshold)
    rng = random.Random(seed)
    scored = [(sum(f.weights), rng.random(), f) for f in merged.matrix]
    scored.sort(key=lambda item: (-item[0], item[1]))
    counts = [0] * matrix.n_variants
    selected: list[int] = []
    for score, _, fragment in scored:
        if score <= 0.0:
            continue
        if any(counts[j] >= coverage for j in range(fragment.start, fragment.end)):
            continue
        for j in range(fragment.start, fragment.end):
            counts[j] += 1
        selected.append(fragment.read_id)
    selected.sort()
    return merged, SelectionResult(
        selected=tuple(selected),
        objective=float("nan"),
        method="hapchat",
        coverage=coverage,
    )
