"""Read support for a phasing's relations, and a threshold with a risk guarantee.

`docs/followup/error-controlled-phasing.md` states the problem. For each pair
of consecutive phased variants in one block, the reads calling both variants
either agree with the called relation (same or opposite haplotype) or not.
The score of a pair is the Wilson lower bound on the agreeing share, with
ties broken by the number of reads. A threshold on the score decides which
pairs are trusted; `learn_then_test` chooses it so that, over units of the
genome treated as independent, the probability that a unit holds a wrong
trusted pair is at most `beta` with confidence `1 - delta`.

The phasing can come from any phaser: only its alleles and blocks are read.
"""

from __future__ import annotations

import math
from bisect import bisect_left
from dataclasses import dataclass
from typing import Iterable, Optional, Sequence

from ._types import Phasing
from .fragments import FragmentMatrix


@dataclass(frozen=True)
class PairSupport:
    """Reads calling both variants of a consecutive phased pair `(first, second)`."""

    first: int
    second: int
    agree: int
    disagree: int

    @property
    def reads(self) -> int:
        return self.agree + self.disagree


def phased_pairs(phasing: Phasing) -> tuple[tuple[int, int], ...]:
    """Consecutive phased variants that share a block, in position order."""
    phased = [j for j, a in enumerate(phasing.alleles) if a is not None]
    return tuple((a, b) for a, b in zip(phased, phased[1:])
                 if phasing.blocks[a] is not None and phasing.blocks[a] == phasing.blocks[b])


def pair_support(matrix: FragmentMatrix, phasing: Phasing) -> tuple[PairSupport, ...]:
    """Agreeing and disagreeing reads for every pair of `phased_pairs(phasing)`."""
    if len(phasing.alleles) != matrix.n_variants:
        raise ValueError("the phasing and the matrix must cover the same variants")
    pairs = phased_pairs(phasing)
    following = dict(pairs)
    agree = dict.fromkeys(following, 0)
    disagree = dict.fromkeys(following, 0)
    for fragment in matrix:
        calls = {fragment.start + k: allele for k, allele in enumerate(fragment.alleles)
                 if allele is not None}
        for j, allele in calls.items():
            k = following.get(j)
            if k is None or k not in calls:
                continue
            read_same = allele == calls[k]
            called_same = phasing.alleles[j] == phasing.alleles[k]
            if read_same == called_same:
                agree[j] += 1
            else:
                disagree[j] += 1
    return tuple(PairSupport(a, b, agree[a], disagree[a]) for a, b in pairs)


def wilson_lower(successes: int, trials: int, z: float = 2.0) -> float:
    """Lower end of the Wilson score interval; -1 when there are no trials."""
    if trials < 0 or not 0 <= successes <= trials:
        raise ValueError("need 0 <= successes <= trials")
    if trials == 0:
        return -1.0
    p = successes / trials
    centre = p + z * z / (2 * trials)
    spread = z * math.sqrt(p * (1 - p) / trials + z * z / (4 * trials * trials))
    return (centre - spread) / (1 + z * z / trials)


def score(support: PairSupport, z: float = 2.0) -> tuple[float, int]:
    """The frozen score: (Wilson lower bound on the agreeing share, reads). Higher is trusted."""
    return wilson_lower(support.agree, support.reads, z), support.reads


def binomial_cdf(k: int, n: int, p: float) -> float:
    """P(X <= k) for X ~ Binomial(n, p), summed in log space."""
    if not 0 < p < 1:
        raise ValueError("p must lie strictly between 0 and 1")
    if k < 0:
        return 0.0
    if k >= n:
        return 1.0
    total = 0.0
    for i in range(k + 1):
        log_term = (math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
                    + i * math.log(p) + (n - i) * math.log1p(-p))
        total += math.exp(log_term)
    return min(1.0, total)


@dataclass(frozen=True)
class Units:
    """Per unit, the best score of any pair and of any wrong pair (None if none is wrong)."""

    best: tuple
    best_wrong: tuple

    @classmethod
    def from_pairs(cls, rows: Iterable[tuple[object, tuple, bool]]) -> "Units":
        """`rows` are `(unit id, score, wrong)`."""
        best: dict = {}
        worst: dict = {}
        for unit, key, wrong in rows:
            if unit not in best or key > best[unit]:
                best[unit] = key
            if wrong and (unit not in worst or key > worst[unit]):
                worst[unit] = key
        return cls(tuple(sorted(best.values())), tuple(sorted(worst.values())))

    def at(self, threshold: tuple) -> tuple[int, int]:
        """`(units holding a trusted pair, units holding a trusted wrong pair)`."""
        return (len(self.best) - bisect_left(self.best, threshold),
                len(self.best_wrong) - bisect_left(self.best_wrong, threshold))


def learn_then_test(units: Units, grid: Sequence[tuple], beta: float,
                    delta: float) -> Optional[tuple]:
    """The loosest threshold of `grid` (ordered strict to loose) certified at `beta`.

    Fixed-sequence testing: at each threshold, H0 says the share of units
    with a wrong trusted pair exceeds `beta`. It is rejected when the exact
    binomial p-value `P(X <= k)` under `beta` is at most `delta`, and the
    scan stops at the first threshold where it is not. Returns None when the
    first threshold is not certified.
    """
    if not 0 < delta < 1:
        raise ValueError("delta must lie strictly between 0 and 1")
    chosen = None
    for threshold in grid:
        n, k = units.at(threshold)
        if n == 0 or binomial_cdf(k, n, beta) > delta:
            break
        chosen = threshold
    return chosen


def coverage_grid(scores: Sequence[tuple], start: float = 0.5, step: float = 0.0025
                  ) -> tuple[tuple, ...]:
    """Thresholds at which `start`, `start + step`, ..., 1 of the given pairs are trusted."""
    ranked = sorted(scores, reverse=True)
    if not ranked:
        return ()
    count = int(round((1 - start) / step))
    shares = [start + i * step for i in range(count + 1)]
    return tuple(ranked[max(0, min(len(ranked), int(round(s * len(ranked)))) - 1)]
                 for s in shares)


def critical_count(n: int, beta: float, delta: float) -> int:
    """Largest `k` with `P(Bin(n, beta) <= k) <= delta`, or -1: the most wrong units a
    threshold with `n` trusted units may hold and still pass `learn_then_test`'s test."""
    from statistics import NormalDist

    def pmf(i: int) -> float:
        return math.exp(math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
                        + i * math.log(beta) + (n - i) * math.log1p(-beta))

    # Start at the normal approximation, then walk to the exact boundary.
    sd = math.sqrt(n * beta * (1 - beta))
    k = min(n, max(-1, int(n * beta + NormalDist().inv_cdf(delta) * sd)))
    c = binomial_cdf(k, n, beta)
    while k < n and c + pmf(k + 1) <= delta:
        k += 1
        c += pmf(k)
    while k >= 0 and c > delta:
        c -= pmf(k)
        k -= 1
    return k


def certification_probability(n: int, r: float, beta: float, delta: float) -> float:
    """`P(Bin(n, r) <= critical_count(n, beta, delta))`: the chance that a threshold with `n`
    independent units, each wrong with probability `r`, passes. `docs/design/l4-calibration-planner.md`."""
    k = critical_count(n, beta, delta)
    if k < 0:
        return 0.0
    if r <= 0:
        return 1.0
    return binomial_cdf(k, n, r) if r < 1 else float(k >= n)


def units_needed(r: float, beta: float, delta: float, power: float,
                 limit: int = 10**8) -> Optional[int]:
    """Units for certification probability `power`: the `n` found by bisection, to 0.1%,
    between the `r = 0` floor and a doubling bound, on the predicate "probability at least
    `power` at `n` and at 50 evenly spaced points of `[n, 2n]`". The window guards against
    the sawtooth of `critical_count`. None if no `n` up to `limit` qualifies."""
    if r >= beta:
        return None

    def ok(n: int) -> bool:
        return all(certification_probability(n + (n * i) // 50, r, beta, delta) >= power
                   for i in range(51))

    lo = math.ceil(math.log(1 / delta) / -math.log1p(-beta))  # nothing below certifies
    if ok(lo):
        return lo
    hi = 2 * lo
    while not ok(hi):
        lo, hi = hi, 2 * hi
        if hi > limit:
            return None
    while hi - lo > max(1, lo // 1000):
        mid = (lo + hi) // 2
        lo, hi = (lo, mid) if ok(mid) else (mid, hi)
    return hi


def learn_then_test_by_group(rows: Iterable[tuple[object, object, tuple, bool]], beta: float,
                             delta: float, start: float = 0.5, step: float = 0.0025
                             ) -> dict:
    """`learn_then_test` run on each group's own rows; None for a group that certifies nothing.

    `rows` are `(group, unit id, score, wrong)`. Each group gets its own grid,
    `coverage_grid` over its scores, and its own units, so one group's errors
    never move another group's threshold. `docs/design/stratified-confidence.md`
    states the procedure and its guarantee.
    """
    grouped: dict = {}
    for group, unit, key, wrong in rows:
        grouped.setdefault(group, []).append((unit, key, wrong))
    return {group: learn_then_test(Units.from_pairs(members),
                                   coverage_grid([key for _, key, _ in members], start, step),
                                   beta, delta)
            for group, members in grouped.items()}


@dataclass(frozen=True)
class ChainCounts:
    """Relations between labelled variants of one block, and which of them a threshold trusts.

    A relation `(i, j)` is trusted when every consecutive link from `i` to `j`
    is trusted, and wrong when the labels of `i` and `j` differ.
    `trusted_span` sums, over trusted relations, the number of units their
    chain can cross, counted as `bins[j] - bins[i] + 1`.
    """

    pairs: int
    wrong: int
    trusted: int
    wrong_trusted: int
    trusted_span: int


def chain_counts(links: Sequence[bool], labels: Sequence[Optional[int]],
                 bins: Sequence[int]) -> ChainCounts:
    """Count the relations of `docs/design/chain-relations.md` in one block, in linear time.

    `labels[k]` is the phaser's first-haplotype allele XOR the truth's at
    variant `k`, or None where the truth has no phase; such variants are
    links in chains but never endpoints. `links[k]` says whether the link
    from variant `k` to `k + 1` is trusted. `bins` must not decrease.
    """
    if len(links) != max(0, len(labels) - 1) or len(bins) != len(labels):
        raise ValueError("need one link between each two consecutive variants")
    totals = [0, 0, 0, 0, 0]

    def close(counts: list, pair_index: int, wrong_index: int) -> None:
        n = counts[0] + counts[1]
        totals[pair_index] += n * (n - 1) // 2
        totals[wrong_index] += counts[0] * counts[1]

    whole = [0, 0]
    seg, seg_bins = [0, 0], [0, 0]  # labelled variants in the trusted segment; sum and count of bins
    for k, label in enumerate(labels):
        if k and not links[k - 1]:
            close(seg, 2, 3)
            seg, seg_bins = [0, 0], [0, 0]
        if label is None:
            continue
        # Each earlier labelled variant `i` of the segment adds bins[k] - bins[i] + 1.
        totals[4] += seg_bins[1] * (bins[k] + 1) - seg_bins[0]
        seg[label] += 1
        seg_bins[0] += bins[k]
        seg_bins[1] += 1
        whole[label] += 1
    close(seg, 2, 3)
    close(whole, 0, 1)
    return ChainCounts(*totals)
