"""Paired comparison of two phasings at the same boundaries.

M5 showed that a ratio of switch-error rates misbehaves when two pipelines
share their switch errors: the ratio becomes a ratio of pair counts, and its
interval can exclude 1 for that reason alone. The comparison here is paired
by boundary instead. A pair is two consecutive variants that the truth and
both phasings phase, lying in one block of each phasing. At each pair each
phasing is right or wrong about the relation. Pairs where both are right or
both are wrong carry no information about which is better and drop out;
the discordant pairs are tested with an exact two-sided sign test
(McNemar's test in its exact form).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Sequence

from ._types import Phasing


@dataclass(frozen=True)
class PairedResult:
    pairs: int
    only_first_wrong: int
    only_second_wrong: int
    both_wrong: int
    p_value: float

    @property
    def difference(self) -> float:
        """First's error rate minus second's, on the common pairs."""
        if not self.pairs:
            return 0.0
        return (self.only_first_wrong - self.only_second_wrong) / self.pairs


def common_pairs(
    truth: Sequence[Optional[int]], first: Phasing, second: Phasing
) -> list[tuple[int, bool, bool]]:
    """`(second variant of the pair, first wrong, second wrong)` for every common pair."""
    if not len(truth) == len(first) == len(second):
        raise ValueError("truth and phasings must cover the same variants")
    usable = [j for j in range(len(truth))
              if truth[j] is not None and first.alleles[j] is not None
              and second.alleles[j] is not None]
    rows = []
    for p, q in zip(usable, usable[1:]):
        if first.blocks[p] != first.blocks[q] or second.blocks[p] != second.blocks[q]:
            continue
        relation = truth[p] ^ truth[q]
        rows.append((q, (first.alleles[p] ^ first.alleles[q]) != relation,
                     (second.alleles[p] ^ second.alleles[q]) != relation))
    return rows


def sign_test(successes: int, trials: int) -> float:
    """Exact two-sided binomial test of `successes` in `trials` against one half."""
    if trials == 0:
        return 1.0
    low = min(successes, trials - successes)
    tail = sum(math.comb(trials, k) for k in range(low + 1)) / 2 ** trials
    return min(1.0, 2 * tail)


def paired_comparison(
    truth: Sequence[Optional[int]], first: Phasing, second: Phasing
) -> PairedResult:
    rows = common_pairs(truth, first, second)
    only_first = sum(1 for _, a, b in rows if a and not b)
    only_second = sum(1 for _, a, b in rows if b and not a)
    both = sum(1 for _, a, b in rows if a and b)
    return PairedResult(len(rows), only_first, only_second, both,
                        sign_test(only_first, only_first + only_second))


def window_discordance(
    rows: Sequence[tuple[int, bool, bool]], positions: Sequence[int], window: int
) -> tuple[tuple[int, int, int], ...]:
    """Per window of `window` bp: `(pairs, only first wrong, only second wrong)`.

    Windows start at the first variant's position, as in `phasebudget.bootstrap`,
    and a pair belongs to the window of its second variant. The tuples feed
    `bootstrap.statistic_interval` for an interval on the difference.
    """
    if window <= 0:
        raise ValueError("window must be positive")
    if not len(positions):
        return ()
    first = positions[0]
    n = (positions[-1] - first) // window + 1
    counts = [[0, 0, 0] for _ in range(n)]
    for q, a, b in rows:
        k = (positions[q] - first) // window
        counts[k][0] += 1
        counts[k][1] += a and not b
        counts[k][2] += b and not a
    return tuple(tuple(c) for c in counts)
