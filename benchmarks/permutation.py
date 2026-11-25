"""A window-level sign-flip test beside the pair-level sign test.

The exact sign test of `phasebudget.paired` treats every discordant pair as
an independent coin. Switch errors cluster in space, so neighbouring
discordant pairs are not independent, and the sign test can be too
confident. The test here keeps each window's discordant pairs together and
randomizes only which pipeline is called "first" in each window.
"""

from __future__ import annotations

import random
from collections import Counter
from dataclasses import dataclass
from typing import Sequence

EXACT_LIMIT = 20


@dataclass(frozen=True)
class FlipResult:
    difference: int
    """Pairs where only the first is wrong minus pairs where only the second is, summed."""
    windows: int
    """Windows with a non-zero difference; the others do not move under a flip."""
    p_value: float
    exact: bool


def window_sign_flip_test(
    windows: Sequence[tuple[int, int, int]], resamples: int = 20000, seed: int = 0
) -> FlipResult:
    """Two-sided window sign-flip test of the paired difference in switch errors.

    `windows` is the output of `paired.window_discordance`: per window
    `(pairs, only first wrong, only second wrong)`. Let `d_k` be `only first
    - only second` in window `k` and `D` the sum of the `d_k`. The null
    hypothesis is that within each window the labels "first" and "second"
    are exchangeable, independently across windows, so each `d_k` is as
    likely as `-d_k`. The p-value is the probability, over independent
    uniform signs `s_k`, that `|sum s_k d_k| >= |D|`.

    Windows with `d_k = 0` contribute nothing whatever their sign. With at
    most `EXACT_LIMIT` non-zero windows the null distribution is computed
    exactly (by convolving the two-point distributions, which gives the same
    counts as enumerating all sign patterns). Otherwise `resamples` random
    sign patterns are drawn with `random.Random(seed)`, and the p-value is
    `(1 + hits) / (1 + resamples)`, which is never zero.

    When every window holds at most one discordant pair this is the exact
    sign test of `paired.sign_test` on the pooled counts.
    """
    if resamples < 1:
        raise ValueError("resamples must be positive")
    d = [first - second for _, first, second in windows if first != second]
    total = sum(d)
    target = abs(total)
    if len(d) <= EXACT_LIMIT:
        counts = Counter({0: 1})
        for x in d:
            step = Counter()
            for s, c in counts.items():
                step[s + x] += c
                step[s - x] += c
            counts = step
        hits = sum(c for s, c in counts.items() if abs(s) >= target)
        return FlipResult(total, len(d), hits / 2 ** len(d), True)
    rng = random.Random(seed)
    hits = 0
    for _ in range(resamples):
        bits = rng.getrandbits(len(d))
        flipped = sum(x if bits >> k & 1 else -x for k, x in enumerate(d))
        hits += abs(flipped) >= target
    return FlipResult(total, len(d), (1 + hits) / (1 + resamples), False)
