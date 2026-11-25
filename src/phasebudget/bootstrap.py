"""Intervals for switch error on one chromosome, by resampling windows.

A real-data run has one instance, so replicates cannot give an interval.
`docs/design/real-data-plan.md` registers a block bootstrap instead: cut the
chromosome into fixed windows of base pairs, resample windows with
replacement, and recompute each statistic from the resampled switch and pair
counts. Windows are resampled jointly across pipelines, so a ratio compares
the pipelines on the same resampled genome.

A comparable pair is assigned to the window holding its second variant,
which is also where `metrics.switch_positions` records a switch. The window
counts therefore sum to `len(switch_positions)` and `comparable_pairs`, and
`tests/test_bootstrap.py` checks that against the metrics module.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Optional, Sequence

from ._types import Phasing
from .metrics import _compared_members


@dataclass(frozen=True)
class Interval:
    """A point estimate with a percentile interval."""

    estimate: float
    low: float
    high: float
    resamples: int
    undefined: int


def window_counts(
    truth: Sequence[Optional[int]],
    phasing: Phasing,
    positions: Sequence[int],
    window: int,
) -> tuple[tuple[int, int], ...]:
    """`(switches, pairs)` in each window of `window` base pairs.

    Window `k` holds positions in `[first + k * window, first + (k + 1) *
    window)`, where `first` is the first variant's position, so every call
    on the same `positions` cuts the same windows.
    """
    if window <= 0:
        raise ValueError("window must be positive")
    if len(positions) != len(phasing):
        raise ValueError("positions must give one position per variant")
    if not len(positions):
        return ()
    first = positions[0]
    n_windows = (positions[-1] - first) // window + 1
    switches = [0] * n_windows
    pairs = [0] * n_windows
    for members in _compared_members(truth, phasing):
        for left, right in zip(members, members[1:]):
            k = (positions[right] - first) // window
            pairs[k] += 1
            predicted = phasing.alleles[left] ^ phasing.alleles[right]
            if predicted != truth[left] ^ truth[right]:
                switches[k] += 1
    return tuple(zip(switches, pairs))


def _rate(counts: Sequence[tuple[int, int]], picks: Sequence[int]) -> Optional[float]:
    pairs = sum(counts[k][1] for k in picks)
    if pairs == 0:
        return None
    return sum(counts[k][0] for k in picks) / pairs


def _ratio(numerator, denominator, picks) -> Optional[float]:
    top = _rate(numerator, picks)
    bottom = _rate(denominator, picks)
    if top is None or bottom is None or bottom == 0.0:
        return None
    return top / bottom


def _percentile(values: Sequence[float], q: float) -> float:
    """The `q` quantile by linear interpolation between order statistics."""
    ordered = sorted(values)
    position = q * (len(ordered) - 1)
    low = int(position)
    high = min(low + 1, len(ordered) - 1)
    return ordered[low] + (position - low) * (ordered[high] - ordered[low])


def _interval(statistic, n_windows: int, resamples: int, seed: int, level: float) -> Interval:
    if not 0.0 < level < 1.0:
        raise ValueError("level must lie strictly between 0 and 1")
    if resamples < 1:
        raise ValueError("resamples must be positive")
    everything = range(n_windows)
    estimate = statistic(everything)
    if estimate is None:
        raise ValueError("the statistic is undefined on the full data")
    rng = random.Random(seed)
    values = []
    undefined = 0
    for _ in range(resamples):
        picks = [rng.randrange(n_windows) for _ in range(n_windows)]
        value = statistic(picks)
        if value is None:
            undefined += 1
        else:
            values.append(value)
    if not values:
        raise ValueError("the statistic is undefined on every resample")
    tail = (1.0 - level) / 2.0
    return Interval(estimate, _percentile(values, tail), _percentile(values, 1.0 - tail),
                    resamples, undefined)


def rate_interval(
    counts: Sequence[tuple[int, int]],
    *,
    resamples: int = 2000,
    seed: int = 0,
    level: float = 0.95,
) -> Interval:
    """Switch error rate with a window-bootstrap percentile interval."""
    return _interval(lambda picks: _rate(counts, picks), len(counts), resamples, seed, level)


def ratio_interval(
    numerator: Sequence[tuple[int, int]],
    denominator: Sequence[tuple[int, int]],
    *,
    resamples: int = 2000,
    seed: int = 0,
    level: float = 0.95,
) -> Interval:
    """`rate(numerator) / rate(denominator)`, with windows resampled jointly.

    Both count sequences must come from `window_counts` on the same positions
    and window. A resample where the denominator has no switches, or no
    pairs, leaves the ratio undefined; it is counted in `undefined` and left
    out of the percentiles, and the report must state that count.
    """
    if len(numerator) != len(denominator):
        raise ValueError("the two count sequences must cover the same windows")
    return _interval(lambda picks: _ratio(numerator, denominator, picks), len(numerator),
                     resamples, seed, level)


def statistic_interval(
    statistic,
    n_windows: int,
    *,
    resamples: int = 2000,
    seed: int = 0,
    level: float = 0.95,
) -> Interval:
    """Any statistic of the windows, with the same resampling as the ratios.

    `statistic` takes a sequence of window indices, repeats allowed, and
    returns a number or None where it is undefined. The rate and ratio
    intervals above are this function with their own statistics, so a
    statistic that combines several pipelines resamples the same windows for
    all of them.
    """
    return _interval(statistic, n_windows, resamples, seed, level)
