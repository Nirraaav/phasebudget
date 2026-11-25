"""Posterior margins: the sum-product counterpart of the MEC margin.

`docs/design/posterior-margins.md` defines them. The forward and backward
passes are those of `phasebudget.dp` on the same states, with `(min, +)`
replaced by `(logsumexp, +)` over minus the cost. The oracle is
`exact.bruteforce_posterior_margins`.
"""

from __future__ import annotations

import math
from typing import Iterable

from . import _native
from ._types import Margin
from .dp import CoverageDP
from .fragments import FragmentMatrix

NEG_INF = float("-inf")


def logsumexp(values: Iterable[float]) -> float:
    """`log(sum(exp(v)))`, stable, and `-inf` for an empty or all `-inf` input."""
    values = list(values)
    top = max(values, default=NEG_INF)
    if top == NEG_INF:
        return NEG_INF
    return top + math.log(sum(math.exp(v - top) for v in values))


def _grouped(table, bits, dp, *, by_allele: bool):
    """Log-sums of a forward or backward table by shared-read code."""
    groups: dict = {}
    for state, value in enumerate(table):
        if value == NEG_INF:
            continue
        code = dp._project(state >> 1, bits)
        key = (code, state & 1) if by_allele else code
        groups.setdefault(key, []).append(value)
    return {key: logsumexp(items) for key, items in groups.items()}


def boundary_log_sums(
    matrix: FragmentMatrix, *, backend: str = "auto"
) -> list[tuple[float, float]]:
    """`(log Z_0(b), log Z_1(b))` for every boundary `b`.

    `backend` is "python", "native", or "auto", which is native when the
    library is built. The two agree to rounding.
    """
    if backend not in ("auto", "python", "native"):
        raise ValueError(f"unknown backend {backend!r}")
    dp = CoverageDP(matrix, backend="python")
    n = matrix.n_variants
    if n < 2:
        return []
    if backend == "native" or (backend == "auto" and _native.available()):
        offsets, active, c0, c1 = dp.native_inputs()
        return _native.posterior(n, offsets, active, c0, c1, len(matrix))
    own = [[-cost for cost in dp._state_costs(j)] for j in range(n)]

    forward = [own[0]]
    for j in range(1, n):
        shared = dp._shared(j - 1, j)
        left = [pair[0] for pair in shared]
        right = [pair[1] for pair in shared]
        incoming = _grouped(forward[j - 1], left, dp, by_allele=False)
        forward.append([
            value + incoming.get(dp._project(state >> 1, right), NEG_INF)
            for state, value in enumerate(own[j])
        ])

    backward = [None] * n
    backward[n - 1] = own[n - 1]
    for j in range(n - 2, -1, -1):
        shared = dp._shared(j, j + 1)
        left = [pair[0] for pair in shared]
        right = [pair[1] for pair in shared]
        outgoing = _grouped(backward[j + 1], right, dp, by_allele=False)
        backward[j] = [
            value + outgoing.get(dp._project(state >> 1, left), NEG_INF)
            for state, value in enumerate(own[j])
        ]

    result = []
    for b in range(n - 1):
        shared = dp._shared(b, b + 1)
        left = [pair[0] for pair in shared]
        right = [pair[1] for pair in shared]
        heads = _grouped(forward[b], left, dp, by_allele=True)
        tails = _grouped(backward[b + 1], right, dp, by_allele=True)
        sums = ([], [])
        for (code, allele), head in heads.items():
            for relation in (0, 1):
                tail = tails.get((code, allele ^ relation))
                if tail is not None:
                    sums[relation].append(head + tail)
        result.append((logsumexp(sums[0]), logsumexp(sums[1])))
    return result


def posterior_margins(matrix: FragmentMatrix, *, backend: str = "auto") -> tuple[Margin, ...]:
    """The posterior margin of every boundary, as `Margin` records.

    `margin` is the log posterior odds of the likelier relation, which is in
    the same units as the MEC margin, so the adaptive loop can rank either.
    """
    spanned = matrix.spanned_boundaries()
    margins = []
    for boundary, (zero, one) in enumerate(boundary_log_sums(matrix, backend=backend)):
        relation = 0 if zero >= one else 1
        gap = abs(zero - one) if zero != one else 0.0
        margins.append(
            Margin(boundary=boundary, relation=relation, margin=gap, spanned=spanned[boundary])
        )
    return tuple(margins)
