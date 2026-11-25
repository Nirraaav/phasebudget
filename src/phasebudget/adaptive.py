"""Solver-guided adaptive read selection.

One round is: select a coverage-feasible subset, solve wMEC on it, read the
margin of every adjacent phase decision off the dynamic program, raise the
weight of the reads that span the fragile boundaries, and select again. The
loop has a fixed iteration cap and keeps the best solution it has seen, so
reweighting cannot make the reported answer worse than the first round.

Solutions from different subsets are not comparable by their own objective
values, because a smaller subset has fewer reads to disagree with. Every
candidate is therefore scored by its wMEC cost against the full matrix, which
is the same yardstick for all of them.

The margins that steer the loop are measured on the selected subset, which is
the thing selection just changed. `margin_source="full"` measures them on the
whole matrix instead. That is an oracle, affordable only below the state limit,
and it exists to test whether the subset signal is the fragility of the phasing
or the thinness of the subset. `docs/design/margin-source.md` states the
confound and what distinguishes the two.

`signal` chooses what fragility means. `"margin"` is the MEC margin, the
max-product gap to the best solution with the relation reversed. `"posterior"`
is the posterior margin of `docs/design/posterior-margins.md`, the log odds of
the relation summed over every solution, which ranks switch errors better when
the reads are noisy.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Callable, Optional

from ._types import Margin, SelectionResult, WMECSolution
from ._validation import check_coverage
from .dp import CoverageDP
from .fragments import FragmentMatrix
from .posterior import posterior_margins
from .selection import read_weights, select_optimal, total_weight_score
from .twolevel import orient_segments

#: Where the loop measures the margins that decide which reads to boost.
MARGIN_SOURCES = ("subset", "full")

#: What the loop ranks boundaries by. "shuffled" is a control: the posterior
#: margins permuted at random across spanned boundaries, so the loop boosts as
#: many boundaries as the posterior would, in places the posterior did not
#: choose.
SIGNALS = ("margin", "posterior", "shuffled")

#: How the loop turns its rounds into one answer.
ACCEPTANCE_RULES = ("global", "stitched", "stitched_oriented")


def _partial_cost(fragment, haplotype, limit: int, side: int) -> float:
    """Weighted mismatches of a fragment against `haplotype` at positions below `limit`."""
    total = 0.0
    for position, allele, weight in fragment.observed():
        if position < limit and allele != haplotype[position] ^ side:
            total += weight
    return total


def stitch_haplotypes(
    matrix: FragmentMatrix, haplotypes, window: int = 25
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """`(haplotype, round chosen per window)`, by `docs/design/adaptive-acceptance.md`.

    Each window of `window` variants takes the round whose haplotype gives the
    reads lying inside it the lowest full-matrix cost; each junction takes the
    orientation that the reads spanning it prefer.
    """
    if window < 1:
        raise ValueError("window must be positive")
    n = matrix.n_variants
    if n == 0 or not haplotypes:
        return (), ()
    windows = [(a, min(a + window, n)) for a in range(0, n, window)]
    inside: list[list[int]] = [[] for _ in windows]
    crossing: list[list[int]] = [[] for _ in windows]
    for index, fragment in enumerate(matrix):
        first, last = fragment.start // window, (fragment.end - 1) // window
        if first == last:
            inside[first].append(index)
        for w in range(first + 1, last + 1):
            crossing[w].append(index)
    choice = []
    for w in range(len(windows)):
        costs = [
            sum(
                min(matrix.fragment_cost(i, h, 0), matrix.fragment_cost(i, h, 1))
                for i in inside[w]
            )
            for h in haplotypes
        ]
        choice.append(min(range(len(haplotypes)), key=lambda r: (costs[r], r)))
    stitched = list(haplotypes[choice[0]][windows[0][0]:windows[0][1]])
    for w in range(1, len(windows)):
        a, b = windows[w]
        piece = list(haplotypes[choice[w]][a:b])
        flipped = [1 - allele for allele in piece]
        scores = []
        for candidate in (piece, flipped):
            trial = stitched + candidate
            scores.append(sum(
                min(_partial_cost(matrix[i], trial, b, 0), _partial_cost(matrix[i], trial, b, 1))
                for i in crossing[w]
            ))
        stitched += piece if scores[0] <= scores[1] else flipped
    return tuple(stitched), tuple(choice)


@dataclass(frozen=True)
class AdaptiveStep:
    """What one round of the loop selected, solved, and decided to reweight."""

    iteration: int
    selection: SelectionResult
    solution: WMECSolution
    full_cost: float
    margins: tuple[Margin, ...]
    fragile: tuple[int, ...]
    boosted: tuple[int, ...]
    margin_source: str = "subset"


@dataclass(frozen=True)
class AdaptiveResult:
    """The best round of the loop, plus the trace of every round."""

    solution: WMECSolution
    selection: SelectionResult
    full_cost: float
    iterations: int
    history: tuple[AdaptiveStep, ...]
    margin_source: str = "subset"
    acceptance: str = "global"
    support: tuple[int, ...] = ()

    @property
    def selected(self) -> tuple[int, ...]:
        return self.selection.selected


def fragile_boundaries(
    margins: tuple[Margin, ...],
    *,
    quantile: float = 0.25,
    minimum: int = 1,
    threshold: Optional[float] = None,
) -> tuple[int, ...]:
    """The spanned boundaries with the smallest margins.

    With `threshold`, the rule is absolute instead: every spanned boundary
    whose margin is below it, however many that is, and none if none are.
    That is the scale-free rule; a quantile marks a quarter of a chromosome's
    boundaries fragile and boosts nearly every read.

    Unspanned boundaries are excluded: their margin is zero for a structural
    reason, and no reweighting of the reads that do not exist can change it.
    The threshold is the `quantile` order statistic of the spanned margins, so
    a run where every decision is equally firm still reports the weakest ones
    rather than none.
    """
    if not 0.0 < quantile <= 1.0:
        raise ValueError(f"quantile must lie in (0, 1], got {quantile}")
    spanned = [margin for margin in margins if margin.spanned]
    if threshold is not None:
        return tuple(sorted(m.boundary for m in spanned if m.margin < threshold))
    if not spanned:
        return ()
    ordered = sorted(spanned, key=lambda m: (m.margin, m.boundary))
    count = max(minimum, int(round(quantile * len(ordered))))
    count = min(count, len(ordered))
    return tuple(sorted(margin.boundary for margin in ordered[:count]))


def adaptive_selection(
    matrix: FragmentMatrix,
    coverage: int,
    *,
    iterations: int = 4,
    boost: float = 2.0,
    quantile: float = 0.25,
    score: Callable = total_weight_score,
    selector: Callable[..., SelectionResult] = select_optimal,
    max_active_reads: Optional[int] = None,
    margin_source: str = "subset",
    signal: str = "margin",
    error_probability: Optional[float] = None,
    acceptance: str = "global",
    window: int = 25,
    shuffle_seed: int = 0,
) -> AdaptiveResult:
    """Run the reweighting loop and return its best round.

    `boost` multiplies the weight of a read for every fragile boundary it
    spans, capped at one multiplication per round. `iterations` counts rounds
    including the first, so `iterations=1` reduces to static selection under
    the base score.

    `margin_source` chooses where fragility is measured. `"subset"` reads it
    off the dynamic program that just solved the selected subset. `"full"`
    reads it off one dynamic program over the entire matrix, which is an oracle
    for OQ-10: it raises `CoverageLimitExceeded` whenever the full matrix is
    above the state limit, which is most of the cases selection exists for.
    """
    check_coverage(coverage)
    if iterations < 1:
        raise ValueError("iterations must be at least 1")
    if boost <= 1.0:
        raise ValueError("boost must exceed 1.0 or the loop cannot move")
    if signal not in SIGNALS:
        raise ValueError(f"signal must be one of {SIGNALS}, got {signal!r}")
    if acceptance not in ACCEPTANCE_RULES:
        raise ValueError(f"acceptance must be one of {ACCEPTANCE_RULES}, got {acceptance!r}")
    threshold = None
    if error_probability is not None:
        # A posterior margin m says the called relation is wrong with
        # probability 1 / (1 + e^m), which simulation found calibrated
        # when base qualities are (docs/design/posterior-margins.md); so "fragile" can mean "more
        # likely wrong than error_probability", a rule that does not depend
        # on how many boundaries there are. An MEC margin carries no such
        # meaning, so the rule is refused for it.
        if signal not in ("posterior", "shuffled"):
            raise ValueError("error_probability needs signal='posterior' or its control")
        if not 0.0 < error_probability < 0.5:
            raise ValueError("error_probability must lie in (0, 0.5)")
        threshold = math.log((1.0 - error_probability) / error_probability)
    if margin_source not in MARGIN_SOURCES:
        raise ValueError(
            f"margin_source must be one of {MARGIN_SOURCES}, got {margin_source!r}"
        )

    dp_kwargs = {} if max_active_reads is None else {"max_active_reads": max_active_reads}
    weights = read_weights(matrix, score)
    # Measured once. A full-matrix margin does not depend on the selection, so
    # the fragile set is the same in every round and the loop stops as soon as
    # the boosted weights stop moving the selection.
    shuffler = random.Random(shuffle_seed)

    def fragility(target: FragmentMatrix, dp: Optional[CoverageDP] = None):
        if signal == "posterior":
            return posterior_margins(target)
        if signal == "shuffled":
            margins = posterior_margins(target)
            values = [m.margin for m in margins if m.spanned]
            shuffler.shuffle(values)
            permuted = iter(values)
            return tuple(
                Margin(m.boundary, m.relation, next(permuted), True) if m.spanned else m
                for m in margins
            )
        return (dp or CoverageDP(target, **dp_kwargs)).margins()

    fixed_margins = fragility(matrix) if margin_source == "full" else None
    history: list[AdaptiveStep] = []
    best: Optional[AdaptiveStep] = None
    seen: set[tuple[int, ...]] = set()

    for iteration in range(iterations):
        selection = selector(matrix, coverage, weights)
        subset = matrix.subset(selection.selected)
        dp = CoverageDP(subset, **dp_kwargs)
        solution = dp.solve()
        full_cost = matrix.cost(solution.haplotype)
        margins = fragility(subset, dp) if fixed_margins is None else fixed_margins
        fragile = fragile_boundaries(margins, quantile=quantile, threshold=threshold)
        boosted = tuple(
            fragment.read_id
            for fragment in matrix
            if any(fragment.spans(boundary) for boundary in fragile)
        )
        step = AdaptiveStep(
            iteration=iteration,
            selection=selection,
            solution=solution,
            full_cost=full_cost,
            margins=margins,
            fragile=fragile,
            boosted=boosted,
            margin_source=margin_source,
        )
        history.append(step)
        if best is None or full_cost < best.full_cost:
            best = step
        if selection.selected in seen or not boosted:
            break
        seen.add(selection.selected)
        for read_id in boosted:
            weights[read_id] *= boost

    assert best is not None
    solution, full_cost, support = best.solution, best.full_cost, best.selection.selected
    if acceptance in ("stitched", "stitched_oriented") and len(history) > 1:
        stitched, choice = stitch_haplotypes(
            matrix, [step.solution.haplotype for step in history], window
        )
        stitched_cost = matrix.cost(stitched)
        if stitched_cost < best.full_cost:
            solution = WMECSolution(
                cost=stitched_cost,
                haplotype=stitched,
                assignment=matrix.best_assignment(stitched),
            )
            full_cost = stitched_cost
            support = tuple(sorted({
                read for r in set(choice) for read in history[r].selection.selected
            }))
    if acceptance == "stitched_oriented" and len(history) > 1:
        # Repair of OQ-15: the stitcher can only reorient at window junctions,
        # so reorient every piece between a junction or a block break of the
        # support, jointly and exactly, by every read.
        segment_of = _segments(matrix, support, window)
        oriented = orient_segments(matrix, solution.haplotype, segment_of,
                                   1 + max((s for s in segment_of if s is not None), default=-1))
        oriented_cost = matrix.cost(oriented.haplotype)
        if oriented_cost < full_cost:
            solution = WMECSolution(cost=oriented_cost, haplotype=oriented.haplotype,
                                    assignment=matrix.best_assignment(oriented.haplotype))
            full_cost = oriented_cost
    return AdaptiveResult(
        solution=solution,
        selection=best.selection,
        full_cost=full_cost,
        iterations=len(history),
        history=tuple(history),
        margin_source=margin_source,
        acceptance=acceptance,
        support=support,
    )


def _segments(matrix: FragmentMatrix, support, window: int) -> list[Optional[int]]:
    """Segment index per variant: cut at window junctions and at the support's breaks."""
    subset = matrix.subset(support)
    # A covered variant no support read calls is unphased in `blocks()`, but it is
    # inside its block's span: it must ride with that block, not break the segment.
    blocks = list(subset.blocks())
    coverage = subset.coverage()
    for j in range(len(blocks)):
        if blocks[j] is None and coverage[j] > 0 and j > 0:
            blocks[j] = blocks[j - 1]
    spanned = subset.spanned_boundaries()
    segment_of: list[Optional[int]] = []
    current = -1
    previous_block = None
    for j, block in enumerate(blocks):
        if block is None:
            segment_of.append(None)
            previous_block = None
            continue
        if (current < 0 or previous_block is None or j % window == 0
                or not spanned[j - 1] or block != previous_block):
            current += 1
        segment_of.append(current)
        previous_block = block
    return segment_of


def static_selection(
    matrix: FragmentMatrix,
    coverage: int,
    *,
    score: Callable = total_weight_score,
    selector: Callable[..., SelectionResult] = select_optimal,
    max_active_reads: Optional[int] = None,
) -> AdaptiveResult:
    """One round of the same pipeline, which is the comparison baseline."""
    return adaptive_selection(
        matrix,
        coverage,
        iterations=1,
        score=score,
        selector=selector,
        max_active_reads=max_active_reads,
    )
