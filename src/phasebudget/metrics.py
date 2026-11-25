"""Phasing accuracy and completeness metrics.

Definitions are written out rather than assumed, because switch error, flip
error, and the N50 family are all reported under more than one convention.
`docs/design/metrics.md` states the convention used here and where it differs
from the ones in common use.
"""

from __future__ import annotations

from typing import Optional, Sequence

from ._types import Allele, Phasing, WMECSolution
from .fragments import FragmentMatrix


def phasing_from_solution(
    matrix: FragmentMatrix,
    solution: WMECSolution,
    *,
    full_matrix: Optional[FragmentMatrix] = None,
) -> Phasing:
    """Turn a wMEC solution into a phase call with block structure.

    Blocks come from the matrix that was solved, since connectivity is a
    property of the selected reads. A variant no selected read covers is
    reported unphased even when the solver assigned it an allele. Passing
    `full_matrix` changes nothing about the blocks and is accepted only so
    callers can assert that the two matrices agree on variant count.
    """
    if full_matrix is not None and full_matrix.n_variants != matrix.n_variants:
        raise ValueError("the two matrices disagree on the number of variants")
    blocks = matrix.blocks()
    alleles: list[Allele] = [
        solution.haplotype[j] if blocks[j] is not None else None
        for j in range(matrix.n_variants)
    ]
    return Phasing(alleles=tuple(alleles), blocks=blocks)


def _block_positions(phasing: Phasing) -> dict[int, list[int]]:
    blocks: dict[int, list[int]] = {}
    for j, (allele, block) in enumerate(zip(phasing.alleles, phasing.blocks)):
        if allele is None or block is None:
            continue
        blocks.setdefault(block, []).append(j)
    return {block: members for block, members in blocks.items() if len(members) > 1}


def _compared_members(
    truth: Sequence[Optional[int]], phasing: Phasing
) -> list[list[int]]:
    """Per block, the phased variants the truth also phases, in order.

    A real truth set leaves some variants unphased, marked `None`. Those are
    skipped, and the next pair compared is the next two variants of the block
    that both phasings call, which is how `whatshap compare` counts. With a
    complete truth this is `_block_positions` unchanged.
    """
    if len(truth) != len(phasing):
        raise ValueError("truth and phasing must cover the same variants")
    compared = []
    for members in _block_positions(phasing).values():
        kept = [j for j in members if truth[j] is not None]
        if len(kept) > 1:
            compared.append(kept)
    return compared


def phased_fraction(phasing: Phasing) -> float:
    """Fraction of variants phased relative to at least one other variant.

    A variant alone in its block carries an allele but no phase information,
    so it does not count here. Counting it would let a method raise its phased
    fraction by emitting singleton blocks.
    """
    if not len(phasing):
        return 0.0
    phased = sum(len(members) for members in _block_positions(phasing).values())
    return phased / len(phasing)


def switch_positions(truth: Sequence[Optional[int]], phasing: Phasing) -> tuple[int, ...]:
    """Indices `i` where the phase relation to the previous phased variant flips.

    Positions are indices into the phasing, and only pairs inside one block are
    compared. A block phased in full agreement or in full opposition to the
    truth contributes nothing, which is the point of measuring switches rather
    than allele identity. Variants the truth leaves as `None` are skipped.
    """
    switches: list[int] = []
    for members in _compared_members(truth, phasing):
        for left, right in zip(members, members[1:]):
            predicted = phasing.alleles[left] ^ phasing.alleles[right]
            actual = truth[left] ^ truth[right]
            if predicted != actual:
                switches.append(right)
    return tuple(sorted(switches))


def comparable_pairs(
    phasing: Phasing, truth: Optional[Sequence[Optional[int]]] = None
) -> int:
    """Number of adjacent phased pairs inside a block, the switch denominator.

    With `truth`, only pairs of variants the truth also phases are counted.
    """
    if truth is None:
        return sum(len(members) - 1 for members in _block_positions(phasing).values())
    return sum(len(members) - 1 for members in _compared_members(truth, phasing))


def switch_error_rate(truth: Sequence[Optional[int]], phasing: Phasing) -> float:
    denominator = comparable_pairs(phasing, truth)
    if denominator == 0:
        return 0.0
    return len(switch_positions(truth, phasing)) / denominator


def flip_errors(truth: Sequence[int], phasing: Phasing) -> tuple[int, ...]:
    """Variants isolated by two consecutive switches, counted once.

    Two switches at consecutive comparable pairs put one variant on the wrong
    haplotype and return the phase afterwards. Reporting that as two switch
    errors overstates the damage, so it is reported here as one flip.
    """
    consumed = flips_and_partners(truth, phasing)
    return tuple(consumed[index] for index in range(0, len(consumed), 2))


def switch_errors_excluding_flips(truth: Sequence[int], phasing: Phasing) -> int:
    """Switch count after every flip has been removed as a pair of switches."""
    switches = switch_positions(truth, phasing)
    flips = set(flips_and_partners(truth, phasing))
    return sum(1 for position in switches if position not in flips)


def flips_and_partners(truth: Sequence[int], phasing: Phasing) -> tuple[int, ...]:
    """The switch positions consumed by flips, both members of each pair."""
    switches = set(switch_positions(truth, phasing))
    consumed: list[int] = []
    for members in _compared_members(truth, phasing):
        index = 1
        while index < len(members) - 1:
            here, following = members[index], members[index + 1]
            if here in switches and following in switches:
                consumed.extend((here, following))
                index += 2
            else:
                index += 1
    return tuple(sorted(consumed))


def n50(lengths: Sequence[float]) -> float:
    """Largest `L` such that blocks of length at least `L` cover half the total."""
    positive = sorted((value for value in lengths if value > 0), reverse=True)
    total = sum(positive)
    if total == 0:
        return 0.0
    running = 0.0
    for value in positive:
        running += value
        if running * 2 >= total:
            return value
    return positive[-1]


def block_lengths(
    phasing: Phasing, positions: Optional[Sequence[int]] = None
) -> tuple[float, ...]:
    """Span of every block, in base pairs if `positions` is given else variants."""
    lengths = []
    for members in _block_positions(phasing).values():
        if positions is None:
            lengths.append(float(len(members)))
        else:
            lengths.append(float(positions[members[-1]] - positions[members[0]] + 1))
    return tuple(lengths)


def block_n50(phasing: Phasing, positions: Optional[Sequence[int]] = None) -> float:
    return n50(block_lengths(phasing, positions))


def an50_lengths(
    phasing: Phasing, positions: Optional[Sequence[int]] = None
) -> tuple[float, ...]:
    """Block spans scaled by the fraction of the span that is phased.

    A block that skips half the variants inside its own span is credited with
    half its span, which stops a sparse block from claiming continuity it does
    not have. Lists from several chromosomes can be joined and passed to `n50`.
    """
    adjusted = []
    for members in _block_positions(phasing).values():
        first, last = members[0], members[-1]
        if positions is None:
            span = float(last - first + 1)
        else:
            span = float(positions[last] - positions[first] + 1)
        covered = last - first + 1
        adjusted.append(span * len(members) / covered)
    return tuple(adjusted)


def an50(phasing: Phasing, positions: Optional[Sequence[int]] = None) -> float:
    """N50 over `an50_lengths`."""
    return n50(an50_lengths(phasing, positions))


def qan50_lengths(
    truth: Sequence[int],
    phasing: Phasing,
    positions: Optional[Sequence[int]] = None,
) -> tuple[float, ...]:
    """`an50_lengths` after every block has been cut at each switch error."""
    switches = set(switch_positions(truth, phasing))
    alleles = list(phasing.alleles)
    blocks: list[Optional[int]] = list(phasing.blocks)
    next_id = max((block for block in blocks if block is not None), default=-1) + 1
    for members in _block_positions(phasing).values():
        current = blocks[members[0]]
        for position in members[1:]:
            if position in switches:
                current = next_id
                next_id += 1
            blocks[position] = current
    return an50_lengths(Phasing(tuple(alleles), tuple(blocks)), positions)


def qan50(
    truth: Sequence[int],
    phasing: Phasing,
    positions: Optional[Sequence[int]] = None,
) -> float:
    """AN50 after every block has been cut at each switch error.

    This is the quantity that punishes a long wrong block, and the reason a
    method cannot buy N50 by extending blocks through low-confidence regions.
    """
    return n50(qan50_lengths(truth, phasing, positions))


def summarize(
    truth: Sequence[int],
    phasing: Phasing,
    positions: Optional[Sequence[int]] = None,
) -> dict[str, float]:
    """Every metric in one dictionary, for benchmark tables."""
    return {
        "phased_fraction": phased_fraction(phasing),
        "comparable_pairs": float(comparable_pairs(phasing, truth)),
        "switch_errors": float(len(switch_positions(truth, phasing))),
        "switch_error_rate": switch_error_rate(truth, phasing),
        "flip_errors": float(len(flip_errors(truth, phasing))),
        "switch_errors_excluding_flips": float(
            switch_errors_excluding_flips(truth, phasing)
        ),
        "blocks": float(len(_block_positions(phasing))),
        "block_n50": block_n50(phasing, positions),
        "an50": an50(phasing, positions),
        "qan50": qan50(truth, phasing, positions),
    }
