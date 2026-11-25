"""Synthetic instances with known ground truth.

Two generators. `simulate` takes read span in variants, so heterozygosity does
not exist as a quantity in it. `simulate_genomic` places variants along a
genome at a given density and gives reads a physical length, which makes the
expected variants per read a consequence of `heterozygosity * read_length`
rather than a parameter. `docs/design/genomic-simulator.md` states why the
first cannot express a density sweep and the second can.

Coverage, sequencing error, and base quality vary independently in both. Every
generator takes a seed and is reproducible from it alone.
"""

from __future__ import annotations

import bisect
import math
import random
from dataclasses import dataclass
from typing import Optional, Sequence

from ._types import Allele, Fragment
from .fragments import FragmentMatrix


@dataclass(frozen=True)
class SimulatedInstance:
    """A fragment matrix with the haplotype and read origins that produced it."""

    matrix: FragmentMatrix
    haplotype: tuple[int, ...]
    origin: tuple[int, ...]
    positions: tuple[int, ...]
    error_rates: tuple[float, ...]
    homozygous: tuple[int, ...] = ()
    systematic: tuple[int, ...] = ()

    @property
    def n_variants(self) -> int:
        return self.matrix.n_variants

    def complement(self) -> tuple[int, ...]:
        return tuple(1 - allele for allele in self.haplotype)

    def second_haplotype(self) -> tuple[int, ...]:
        """The other chromosome copy, which is the complement away from `homozygous`.

        A site in `homozygous` was called heterozygous and is not, which is the
        genotyping error OQ-7 is about. There the two copies agree, so the site
        carries an allele and no phase.
        """
        miscalled = set(self.homozygous)
        return tuple(
            allele if index in miscalled else 1 - allele
            for index, allele in enumerate(self.haplotype)
        )

    def heterozygous_sites(self) -> tuple[int, ...]:
        """Indices where the two copies really do differ."""
        miscalled = set(self.homozygous)
        return tuple(
            index for index in range(self.n_variants) if index not in miscalled
        )


def log_odds_weight(error_rate: float) -> float:
    """Confidence weight of one observed allele, in log-odds.

    An allele called with error probability `e` costs `log((1-e)/e)` to
    correct. The weight is the standard wMEC one: it is zero when the call is a
    coin flip and grows without bound as the call becomes certain.
    """
    if not 0.0 < error_rate < 1.0:
        raise ValueError(f"error_rate must lie strictly in (0, 1), got {error_rate}")
    return math.log((1.0 - error_rate) / error_rate)


def phred_to_error(quality: float) -> float:
    return 10.0 ** (-quality / 10.0)


COIN_FLIP_PHRED = 10.0 * math.log10(2.0)
"""Phred score of a call with error probability 0.5, which weighs zero."""


def jittered_error(
    rng: random.Random, error_rate: float, spread: float
) -> float:
    """An error probability drawn from a phred score jittered around `error_rate`.

    The score is floored at the coin-flip value, so the drawn probability never
    exceeds 0.5 and `log_odds_weight` never returns a negative number. A call
    worse than a coin flip is a mislabelled allele rather than a low-quality
    one, and the weight model has no place to put it: a negative weight would
    pay the objective to corrupt that allele, and `check_fragment` rejects it.
    """
    if spread <= 0.0:
        return error_rate
    phred = -10.0 * math.log10(error_rate)
    return phred_to_error(max(COIN_FLIP_PHRED, rng.gauss(phred, spread)))


def simulate(
    n_variants: int,
    *,
    coverage: float = 10.0,
    read_span: int = 6,
    error_rate: float = 0.02,
    seed: int = 0,
    span_jitter: int = 0,
    gap_rate: float = 0.0,
    quality_spread: float = 0.0,
    weight_model: str = "log_odds",
    positions_spacing: int = 1000,
) -> SimulatedInstance:
    """Draw reads from two complementary haplotypes.

    `coverage` is the expected number of reads covering a variant, `read_span`
    the number of variants a read touches, `gap_rate` the chance that an
    interior variant inside a span goes uncalled, and `quality_spread` the
    width in phred units of the per-base quality jitter around `error_rate`.
    With `weight_model="uniform"` every observation weighs 1.0, which turns the
    objective into unweighted MEC.
    """
    if n_variants <= 0:
        raise ValueError("n_variants must be positive")
    if read_span <= 0:
        raise ValueError("read_span must be positive")
    if coverage < 0:
        raise ValueError("coverage must be non-negative")
    if weight_model not in ("log_odds", "uniform"):
        raise ValueError(f"unknown weight model {weight_model!r}")
    rng = random.Random(seed)
    haplotype = tuple(rng.randint(0, 1) for _ in range(n_variants))
    n_reads = max(0, round(coverage * n_variants / read_span))

    fragments: list[Fragment] = []
    origins: list[int] = []
    rates: list[float] = []
    next_id = 0
    for _ in range(n_reads):
        span = read_span
        if span_jitter:
            span = max(1, read_span + rng.randint(-span_jitter, span_jitter))
        start = rng.randint(0, max(0, n_variants - 1))
        stop = min(n_variants, start + span)
        side = rng.randint(0, 1)
        alleles: list[Allele] = []
        weights: list[float] = []
        for position in range(start, stop):
            interior = position not in (start, stop - 1)
            if interior and rng.random() < gap_rate:
                alleles.append(None)
                weights.append(0.0)
                continue
            local = jittered_error(rng, error_rate, quality_spread)
            true_allele = haplotype[position] ^ side
            observed = 1 - true_allele if rng.random() < local else true_allele
            alleles.append(observed)
            weights.append(1.0 if weight_model == "uniform" else log_odds_weight(local))
            rates.append(local)
        while alleles and alleles[0] is None:
            alleles.pop(0)
            weights.pop(0)
            start += 1
        while alleles and alleles[-1] is None:
            alleles.pop()
            weights.pop()
        if not alleles:
            continue
        fragments.append(
            Fragment(
                read_id=next_id,
                start=start,
                alleles=tuple(alleles),
                weights=tuple(weights),
            )
        )
        origins.append(side)
        next_id += 1

    positions = tuple(
        index * positions_spacing + rng.randint(0, positions_spacing - 1)
        for index in range(n_variants)
    )
    return SimulatedInstance(
        matrix=FragmentMatrix(fragments, n_variants),
        haplotype=haplotype,
        origin=tuple(origins),
        positions=tuple(sorted(positions)),
        error_rates=tuple(rates),
    )


def simulate_connected(
    n_variants: int,
    *,
    coverage: float = 10.0,
    read_span: int = 6,
    error_rate: float = 0.02,
    seed: int = 0,
    attempts: int = 200,
    **kwargs,
) -> SimulatedInstance:
    """Retry `simulate` until every variant lies in one phase block.

    A disconnected instance is a legitimate case and the metrics handle it, but
    a benchmark that varies coverage should not silently vary connectivity too.
    Raises `RuntimeError` when the parameters cannot produce one block.
    """
    for attempt in range(attempts):
        instance = simulate(
            n_variants,
            coverage=coverage,
            read_span=read_span,
            error_rate=error_rate,
            seed=seed + attempt,
            **kwargs,
        )
        blocks = instance.matrix.blocks()
        # A variant no read calls is unphased (None) but does not break the block.
        if all(instance.matrix.coverage()) and all(block in (0, None) for block in blocks):
            return instance
    raise RuntimeError(
        f"no single-block instance in {attempts} attempts at coverage "
        f"{coverage} and span {read_span}"
    )


def fragments_from_assignment(
    haplotype: Sequence[int],
    spans: Sequence[tuple[int, int]],
    sides: Sequence[int],
    *,
    weights: Optional[Sequence[float]] = None,
) -> FragmentMatrix:
    """Build an error-free matrix directly, for tests that need a fixed shape."""
    if len(spans) != len(sides):
        raise ValueError("one side is needed per span")
    fragments = []
    for index, ((start, stop), side) in enumerate(zip(spans, sides)):
        alleles = tuple(haplotype[j] ^ side for j in range(start, stop))
        row_weights = tuple(
            1.0 if weights is None else float(weights[index]) for _ in alleles
        )
        fragments.append(Fragment(index, start, alleles, row_weights))
    return FragmentMatrix(fragments, len(haplotype))


def geometric_gap(rng: random.Random, rate: float) -> int:
    """A waiting time in `{1, 2, ...}` with mean `1 / rate`.

    Inverse transform on the geometric distribution, which is the gap between
    successive successes of a Bernoulli process at probability `rate`. This is
    the distance between heterozygous sites under the standard independent-sites
    model, and `random.Random` has no geometric variate of its own.
    """
    if not 0.0 < rate <= 1.0:
        raise ValueError(f"rate must lie in (0, 1], got {rate}")
    if rate == 1.0:
        return 1
    return 1 + int(math.log(1.0 - rng.random()) / math.log(1.0 - rate))


def simulate_genomic(
    n_variants: int,
    *,
    coverage: float = 10.0,
    read_length: int = 10000,
    heterozygosity: float = 1e-3,
    error_rate: float = 0.02,
    seed: int = 0,
    gap_rate: float = 0.0,
    quality_spread: float = 0.0,
    weight_model: str = "log_odds",
    first_position: int = 1,
    systematic_rate: float = 0.0,
    systematic_strength: float = 0.5,
    burst_rate: float = 0.0,
    burst_length: int = 5,
    burst_error: float = 0.35,
) -> SimulatedInstance:
    """Draw reads of a fixed physical length over variants at a given density.

    `heterozygosity` is the per-base probability that a site is a heterozygous
    variant, so successive variant positions are `geometric_gap` apart and the
    expected number of variants a read touches is
    `heterozygosity * read_length`. `coverage` is the expected number of reads
    over a base.

    Read start points run from `read_length - 1` bases before the first variant
    to the last variant, so every variant has exactly `read_length` start
    positions that would cover it and expected per-variant coverage is flat.
    A read covering no variant is dropped; one covering a single variant is
    kept, because it occupies coverage without carrying a phase link.

    `systematic_rate`, `systematic_strength`, `burst_rate`, `burst_length`
    and `burst_error` add errors the weights do not report, as
    `docs/design/genomic-simulator.md` defines; at their defaults no extra
    random numbers are drawn and earlier instances are reproduced exactly.
    """
    if n_variants <= 0:
        raise ValueError("n_variants must be positive")
    if read_length <= 0:
        raise ValueError("read_length must be positive")
    if coverage < 0:
        raise ValueError("coverage must be non-negative")
    if not 0.0 < heterozygosity <= 1.0:
        raise ValueError(
            f"heterozygosity must lie in (0, 1], got {heterozygosity}"
        )
    if weight_model not in ("log_odds", "uniform"):
        raise ValueError(f"unknown weight model {weight_model!r}")
    if not 0.0 <= systematic_rate <= 1.0 or not 0.0 <= systematic_strength <= 1.0:
        raise ValueError("systematic_rate and systematic_strength must lie in [0, 1]")
    if not 0.0 <= burst_rate <= 1.0 or not 0.0 <= burst_error <= 1.0 or burst_length < 1:
        raise ValueError("burst_rate and burst_error must lie in [0, 1], burst_length >= 1")

    rng = random.Random(seed)
    positions = [first_position]
    for _ in range(n_variants - 1):
        positions.append(positions[-1] + geometric_gap(rng, heterozygosity))
    haplotype = tuple(rng.randint(0, 1) for _ in range(n_variants))
    biased: dict[int, int] = {}
    if systematic_rate > 0.0:
        for index in range(n_variants):
            if rng.random() < systematic_rate:
                biased[index] = rng.randint(0, 1)

    lowest = positions[0] - read_length + 1
    highest = positions[-1]
    window = highest - lowest + 1
    n_reads = max(0, round(coverage * window / read_length))

    fragments: list[Fragment] = []
    origins: list[int] = []
    rates: list[float] = []
    next_id = 0
    for _ in range(n_reads):
        begin = rng.randint(lowest, highest)
        first = bisect.bisect_left(positions, begin)
        stop = bisect.bisect_left(positions, begin + read_length)
        if first == stop:
            continue
        side = rng.randint(0, 1)
        burst = range(0)
        if burst_rate > 0.0 and rng.random() < burst_rate:
            begin_burst = rng.randint(first, max(first, stop - burst_length))
            burst = range(begin_burst, begin_burst + burst_length)
        alleles: list[Allele] = []
        weights: list[float] = []
        start = first
        for index in range(first, stop):
            interior = index not in (first, stop - 1)
            if interior and rng.random() < gap_rate:
                alleles.append(None)
                weights.append(0.0)
                continue
            local = jittered_error(rng, error_rate, quality_spread)
            true_allele = haplotype[index] ^ side
            actual = burst_error if index in burst else local
            observed = 1 - true_allele if rng.random() < actual else true_allele
            if index in biased and rng.random() < systematic_strength:
                observed = biased[index]
            alleles.append(observed)
            weights.append(1.0 if weight_model == "uniform" else log_odds_weight(local))
            rates.append(local)
        while alleles and alleles[0] is None:
            alleles.pop(0)
            weights.pop(0)
            start += 1
        while alleles and alleles[-1] is None:
            alleles.pop()
            weights.pop()
        if not alleles:
            continue
        fragments.append(
            Fragment(
                read_id=next_id,
                start=start,
                alleles=tuple(alleles),
                weights=tuple(weights),
            )
        )
        origins.append(side)
        next_id += 1

    return SimulatedInstance(
        matrix=FragmentMatrix(fragments, n_variants),
        haplotype=haplotype,
        origin=tuple(origins),
        positions=tuple(positions),
        error_rates=tuple(rates),
        systematic=tuple(sorted(biased)),
    )


def simulate_genomic_connected(
    n_variants: int,
    *,
    seed: int = 0,
    attempts: int = 200,
    **kwargs,
) -> SimulatedInstance:
    """Retry `simulate_genomic` until every variant lies in one phase block.

    The same guard as `simulate_connected` and for the same reason: a sweep
    over heterozygosity changes connectivity fastest of all the parameters,
    since lowering the density is exactly what disconnects the blocks. A sweep
    that wants to separate accuracy from connectivity has to hold one fixed.
    Raises `RuntimeError` when the parameters cannot produce one block.
    """
    for attempt in range(attempts):
        instance = simulate_genomic(n_variants, seed=seed + attempt, **kwargs)
        blocks = instance.matrix.blocks()
        if all(instance.matrix.coverage()) and all(block in (0, None) for block in blocks):
            return instance
    raise RuntimeError(
        f"no single-block instance in {attempts} attempts at "
        f"heterozygosity {kwargs.get('heterozygosity', 1e-3)} and read length "
        f"{kwargs.get('read_length', 10000)}"
    )


def simulate_with_genotype_errors(
    n_variants: int,
    *,
    coverage: float = 8.0,
    read_span: int = 5,
    error_rate: float = 0.02,
    homozygous_rate: float = 0.0,
    seed: int = 0,
    weight_model: str = "log_odds",
    positions_spacing: int = 1000,
) -> SimulatedInstance:
    """Reads from two copies that agree at a fraction of the called sites.

    Every column is presented to the solver as a heterozygous variant, but a
    site drawn into `homozygous` really carries the same allele on both copies.
    That is a genotyping error, and it is the case the complementary model in
    `phasebudget.dp` cannot represent: it will insist the two copies differ
    there and pay for the correction somewhere. `docs/design/genotype-errors.md`
    states the comparison this generator exists for.

    Kept deliberately small. The general-model oracle in
    `exact.solve_bipartition_bruteforce` enumerates read bipartitions, so the
    comparison is only available at a handful of reads.
    """
    if n_variants <= 0:
        raise ValueError("n_variants must be positive")
    if read_span <= 0:
        raise ValueError("read_span must be positive")
    if coverage < 0:
        raise ValueError("coverage must be non-negative")
    if not 0.0 <= homozygous_rate <= 1.0:
        raise ValueError(
            f"homozygous_rate must lie in [0, 1], got {homozygous_rate}"
        )
    if weight_model not in ("log_odds", "uniform"):
        raise ValueError(f"unknown weight model {weight_model!r}")
    if not 0.0 <= error_rate < 1.0:
        raise ValueError(f"error_rate must lie in [0, 1), got {error_rate}")
    if error_rate == 0.0 and weight_model == "log_odds":
        # log((1-e)/e) diverges at e = 0. An error-free instance has to say
        # what a correction costs, and uniform weights are that answer.
        raise ValueError(
            "error_rate 0.0 has no finite log-odds weight; pass "
            'weight_model="uniform" for an error-free instance'
        )

    rng = random.Random(seed)
    first = tuple(rng.randint(0, 1) for _ in range(n_variants))
    homozygous = tuple(
        index for index in range(n_variants) if rng.random() < homozygous_rate
    )
    miscalled = set(homozygous)
    second = tuple(
        allele if index in miscalled else 1 - allele
        for index, allele in enumerate(first)
    )
    copies = (first, second)
    n_reads = max(0, round(coverage * n_variants / read_span))

    fragments: list[Fragment] = []
    origins: list[int] = []
    rates: list[float] = []
    next_id = 0
    for _ in range(n_reads):
        start = rng.randint(0, max(0, n_variants - 1))
        stop = min(n_variants, start + read_span)
        side = rng.randint(0, 1)
        alleles: list[Allele] = []
        weights: list[float] = []
        for position in range(start, stop):
            true_allele = copies[side][position]
            observed = (
                1 - true_allele if rng.random() < error_rate else true_allele
            )
            alleles.append(observed)
            weights.append(
                1.0 if weight_model == "uniform" else log_odds_weight(error_rate)
            )
            rates.append(error_rate)
        if not alleles:
            continue
        fragments.append(
            Fragment(
                read_id=next_id,
                start=start,
                alleles=tuple(alleles),
                weights=tuple(weights),
            )
        )
        origins.append(side)
        next_id += 1

    positions = tuple(index * positions_spacing for index in range(n_variants))
    return SimulatedInstance(
        matrix=FragmentMatrix(fragments, n_variants),
        haplotype=first,
        origin=tuple(origins),
        positions=positions,
        error_rates=tuple(rates),
        homozygous=homozygous,
    )
