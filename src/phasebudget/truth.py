"""A truth phasing for a trio child, by transmission from the parents' calls.

`docs/design/truth-sets.md` states the rules; this module follows them. The
result uses no reads and no phasing tool, which is why it exists: the phased
GIAB files for HG002 were made partly with WhatsHap or are unevaluated.
"""

from __future__ import annotations

import bisect
import gzip
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Union

PathLike = Union[str, Path]


@dataclass(frozen=True)
class Call:
    """One VCF record for the first sample, reduced to what phasing needs."""

    position: int
    ref: str
    alt: tuple[str, ...]
    genotype: tuple[Optional[int], ...]
    phased: bool
    phase_set: Optional[str]

    @property
    def is_biallelic_snv(self) -> bool:
        return len(self.ref) == 1 and len(self.alt) == 1 and len(self.alt[0]) == 1

    def bases(self) -> Optional[tuple[str, ...]]:
        """The genotype as bases, or None when an allele is missing."""
        alleles = (self.ref,) + self.alt
        if any(allele is None or allele >= len(alleles) for allele in self.genotype):
            return None
        return tuple(alleles[allele] for allele in self.genotype)  # type: ignore[index]


def _open(path: PathLike):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


def _parse_genotype(text: str) -> tuple[tuple[Optional[int], ...], bool]:
    phased = "|" in text
    parts = text.replace("|", "/").split("/")
    return tuple(None if part == "." else int(part) for part in parts), phased


def read_calls(path: PathLike, chromosome: str) -> dict[int, Call]:
    """Calls of the first sample on one chromosome, keyed by 1-based position."""
    calls: dict[int, Call] = {}
    with _open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] != chromosome or len(fields) < 10:
                continue
            keys = fields[8].split(":")
            values = fields[9].split(":")
            sample = dict(zip(keys, values))
            genotype, phased = _parse_genotype(sample.get("GT", "."))
            phase_set = sample.get("PS")
            position = int(fields[1])
            calls[position] = Call(
                position=position,
                ref=fields[3],
                alt=tuple(fields[4].split(",")),
                genotype=genotype,
                phased=phased,
                phase_set=None if phase_set in (None, ".") else phase_set,
            )
    return calls


class Regions:
    """BED intervals on one chromosome; BED is 0-based and half-open."""

    def __init__(self, intervals: Iterable[tuple[int, int]]) -> None:
        merged: list[list[int]] = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        self._starts = [start for start, _ in merged]
        self._ends = [end for _, end in merged]

    @classmethod
    def from_bed(cls, path: PathLike, chromosome: str) -> "Regions":
        intervals = []
        with _open(path) as handle:
            for line in handle:
                fields = line.split()
                if len(fields) >= 3 and fields[0] == chromosome:
                    intervals.append((int(fields[1]), int(fields[2])))
        return cls(intervals)

    def contains(self, position: int) -> bool:
        """True when the 1-based `position` lies inside an interval."""
        index = bisect.bisect_right(self._starts, position - 1) - 1
        return index >= 0 and position - 1 < self._ends[index]


def _parent_bases(
    call: Optional[Call], regions: Regions, position: int, ref: str
) -> Union[tuple[str, ...], str]:
    """The parent's two bases, or a reason the site cannot use this parent."""
    if call is None:
        return (ref, ref) if regions.contains(position) else "parent outside benchmark"
    if call.ref != ref:
        return "parent representation differs"
    bases = call.bases()
    if bases is None or len(bases) != 2:
        return "parent genotype missing"
    return bases


def phase_by_transmission(
    child: dict[int, Call],
    father: dict[int, Call],
    mother: dict[int, Call],
    father_regions: Regions,
    mother_regions: Regions,
) -> tuple[dict[int, tuple[int, int]], Counter]:
    """`(phased, reasons)`: paternal and maternal allele index per child site.

    Only biallelic SNVs with genotype `0/1` in the child are considered. Each
    one lands in `phased` or adds one to a reason in `reasons`, so the counts
    account for every candidate site.
    """
    phased: dict[int, tuple[int, int]] = {}
    reasons: Counter = Counter()
    for position, call in sorted(child.items()):
        if not call.is_biallelic_snv or sorted(call.genotype, key=lambda a: -1 if a is None else a) != [0, 1]:
            continue
        reasons["candidate"] += 1
        ref, alt = call.ref, call.alt[0]
        paternal = _parent_bases(father.get(position), father_regions, position, ref)
        maternal = _parent_bases(mother.get(position), mother_regions, position, ref)
        if isinstance(paternal, str):
            reasons[paternal] += 1
            continue
        if isinstance(maternal, str):
            reasons[maternal] += 1
            continue
        if any(base not in (ref, alt) for base in paternal + maternal):
            reasons["inconsistent"] += 1
            continue
        options = [
            (p, m) for p, m in ((ref, alt), (alt, ref)) if p in paternal and m in maternal
        ]
        if len(options) == 2:
            reasons["ambiguous, both parents heterozygous"] += 1
        elif not options:
            reasons["inconsistent"] += 1
        else:
            p, m = options[0]
            phased[position] = (0 if p == ref else 1, 0 if m == ref else 1)
            reasons["phased"] += 1
    return phased, reasons


@dataclass(frozen=True)
class Comparison:
    """Per-site agreement of a phased file with a transmission truth."""

    compared: int
    discordant: int
    not_oriented: int
    unphased_or_absent: int

    @property
    def discordance(self) -> float:
        return self.discordant / self.compared if self.compared else float("nan")


def compare_phasing(
    truth: dict[int, tuple[int, int]],
    other: dict[int, Call],
    oriented_sets: Iterable[str],
) -> Comparison:
    """Compare paternal alleles at the sites both phase, in paternal|maternal sets."""
    oriented = set(oriented_sets)
    compared = discordant = not_oriented = missing = 0
    for position, (paternal, _) in truth.items():
        call = other.get(position)
        if (
            call is None
            or not call.phased
            or not call.is_biallelic_snv
            or tuple(call.genotype) not in ((0, 1), (1, 0))
        ):
            missing += 1
            continue
        if call.phase_set not in oriented:
            not_oriented += 1
            continue
        compared += 1
        discordant += call.genotype[0] != paternal
    return Comparison(compared, discordant, not_oriented, missing)


def write_vcf(
    path: PathLike, chromosome: str, child: dict[int, Call], phased: dict[int, tuple[int, int]], sample: str
) -> None:
    """A minimal phased VCF of the transmission truth, paternal|maternal."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype, paternal|maternal">\n')
        handle.write('##FORMAT=<ID=PS,Number=1,Type=String,Description="Phase set; PATMAT is parent-oriented">\n')
        handle.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        for position in sorted(phased):
            call = child[position]
            paternal, maternal = phased[position]
            handle.write(
                f"{chromosome}\t{position}\t.\t{call.ref}\t{call.alt[0]}\t.\tPASS\t.\tGT:PS\t"
                f"{paternal}|{maternal}:PATMAT\n"
            )
