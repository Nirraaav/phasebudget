"""Run WhatsHap's own core on phasebudget's fragment matrices.

WhatsHap is not a dependency of phasebudget; this module is imported only by
`hg002_pipelines.py` and `tests/test_whatshap_port.py`, which skip without it. It converts a `FragmentMatrix` into a WhatsHap `ReadSet`,
runs WhatsHap's exact dynamic program (`PedigreeDPTable`, one individual,
every site heterozygous, genotypes trusted) and its deployed read selection
(`readselection`), and maps the results back.

WhatsHap takes integer qualities, so weights are scaled and rounded by
`integer_weights`; comparisons are made on the matrix with those integer
weights, so both sides solve the same instance exactly. Genomic positions
default to 10, 20, 30, ...; a simulated instance's own `positions` can be
passed instead. Only the ordering of tied reads in the selection depends on
them.
"""

from __future__ import annotations

from typing import Optional, Sequence

from phasebudget._types import Fragment
from phasebudget.fragments import FragmentMatrix

SAMPLE = "sample"


def integer_weights(matrix: FragmentMatrix, scale: float = 10.0) -> FragmentMatrix:
    """The same matrix with every weight replaced by `round(scale * weight)`, at least 1."""
    fragments = []
    for fragment in matrix:
        weights = tuple(
            0.0 if allele is None else float(max(1, round(scale * weight)))
            for allele, weight in zip(fragment.alleles, fragment.weights)
        )
        fragments.append(Fragment(fragment.read_id, fragment.start, fragment.alleles, weights))
    return FragmentMatrix(fragments, matrix.n_variants)


def phase_informative(matrix: FragmentMatrix) -> FragmentMatrix:
    """The reads calling at least two variants, which is what WhatsHap selects among.

    `whatshap phase` drops the others before read selection, and
    `readselection` refuses them.
    """
    from phasebudget.selection import is_phase_informative

    return matrix.subset([f.read_id for f in matrix if is_phase_informative(f)])


def positions_of(matrix: FragmentMatrix, positions: Optional[Sequence[int]] = None) -> list[int]:
    """Genomic positions for WhatsHap: `positions` if given, else 10, 20, 30, ..."""
    if positions is not None:
        return list(positions)
    return [10 * (j + 1) for j in range(matrix.n_variants)]


def to_readset(matrix: FragmentMatrix, positions: Optional[Sequence[int]] = None):
    """A WhatsHap ReadSet with one Read per fragment, in matrix order."""
    from whatshap.core import NumericSampleIds, Read, ReadSet

    ids = NumericSampleIds()
    sample_id = ids[SAMPLE]
    positions = positions_of(matrix, positions)
    readset = ReadSet()
    for index, fragment in enumerate(matrix):
        read = Read(f"r{index}", 60, 0, sample_id)
        for offset, allele in enumerate(fragment.alleles):
            if allele is None:
                continue
            read.add_variant(
                position=positions[fragment.start + offset],
                allele=int(allele),
                quality=int(round(fragment.weights[offset])),
            )
        read.sort()
        readset.add(read)
    return readset, ids


def whatshap_cost(matrix: FragmentMatrix, positions: Optional[Sequence[int]] = None) -> int:
    """WhatsHap's optimal cost for one heterozygous individual with trusted genotypes."""
    from whatshap.core import Genotype, Pedigree, PedigreeDPTable

    readset, ids = to_readset(matrix, positions)
    readset.sort()
    pedigree = Pedigree(ids)
    pedigree.add_individual(SAMPLE, [Genotype([0, 1]) for _ in range(matrix.n_variants)])
    table = PedigreeDPTable(
        readset, [1] * matrix.n_variants, pedigree, False, positions_of(matrix, positions)
    )
    return table.get_optimal_cost()


def whatshap_selection(
    matrix: FragmentMatrix, coverage: int, positions: Optional[Sequence[int]] = None
) -> tuple[int, ...]:
    """Read ids of `matrix` that WhatsHap's `readselection` keeps at `coverage`."""
    from whatshap.readselect import readselection

    readset, _ = to_readset(matrix, positions)
    names = {f"r{index}": fragment.read_id for index, fragment in enumerate(matrix)}
    selected = readselection(readset, coverage)
    return tuple(sorted(names[readset[i].name] for i in selected))
