"""Phasing by transmission, one site per rule of docs/design/truth-sets.md."""

import pytest

from phasebudget.truth import Regions, compare_phasing, phase_by_transmission, read_calls, write_vcf

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"


def _vcf(tmp_path, name, rows):
    path = tmp_path / name
    body = "".join(
        f"chr20\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT:PS\t{gt}:{ps}\n"
        for pos, ref, alt, gt, ps in rows
    )
    path.write_text(HEADER + body)
    return path


@pytest.fixture
def trio(tmp_path):
    child = _vcf(tmp_path, "child.vcf", [
        (100, "A", "G", "0/1", "."),    # father hom ref by BED, mother 1/1: 0|1
        (200, "C", "T", "0/1", "."),    # father het, mother hom ref by BED: 1|0
        (300, "G", "A", "0/1", "."),    # both parents het: ambiguous
        (400, "T", "C", "0/1", "."),    # both parents 1/1: inconsistent
        (500, "A", "C", "0/1", "."),    # father outside his BED, no call: unknown
        (600, "C", "G", "0/1", "."),    # father call has another REF: skipped
        (700, "A", "C,G", "1/2", "."),  # multiallelic: not a candidate
        (800, "A", "G", "1/1", "."),    # homozygous: not a candidate
        (900, "A", "AT", "0/1", "."),   # indel: not a candidate
        (1000, "G", "T", "0/1", "."),   # mother carries a third base: inconsistent
    ])
    father = _vcf(tmp_path, "father.vcf", [
        (200, "C", "T", "0/1", "."),
        (300, "G", "A", "0/1", "."),
        (400, "T", "C", "1/1", "."),
        (600, "CA", "C", "0/1", "."),
    ])
    mother = _vcf(tmp_path, "mother.vcf", [
        (100, "A", "G", "1/1", "."),
        (300, "G", "A", "0/1", "."),
        (400, "T", "C", "1/1", "."),
        (1000, "G", "C", "0/1", "."),
    ])
    return (
        read_calls(child, "chr20"),
        read_calls(father, "chr20"),
        read_calls(mother, "chr20"),
        Regions([(0, 450), (550, 2000)]),
        Regions([(0, 2000)]),
    )


def test_each_rule_lands_where_the_design_note_says(trio):
    phased, reasons = phase_by_transmission(*trio)
    assert phased == {100: (0, 1), 200: (1, 0)}
    assert reasons["candidate"] == 7
    assert reasons["phased"] == 2
    assert reasons["ambiguous, both parents heterozygous"] == 1
    assert reasons["inconsistent"] == 2
    assert reasons["parent outside benchmark"] == 1
    assert reasons["parent representation differs"] == 1
    accounted = sum(count for reason, count in reasons.items() if reason != "candidate")
    assert accounted == reasons["candidate"]


def test_bed_intervals_are_zero_based_and_half_open():
    regions = Regions([(10, 20), (15, 30), (40, 41)])
    assert not regions.contains(10)
    assert regions.contains(11)
    assert regions.contains(30)
    assert not regions.contains(31)
    assert regions.contains(41)
    assert not regions.contains(42)


def test_comparison_counts_only_parent_oriented_sets(tmp_path):
    truth = {100: (0, 1), 200: (1, 0), 300: (0, 1), 400: (1, 0)}
    other = read_calls(_vcf(tmp_path, "other.vcf", [
        (100, "A", "G", "0|1", "PATMAT"),  # agrees
        (200, "C", "T", "0|1", "PATMAT"),  # disagrees
        (300, "G", "A", "1|0", "12345"),   # local phase set: not oriented
        (400, "T", "C", "0/1", "."),       # unphased
    ]), "chr20")
    result = compare_phasing(truth, other, {"PATMAT"})
    assert (result.compared, result.discordant, result.not_oriented, result.unphased_or_absent) == (2, 1, 1, 1)
    assert result.discordance == 0.5


def test_written_vcf_reads_back_as_the_same_phasing(trio, tmp_path):
    child = trio[0]
    phased = {100: (0, 1), 200: (1, 0)}
    path = tmp_path / "truth.vcf.gz"
    write_vcf(path, "chr20", child, phased, "HG002")
    back = read_calls(path, "chr20")
    assert {p: tuple(c.genotype) for p, c in back.items()} == phased
    assert all(c.phased and c.phase_set == "PATMAT" for c in back.values())
