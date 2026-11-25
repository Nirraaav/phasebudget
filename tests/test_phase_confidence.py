"""The BAM-and-VCF scorer on a hand-built BAM where every read's alleles are known."""

import importlib.util
import random
from pathlib import Path

import pytest

pysam = pytest.importorskip("pysam")
pytest.importorskip("whatshap")

ROOT = Path(__file__).parents[1]
CHROM = "chr20"
SNVS = [(300, "A", "G"), (700, "C", "T"), (1100, "G", "A"), (1500, "T", "C")]  # 1-based


def _load():
    path = ROOT / "benchmarks" / "phase_confidence.py"
    spec = importlib.util.spec_from_file_location("phase_confidence", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _inputs(tmp_path, genotypes):
    rng = random.Random(5)
    bases = [rng.choice("ACGT") for _ in range(2000)]
    for position, ref, _ in SNVS:
        bases[position - 1] = ref
    reference = "".join(bases)
    fasta = tmp_path / "ref.fa"
    fasta.write_text(f">{CHROM}\n" + "\n".join(reference[i:i + 60]
                                               for i in range(0, len(reference), 60)) + "\n")
    pysam.faidx(str(fasta))
    lines = ["##fileformat=VCFv4.2", f"##contig=<ID={CHROM},length={len(reference)}>",
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002"]
    for (position, ref, alt), (gt, ps) in zip(SNVS, genotypes):
        lines.append(f"{CHROM}\t{position}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT:PS\t{gt}:{ps}")
    vcf = tmp_path / "phased.vcf"
    vcf.write_text("\n".join(lines) + "\n")
    pysam.tabix_index(str(vcf), preset="vcf", force=True)
    # (name, start, end, alleles at the SNVs it covers, in order)
    reads = [("h0_long", 100, 1700, [0, 0, 0, 0]), ("h1_long", 100, 1700, [1, 1, 1, 1]),
             ("h0_mid", 500, 1300, [0, 0]), ("chimera", 200, 900, [0, 1])]
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CHROM, "LN": len(reference)}],
              "RG": [{"ID": "rg", "SM": "HG002"}]}
    bam = tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for name, start, end, alleles in sorted(reads, key=lambda r: r[1]):
            sequence = list(reference[start:end])
            covered = [s for s in SNVS if start < s[0] <= end]
            for (position, _, alt), allele in zip(covered, alleles):
                if allele == 1:
                    sequence[position - 1 - start] = alt
            record = pysam.AlignedSegment(out.header)
            record.query_name = name
            record.query_sequence = "".join(sequence)
            record.reference_id = 0
            record.reference_start = start
            record.cigarstring = f"{end - start}M"
            record.mapping_quality = 60
            record.query_qualities = pysam.qualitystring_to_array("I" * (end - start))
            record.set_tag("RG", "rg")
            out.write(record)
    pysam.index(str(bam))
    return str(fasta), str(vcf) + ".gz", str(bam)


def test_support_counts_each_read_once_per_pair(tmp_path):
    fasta, vcf, bam = _inputs(tmp_path, [("0|1", 300)] * 4)
    rows = _load().score_chromosome(bam, vcf, fasta, "HG002", CHROM)
    support = {(a, b): (agree, disagree) for a, b, agree, disagree, _ in rows}
    # The chimera calls 0 at 300 and 1 at 700, against the called same-haplotype relation.
    assert support == {(300, 700): (2, 1), (700, 1100): (3, 0), (1100, 1500): (2, 0)}
    scores = {(a, b): s for a, b, _, _, s in rows}
    assert scores[(300, 700)] < scores[(1100, 1500)] < scores[(700, 1100)]


def test_phase_sets_and_unphased_genotypes_split_pairs(tmp_path):
    genotypes = [("0|1", 300), ("1|0", 300), ("0/1", 300), ("0|1", 1500)]
    fasta, vcf, bam = _inputs(tmp_path, genotypes)
    rows = _load().score_chromosome(bam, vcf, fasta, "HG002", CHROM)
    # 1100 is unphased and 1500 opens a new phase set: only (300, 700) remains, now with
    # the opposite-haplotype relation, which the two long reads contradict.
    assert [(a, b, agree, disagree) for a, b, agree, disagree, _ in rows] == [(300, 700, 1, 2)]
