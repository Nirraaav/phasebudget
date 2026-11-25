"""The fragment extractor on a BAM built by hand, where every allele is known."""

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
    path = ROOT / "benchmarks" / "extract_fragments.py"
    spec = importlib.util.spec_from_file_location("extract_fragments", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _reference():
    rng = random.Random(5)
    bases = [rng.choice("ACGT") for _ in range(2000)]
    for position, ref, _ in SNVS:
        bases[position - 1] = ref
    return "".join(bases)


def _write_inputs(tmp_path):
    reference = _reference()
    fasta = tmp_path / "ref.fa"
    fasta.write_text(f">{CHROM}\n" + "\n".join(reference[i:i + 60]
                                               for i in range(0, len(reference), 60)) + "\n")
    pysam.faidx(str(fasta))

    vcf = tmp_path / "calls.vcf"
    lines = ["##fileformat=VCFv4.2", f"##contig=<ID={CHROM},length={len(reference)}>",
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002"]
    for position, ref, alt in SNVS:
        lines.append(f"{CHROM}\t{position}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t0/1")
    lines.append(f"{CHROM}\t1800\t.\tA\tAT\t50\tPASS\t.\tGT\t0/1")   # an indel, dropped
    lines.append(f"{CHROM}\t1900\t.\tC\tG\t50\tPASS\t.\tGT\t1/1")    # homozygous, dropped
    vcf.write_text("\n".join(lines) + "\n")
    pysam.tabix_index(str(vcf), preset="vcf", force=True)

    # (name, start, end, alleles at the SNVs it covers in order, mapq)
    reads = [
        ("h0_long", 100, 1700, [0, 0, 0, 0], 60),
        ("h1_long", 100, 1700, [1, 1, 1, 1], 60),
        ("h0_mid", 500, 1300, [0, 0], 60),       # covers 700 and 1100
        ("single", 200, 500, [1], 60),           # covers 300 only
        ("lowmapq", 100, 1700, [0, 0, 0, 0], 5),  # removed by the mapq filter
    ]
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CHROM, "LN": len(reference)}],
              "RG": [{"ID": "rg", "SM": "HG002"}]}
    bam = tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for name, start, end, alleles, mapq in sorted(reads, key=lambda r: r[1]):
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
            record.mapping_quality = mapq
            record.query_qualities = pysam.qualitystring_to_array("I" * (end - start))
            record.set_tag("RG", "rg")
            out.write(record)
    pysam.index(str(bam))
    return str(fasta), str(vcf) + ".gz", str(bam)


def test_extraction_recovers_the_alleles_that_were_written(tmp_path):
    module = _load()
    fasta, vcf, bam = _write_inputs(tmp_path)
    snvs = str(tmp_path / "het.vcf.gz")
    assert module.het_snvs(vcf, "HG002", CHROM, None, snvs) == 4
    matrix, positions, names = module.read_fragments(bam, snvs, fasta, CHROM, "HG002")
    assert positions == [p - 1 for p, _, _ in SNVS]
    by_name = {name: fragment for name, fragment in zip(names, matrix)}
    assert set(by_name) == {"h0_long", "h1_long", "h0_mid", "single"}
    assert by_name["h0_long"].alleles == (0, 0, 0, 0)
    assert by_name["h1_long"].alleles == (1, 1, 1, 1)
    assert (by_name["h0_mid"].start, by_name["h0_mid"].alleles) == (1, (0, 0))
    assert (by_name["single"].start, by_name["single"].alleles) == (0, (1,))
    for fragment in matrix:
        assert all(weight > 0 for weight in fragment.weights)


def test_the_bed_restricts_the_variants(tmp_path):
    module = _load()
    _, vcf, _ = _write_inputs(tmp_path)
    bed = tmp_path / "regions.bed"
    bed.write_text(f"{CHROM}\t600\t1200\n")   # 0-based: holds 700 and 1100
    assert module.het_snvs(vcf, "HG002", CHROM, str(bed), str(tmp_path / "x.vcf.gz")) == 2


def test_saved_fragments_load_back_unchanged(tmp_path):
    module = _load()
    fasta, vcf, bam = _write_inputs(tmp_path)
    snvs = str(tmp_path / "het.vcf.gz")
    module.het_snvs(vcf, "HG002", CHROM, None, snvs)
    matrix, positions, names = module.read_fragments(bam, snvs, fasta, CHROM, "HG002")
    path = str(tmp_path / "fragments.json.gz")
    module.save_fragments(path, matrix, positions, names, {"test": True})
    loaded, loaded_positions, loaded_names, provenance = module.load_fragments(path)
    assert list(loaded) == list(matrix)
    assert (loaded_positions, loaded_names, provenance) == (positions, names, {"test": True})


def test_phred_qualities_become_log_odds():
    import math

    from phasebudget import FragmentMatrix

    module = _load()
    matrix = FragmentMatrix.from_rows(["01-1"], weights=[[30, 20, 0, 3]])
    converted = module.phred_to_log_odds(matrix)
    weights = converted[0].weights
    assert weights[0] == pytest.approx(math.log(999))
    assert weights[1] == pytest.approx(math.log(99))
    assert weights[2] == 0.0          # the gap
    assert weights[3] == 0.0          # phred 3 is worse than a coin flip's 3.01
    assert converted[0].alleles == matrix[0].alleles


def test_indels_are_kept_on_request_and_phase_is_removed(tmp_path):
    module = _load()
    _, vcf, _ = _write_inputs(tmp_path)
    phased = tmp_path / "phased.vcf"
    lines = [line.replace("\t0/1", "\t1|0") for line in pysam.BGZFile(vcf).read().decode().splitlines()]
    phased.write_text("\n".join(lines) + "\n")
    pysam.tabix_index(str(phased), preset="vcf", force=True)
    out = str(tmp_path / "het.vcf.gz")
    assert module.het_snvs(str(phased) + ".gz", "HG002", CHROM, None, out) == 4
    assert module.het_snvs(str(phased) + ".gz", "HG002", CHROM, None, out, indels=True,
                           unphase=True) == 5
    with pysam.VariantFile(out) as reader:
        records = list(reader)
    assert (records[-1].pos, records[-1].ref, records[-1].alts) == (1800, "A", ("AT",))
    for record in records:
        assert record.samples["HG002"]["GT"] == (0, 1)
        assert not record.samples["HG002"].phased
