"""Score every phased relation of any phaser's VCF by the reads in a BAM.

For each chromosome:
- keep the VCF's biallelic heterozygous SNVs, with their phase;
- read every read's alleles at those SNVs exactly as `whatshap phase` does:
  realignment against the reference, mapping quality 20, no supplementary
  alignments (`extract_fragments.read_fragments`);
- for each pair of consecutive phased SNVs in one phase set (PS), count the
  reads that agree and disagree with the called relation;
- score the pair by `phasebudget.confidence.score`, the Wilson lower bound
  (z = 2) on the agreeing share, ties broken by the number of reads.

The output is a TSV with one row per pair: chrom, pos1 of each SNV, agree,
disagree, score. Only the VCF's GT and PS fields are read, so the phaser can
be anything that writes them. A phased genotype with no PS is put in one
block per chromosome.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/phase_confidence.py
--bam reads.bam --vcf phased.vcf.gz --reference ref.fa --sample HG002 --out pairs.tsv`
"""

from __future__ import annotations

import argparse
import gzip
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from phasebudget._types import Phasing  # noqa: E402
from phasebudget.confidence import pair_support, score  # noqa: E402


def _extractor():
    import importlib.util

    spec = importlib.util.spec_from_file_location("extract_fragments",
                                                  HERE / "extract_fragments.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _records(vcf: str, sample: str, chromosome: str):
    """`(pos1, GT string, PS or None)` for `chromosome`, from an indexed or a plain file.

    Phasers often write plain or plain-gzipped VCFs without an index, which
    htslib cannot fetch from, so the file is read as text.
    """
    opener = gzip.open if vcf.endswith(".gz") else open
    column = None
    with opener(vcf, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                column = fields.index(sample)
                continue
            if fields[0] != chromosome:
                continue
            keys = fields[8].split(":")
            values = dict(zip(keys, fields[column].split(":")))
            ps = values.get("PS")
            yield int(fields[1]), values.get("GT", "."), (None if ps in (None, ".") else ps)


def phasing_from_vcf(vcf: str, sample: str, chromosome: str, positions) -> Phasing:
    """The VCF's first-haplotype alleles and PS blocks at `positions` (0-based)."""
    index = {p: j for j, p in enumerate(positions)}
    alleles: list = [None] * len(positions)
    blocks: list = [None] * len(positions)
    for pos1, genotype, ps in _records(vcf, sample, chromosome):
        j = index.get(pos1 - 1)
        if j is None or genotype not in ("0|1", "1|0"):
            continue
        alleles[j] = int(genotype[0])
        blocks[j] = 0 if ps is None else ps
    return Phasing(alleles=tuple(alleles), blocks=tuple(blocks))


def score_chromosome(bam: str, vcf: str, reference: str, sample: str, chromosome: str,
                     ignore_read_groups: bool = False):
    """`[(pos1 first, pos1 second, agree, disagree, score)]` for one chromosome."""
    extractor = _extractor()
    with tempfile.TemporaryDirectory() as tmp:
        snvs = str(Path(tmp) / "het.vcf.gz")
        if extractor.het_snvs(vcf, sample, chromosome, None, snvs) < 2:
            return []
        matrix, positions, _ = extractor.read_fragments(
            bam, snvs, reference, chromosome, sample, ignore_read_groups=ignore_read_groups)
        phasing = phasing_from_vcf(snvs, sample, chromosome, positions)
    return [(positions[s.first] + 1, positions[s.second] + 1, s.agree, s.disagree,
             score(s)[0]) for s in pair_support(matrix, phasing)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--vcf", required=True, help="any phaser's output, with GT and PS")
    parser.add_argument("--reference", required=True,
                        help="a FASTA, or a directory holding CHROM.fa per chromosome")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--chromosomes", nargs="+", help="default: every contig in the VCF")
    parser.add_argument("--ignore-read-groups", action="store_true")
    parser.add_argument("--out", required=True, help="TSV, gzipped if it ends in .gz")
    args = parser.parse_args()
    import pysam

    with pysam.VariantFile(args.vcf) as reader:
        chromosomes = args.chromosomes or list(reader.header.contigs)
    opener = gzip.open if args.out.endswith(".gz") else open
    with opener(args.out, "wt") as out:
        out.write("chrom\tpos1\tpos2\tagree\tdisagree\tscore\n")
        for chrom in chromosomes:
            reference = args.reference
            if Path(reference).is_dir():
                reference = str(Path(reference) / f"{chrom}.fa")
            for row in score_chromosome(args.bam, args.vcf, reference, args.sample, chrom,
                                        args.ignore_read_groups):
                out.write(f"{chrom}\t{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]:.6f}\n")
            print(f"{chrom}: scored", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
