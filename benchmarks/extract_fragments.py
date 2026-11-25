"""Turn an aligned BAM and a VCF into a phasebudget fragment matrix.

This is the fragment step of `docs/design/real-data-plan.md`. It keeps the
biallelic heterozygous SNVs of one sample inside a BED, detects each read's
alleles with WhatsHap 2.8's `PhasedInputReader`, configured as `whatshap
phase` configures it (realignment against the reference, mapping quality 20,
no supplementary alignments), and converts the resulting `ReadSet` into a
`FragmentMatrix`: variant `j` is the `j`-th kept SNV in position order, and a
read's weights are the integer allele qualities WhatsHap reports, unscaled.

Reads calling one variant are kept, for coverage accounting; alleles WhatsHap
could not resolve (anything other than 0 or 1) are dropped. Nothing else is
filtered. The output is a gzipped JSON file that `load_fragments` reads back.

WhatsHap is an optional dependency (`pip install -e ".[whatshap]"`).
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import platform
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Optional, Sequence

from phasebudget._types import Fragment
from phasebudget.fragments import FragmentMatrix

MAX_INDEL = 50  # bases; the Platinum Pedigree truthset keeps indels shorter than 51


def het_snvs(vcf: str, sample: str, chromosome: str, bed: Optional[str], out: str,
             indels: bool = False, unphase: bool = False) -> int:
    """Write the biallelic heterozygous SNVs of `sample` on `chromosome` inside `bed`.

    Returns the number written. The genotype must be 0/1 or 1/0, phased or
    not. Positions in the BED are 0-based half-open, as in the BED format.
    With `indels`, biallelic indels of at most `MAX_INDEL` bases are kept too.
    With `unphase`, written genotypes lose any phase, so no phaser can read the
    truth's phase from its input.
    """
    import pysam

    intervals: list[tuple[int, int]] = []
    if bed is not None:
        with open(bed) as handle:
            for line in handle:
                fields = line.split()
                if fields and fields[0] == chromosome:
                    intervals.append((int(fields[1]), int(fields[2])))
        intervals.sort()

    def inside(position: int) -> bool:  # 0-based
        lo, hi = 0, len(intervals)
        while lo < hi:
            mid = (lo + hi) // 2
            if intervals[mid][1] <= position:
                lo = mid + 1
            else:
                hi = mid
        return lo < len(intervals) and intervals[lo][0] <= position

    if vcf.endswith(".gz"):
        # The chromosome extracts are plain gzip, which htslib cannot open as a
        # VCF, so they are decompressed beside the output first.
        plain = out + ".input.vcf"
        with gzip.open(vcf, "rb") as source, open(plain, "wb") as target:
            shutil.copyfileobj(source, target)
        vcf = plain
    reader = pysam.VariantFile(vcf)
    reader.subset_samples([sample])  # a family truthset has every member; phasers need one
    header = reader.header.copy()
    writer = pysam.VariantFile(out, "wz", header=header)
    kept = 0
    # Streamed rather than fetched, since a decompressed file has no index.
    for record in reader:
        if record.chrom != chromosome:
            continue
        if record.alts is None or len(record.alts) != 1 or record.alts[0].startswith("<"):
            continue
        if indels:
            if abs(len(record.ref) - len(record.alts[0])) > MAX_INDEL or \
                    (len(record.ref) > 1 and len(record.alts[0]) > 1):
                continue
        elif len(record.ref) != 1 or len(record.alts[0]) != 1:
            continue
        genotype = record.samples[sample]["GT"]
        if genotype is None or sorted(a for a in genotype if a is not None) != [0, 1]:
            continue
        if intervals and not inside(record.pos - 1):
            continue
        if unphase:
            record.samples[sample]["GT"] = (0, 1)  # the truth's allele order is its phase
            record.samples[sample].phased = False
        writer.write(record)
        kept += 1
    writer.close()
    if unphase:
        # Keep GT alone: a truth's other FORMAT fields can carry phase (GIAB v4.2.1 has IGT,
        # IPS and PS; amendment 2 of docs/design/stratified-run-plan.md).
        plain = out + ".gtonly.vcf"
        with gzip.open(out, "rt") as source, open(plain, "w") as target:
            for line in source:
                if line.startswith("##FORMAT=") and not line.startswith("##FORMAT=<ID=GT,"):
                    continue
                if not line.startswith("#"):
                    f = line.rstrip("\n").split("\t")
                    line = "\t".join(f[:8] + ["GT", f[9].split(":")[0]]) + "\n"
                target.write(line)
        pysam.tabix_compress(plain, out, force=True)
        Path(plain).unlink()
    pysam.tabix_index(out, preset="vcf", force=True)
    return kept


def matrix_from_readset(readset, positions: Sequence[int]):
    """`(FragmentMatrix, names)` from a WhatsHap ReadSet over `positions`.

    Read `i` of the ReadSet becomes read id `i`. Positions are WhatsHap's,
    0-based.
    """
    index = {position: j for j, position in enumerate(positions)}
    fragments = []
    names = []
    for read in readset:
        calls = sorted(
            (index[variant.position], int(variant.allele), float(variant.quality))
            for variant in read
            if variant.position in index and variant.allele in (0, 1)
        )
        if not calls:
            continue
        start = calls[0][0]
        width = calls[-1][0] - start + 1
        alleles: list = [None] * width
        weights = [0.0] * width
        for j, allele, quality in calls:
            alleles[j - start] = allele
            weights[j - start] = quality
        fragments.append(Fragment(len(fragments), start, tuple(alleles), tuple(weights)))
        names.append(read.name)
    return FragmentMatrix(fragments, len(positions)), names


def read_fragments(bam: str, vcf: str, reference: Optional[str], chromosome: str, sample: str,
                   mapq: int = 20, ignore_read_groups: bool = False, indels: bool = False):
    """`(matrix, positions, names)` for the variants in `vcf`, as `whatshap phase` reads them.

    With a reference, alleles come from WhatsHap's realignment, whose quality
    is a constant 30. With `reference=None` they come from the alignment's
    CIGAR, as `whatshap phase` without `--reference` does, and each allele's
    quality is the read's own base quality at the SNV.
    """
    from whatshap.cli import PhasedInputReader
    from whatshap.core import NumericSampleIds
    from whatshap.vcf import VcfReader

    ids = NumericSampleIds()
    with VcfReader(vcf, only_snvs=not indels) as vcf_reader:
        table = vcf_reader.fetch(chromosome)
    variants = table.variants
    positions = [variant.position for variant in variants]
    with PhasedInputReader([bam], reference, ids, ignore_read_groups=ignore_read_groups,
                           mapq_threshold=mapq, only_snvs=not indels,
                           use_supplementary=False) as reader:
        readset, _ = reader.read(chromosome, variants, sample, read_vcf=False)
    matrix, names = matrix_from_readset(readset, positions)
    return matrix, positions, names


def save_fragments(path: str, matrix: FragmentMatrix, positions: Sequence[int],
                   names: Sequence[str], provenance: dict) -> None:
    payload = {
        "provenance": provenance,
        "n_variants": matrix.n_variants,
        "positions": list(positions),
        "names": list(names),
        "fragments": [[f.start, [a if a is not None else -1 for a in f.alleles],
                       list(f.weights)] for f in matrix],
    }
    with gzip.open(path, "wt") as handle:
        json.dump(payload, handle)


def load_fragments(path: str):
    """`(matrix, positions, names, provenance)` from a file `save_fragments` wrote."""
    with gzip.open(path, "rt") as handle:
        payload = json.load(handle)
    fragments = [
        Fragment(i, start, tuple(None if a < 0 else a for a in alleles), tuple(weights))
        for i, (start, alleles, weights) in enumerate(payload["fragments"])
    ]
    matrix = FragmentMatrix(fragments, payload["n_variants"])
    return matrix, payload["positions"], payload["names"], payload["provenance"]


def phred_to_log_odds(matrix: FragmentMatrix) -> FragmentMatrix:
    """The same matrix with each phred quality `q` replaced by `log((1-e)/e)`, `e = 10^(-q/10)`.

    WhatsHap reports allele qualities on the phred scale; the library's
    weights are natural log-odds (`simulate.log_odds_weight`), the scale the
    posterior and the loop's error-probability rule read. A quality at or
    below the coin-flip phred, about 3, becomes weight zero. Gaps stay zero.
    """
    from phasebudget.simulate import COIN_FLIP_PHRED, log_odds_weight, phred_to_error

    def convert(allele, quality):
        if allele is None or quality <= COIN_FLIP_PHRED:
            return 0.0
        return log_odds_weight(phred_to_error(quality))

    fragments = [
        Fragment(f.read_id, f.start, f.alleles,
                 tuple(convert(a, q) for a, q in zip(f.alleles, f.weights)))
        for f in matrix
    ]
    return FragmentMatrix(fragments, matrix.n_variants)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--bed")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--alleles", choices=["realign", "basequality"], default="realign",
                        help="realign: WhatsHap's default realignment (constant quality 30); "
                             "basequality: alleles from the alignment, weighted by base quality")
    parser.add_argument("--chromosome", default="chr20")
    parser.add_argument("--sample", default="HG002")
    parser.add_argument("--mapq", type=int, default=20)
    parser.add_argument("--ignore-read-groups", action="store_true",
                        help="treat every read as the sample's, for BAMs without read groups")
    parser.add_argument("--indels", action="store_true",
                        help="keep biallelic indels as well as SNVs (extension track B)")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    import whatshap

    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    started = time.time()
    with tempfile.TemporaryDirectory() as scratch:
        snvs = os.path.join(scratch, "het_snvs.vcf.gz")
        kept = het_snvs(args.vcf, args.sample, args.chromosome, args.bed, snvs, args.indels)
        reference = args.reference if args.alleles == "realign" else None
        matrix, positions, names = read_fragments(args.bam, snvs, reference,
                                                  args.chromosome, args.sample, args.mapq,
                                                  args.ignore_read_groups, args.indels)
    provenance = {
        "commit": commit, "whatshap": whatshap.__version__,
        "machine": f"{platform.system()} {platform.release()} {platform.machine()}",
        "bam": Path(args.bam).name, "vcf": Path(args.vcf).name,
        "bed": Path(args.bed).name if args.bed else None,
        "reference": Path(args.reference).name, "chromosome": args.chromosome,
        "sample": args.sample, "mapq": args.mapq,
        "ignore_read_groups": args.ignore_read_groups, "het_snvs": kept, "alleles": args.alleles,
        "indels": args.indels,
        "seconds": round(time.time() - started, 1),
    }
    save_fragments(args.out, matrix, positions, names, provenance)
    informative = sum(1 for f in matrix if f.n_observed >= 2)
    print(json.dumps(provenance))
    print(f"{kept} heterozygous SNVs; {len(matrix)} reads call at least one, "
          f"{informative} call two or more; written to {args.out}")


if __name__ == "__main__":
    main()
