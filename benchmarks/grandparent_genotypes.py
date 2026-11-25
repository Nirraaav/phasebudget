"""Z6 data, `docs/design/z6-grandparent-genotypes.md`: NA12891 and NA12892 at NA12878's sites.

Streams each parent's ONT BAM from the Platinum Pedigree bucket in 10 Mb
regions with `pysam.AlignmentFile.count_coverage`, one connection per 1 Mb
window, and keeps only reference and alternative read counts at NA12878's
heterozygous SNVs (amendment 1 of the statement says why not `samtools
mpileup`). Each region writes a marker, so the run resumes where it stopped.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/grandparent_genotypes.py --workers 3`
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import genome_common as common  # noqa: E402
from confidence_run import LENGTH_MB  # noqa: E402

SAMPLES = ("NA12891", "NA12892")
URL = "https://platinum-pedigree-data.s3.amazonaws.com/data/ont/mapped/GRCh38/{}.minimap2.bam"
OUT = common.ROOT / "data" / "z6"
REGION = 10_000_000
WINDOW = 1_000_000  # count_coverage window: bounded memory, little re-reading of long reads
MIN_DEPTH, HOM, HET_LOW, HET_HIGH = 8, 0.10, 0.25, 0.75
# (min depth, homozygous band, heterozygous low, heterozygous high); R0 is the rule fixed
# before the run, R1 and R2 the stricter rules of the repair in the design statement.
RULES = {"R0": (MIN_DEPTH, HOM, HET_LOW, HET_HIGH), "R1": (15, 0.05, 0.30, 0.70),
         "R2": (20, 0.03, 0.35, 0.65)}
INDEL = re.compile(r"[+-](\d+)")


def sites(chrom: str) -> dict:
    """`{pos1: (ref, alt)}` of NA12878's heterozygous SNVs, written once per chromosome."""
    path = OUT / "sites" / f"{chrom}.vcf.gz"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        extractor = common._module("extract_fragments.py")
        extractor.het_snvs(str(common.PEDIGREE / "CEPH1463.GRCh38.family-truthset.ov.vcf.gz"),
                           "NA12878", chrom, str(common.PEDIGREE / "hq_regions_final.bed"),
                           str(path) + ".part.vcf.gz", unphase=True)
        Path(str(path) + ".part.vcf.gz").rename(path)
        Path(str(path) + ".part.vcf.gz.tbi").rename(str(path) + ".tbi")
        Path(str(path) + ".part.vcf.gz.input.vcf").unlink(missing_ok=True)  # het_snvs' plain copy
    out = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("#"):
                f = line.split("\t", 6)
                out[int(f[1])] = (f[3], f[4])
    return out


def counts(bases: str, alt: str) -> tuple[int, int]:
    """Reference and `alt` reads in an mpileup base column."""
    ref = other = 0
    i = 0
    alt = alt.upper()
    while i < len(bases):
        c = bases[i]
        if c == "^":
            i += 2
            continue
        if c in "+-":
            m = INDEL.match(bases, i)
            i = m.end() + int(m.group(1))
            continue
        if c in ".,":
            ref += 1
        elif c.upper() == alt:
            other += 1
        i += 1
    return ref, other


def usable(read) -> bool:
    """The filters of `samtools mpileup -q 20 --ff UNMAP,SECONDARY,QCFAIL,DUP`."""
    return not (read.is_unmapped or read.is_secondary or read.is_qcfail or read.is_duplicate) \
        and read.mapping_quality >= 20


def genotype(ref: int, alt: int, rule: str = "R0"):
    """`(0, 0)`, `(0, 1)`, `(1, 1)` or None, by a rule of the design statement."""
    min_depth, hom, het_low, het_high = RULES[rule]
    d = ref + alt
    if d < min_depth:
        return None
    f = alt / d
    if f <= hom:
        return (0, 0)
    if f >= 1 - hom:
        return (1, 1)
    if het_low <= f <= het_high:
        return (0, 1)
    return None


def region(job) -> bool:
    import pysam

    sample, chrom, start = job
    done = OUT / sample / "done" / f"{chrom}_{start}.ok"
    if done.exists():
        return True
    known = sites(chrom)
    end = min(start + REGION, LENGTH_MB[chrom] * 1_000_000 + 10_000_000)
    wanted = sorted(p for p in known if start < p <= end)
    out = OUT / sample / f"{chrom}_{start}.tsv.gz"
    out.parent.mkdir(parents=True, exist_ok=True)
    done.parent.mkdir(parents=True, exist_ok=True)
    if not wanted:
        done.write_text("0 sites\n")
        return True
    for attempt in range(1, 7):
        began = time.time()
        here = os.getcwd()
        try:
            os.chdir(out.parent)  # htslib caches the remote index in the working directory
            with gzip.open(str(out) + ".part", "wt") as target:
                for w in range(start, end, WINDOW):
                    inside = [p for p in wanted if w < p <= w + WINDOW]
                    if not inside:
                        continue
                    # A new connection per window: long-lived S3 connections slowed to under
                    # 100 KB/s while fresh ones ran at several MB/s (development notebook, 8 October).
                    with pysam.AlignmentFile(URL.format(sample)) as bam:
                        acgt = bam.count_coverage(chrom, w, w + WINDOW, quality_threshold=0,
                                                  read_callback=usable)
                    for p in inside:
                        ref, alt = known[p]
                        target.write(f"{p}\t{acgt['ACGT'.index(ref)][p - 1 - w]}\t"
                                     f"{acgt['ACGT'.index(alt)][p - 1 - w]}\n")
            Path(str(out) + ".part").rename(out)
            done.write_text(f"{time.time() - began:.0f} s, {len(wanted)} sites\n")
            print(f"{sample} {chrom}:{start}: done in {time.time() - began:.0f} s", flush=True)
            return True
        except Exception:
            print(f"{sample} {chrom}:{start}: attempt {attempt} failed\n{traceback.format_exc()}",
                  flush=True)
            time.sleep(120)
        finally:
            os.chdir(here)
    return False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--chromosomes", nargs="+", default=list(common.CONFIRMATORY))
    args = parser.parse_args()
    jobs = [(s, c, start) for c in sorted(args.chromosomes, key=lambda c: LENGTH_MB[c])
            for start in range(0, LENGTH_MB[c] * 1_000_000 + 1, REGION) for s in SAMPLES]
    with ProcessPoolExecutor(args.workers) as pool:
        list(pool.map(sites, args.chromosomes))  # one writer per chromosome, before any region
        failed = [j for j, ok in zip(jobs, pool.map(region, jobs)) if not ok]
    if failed:
        raise SystemExit(f"failed regions: {len(failed)}")


if __name__ == "__main__":
    main()
