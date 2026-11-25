"""L7 (docs/design/l7-na12878-plan.md): a deterministic ~30x subsample of a remote HiFi BAM as FASTQ.

aria2c could not finish the NA12878 BAM (TLS handshake failures on 22 Nov), so reads are drawn
through htslib, one 10 Mb chunk per job. A read belongs to the chunk holding its start, so no read
is written twice; secondary and supplementary alignments are skipped, and unmapped reads are not
reached (deviation from K1b, which used every read). A read is kept when crc32 of
"<seed>:<name>" over 2^32 is below the fraction, so the subsample does not depend on chunking.

Usage:
  python l7_reads.py plan URL SEED TARGET_X OUT_DIR   # writes OUT_DIR/chunks.tsv and fraction.txt
  python l7_reads.py chunk URL SEED FRACTION CONTIG START END OUT.fq.gz
"""

from __future__ import annotations

import gzip
import sys
import zlib
from pathlib import Path

import pysam

CHUNK = 10_000_000
CONTIGS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def _close(bam) -> None:
    """Closing a remote file can raise after every read succeeded (seen 22 Nov); the data are complete."""
    try:
        bam.close()
    except OSError:
        pass


def keep(name: str, seed: int, fraction: float) -> bool:
    return zlib.crc32(f"{seed}:{name}".encode()) / 2**32 < fraction


def plan(url: str, seed: int, target: float, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    bam = pysam.AlignmentFile(url)
    if True:
        lengths = dict(zip(bam.references, bam.lengths))
        mapped = sum(s.mapped for s in bam.get_index_statistics())
        total, n = 0, 0
        for read in bam.fetch("chr1", 50_000_000, 60_000_000):
            if read.is_secondary or read.is_supplementary:
                continue
            total += read.query_length or read.infer_read_length() or 0
            n += 1
            if n == 20_000:
                break
    _close(bam)
    mean = total / n
    fraction = min(1.0, target * 3.1e9 / (mapped * mean))
    (out / "fraction.txt").write_text(
        f"mapped records {mapped}, mean read length {mean:.0f} over {n} reads, target {target}x, "
        f"fraction {fraction:.4f}, seed {seed}\n{fraction:.6f}\n")
    with open(out / "chunks.tsv", "w") as handle:
        for c in CONTIGS:
            for start in range(0, lengths[c], CHUNK):
                handle.write(f"{c}\t{start}\t{min(start + CHUNK, lengths[c])}\n")


def chunk(url: str, seed: int, fraction: float, contig: str, start: int, end: int, out: Path) -> None:
    part = out.with_name(out.name + ".part")
    written = 0
    bam = pysam.AlignmentFile(url)
    with gzip.open(part, "wt", compresslevel=1) as fq:
        for read in bam.fetch(contig, start, end):
            if read.is_secondary or read.is_supplementary or not start <= read.reference_start < end:
                continue
            if not keep(read.query_name, seed, fraction):
                continue
            seq, qual = read.get_forward_sequence(), read.get_forward_qualities()
            q = "".join(chr(x + 33) for x in qual) if qual is not None else "I" * len(seq)
            fq.write(f"@{read.query_name}\n{seq}\n+\n{q}\n")
            written += 1
    _close(bam)
    part.replace(out)
    print(f"{contig}:{start}-{end} {written}", flush=True)


def _check() -> None:
    assert keep("a", 78, 1.0) and not keep("a", 78, 0.0)
    kept = sum(keep(f"r{i}", 78, 0.4) for i in range(100_000))
    assert 39_000 < kept < 41_000, kept


if __name__ == "__main__":
    if sys.argv[1] == "check":
        _check()
    elif sys.argv[1] == "plan":
        plan(sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), Path(sys.argv[5]))
    else:
        _, _, url, seed, fraction, contig, start, end, out = sys.argv
        chunk(url, int(seed), float(fraction), contig, int(start), int(end), Path(out))
