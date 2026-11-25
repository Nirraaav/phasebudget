#!/usr/bin/env bash
# Round 4, W4: realign HG002's HPRC HiFi reads with one aligner (minimap2 map-hifi) to GRCh38
# no-alt and to T2T-CHM13v2.0, so the two references are compared with the aligner held fixed.
# The reads are the primary records of HPRC's winnowmap GRCh38 BAM (8,480,927 mapped records and
# 1,963 unmapped on 14 Oct). Outputs in data/w4/; per-chromosome CHM13 FASTAs in
# data/reference_chm13/ for the phasers.
# Usage: scripts/realign_hifi.sh   (THREADS per alignment, default 32)
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
py=$root/.venv/bin/python
url=https://human-pangenomics.s3.amazonaws.com/working/HPRC_PLUS/HG002/analysis/aligned_reads/hifi/GRCh38/HG002_aligned_GRCh38_winnowmap.sorted.bam
chm13_url=https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
grch38=$root/data/reference/GRCh38_no_alt_analysis_set.fna
refdir=$root/data/reference_chm13
chm13=$refdir/chm13v2.0.fa
out=$root/data/w4
mkdir -p "$refdir" "$out"

if [ ! -s "$chm13.fai" ]; then
    curl -sfL --retry 5 -o "$chm13.gz.part" "$chm13_url"
    gzip -dc "$chm13.gz.part" > "$chm13.part" && mv "$chm13.part" "$chm13" && rm "$chm13.gz.part"
    "$py" -c "import pysam, sys; pysam.faidx(sys.argv[1])" "$chm13"
fi
for c in chr{1..22}; do
    if [ ! -s "$refdir/$c.fa.fai" ]; then
        "$py" -c "import pysam, sys; pysam.faidx(sys.argv[1], sys.argv[2], '-o', sys.argv[3])" \
            "$chm13" "$c" "$refdir/$c.fa"
        "$py" -c "import pysam, sys; pysam.faidx(sys.argv[1])" "$refdir/$c.fa"
    fi
done

cd "$out"
# One stream per contig: a single remote stream runs at about 2 MB/s here.
if [ ! -s reads.fq.gz ]; then
    mkdir -p fq
    "$py" -c "import pysam, sys; print('\n'.join(pysam.AlignmentFile(sys.argv[1]).references + ('*',)))" "$url" > contigs.txt
    export py url
    xargs -P 24 -I{} bash -c '[ -s "fq/{}.done" ] || { "$py" -c "import pysam, sys; pysam.view(\"-u\", \"-F\", \"0x900\", \"-o\", sys.argv[2] + \".bam\", sys.argv[1], sys.argv[3], catch_stdout=False); pysam.samtools.fastq(\"-0\", sys.argv[2], sys.argv[2] + \".bam\", catch_stdout=False)" "$url" "fq/{}.fq.gz" "{}" && rm "fq/{}.fq.gz.bam" && echo > "fq/{}.done"; }' < contigs.txt
    [ "$(ls fq/*.done | wc -l)" -eq "$(wc -l < contigs.txt)" ]
    while read -r c; do cat "fq/$c.fq.gz"; done < contigs.txt > reads.part.fq.gz
    mv reads.part.fq.gz reads.fq.gz
    rm -r fq
fi
align() {  # name reference
    [ -s "HG002_$1_mm2.bam.bai" ] && return 0
    "$HOME/.local/bin/minimap2" -ax map-hifi -t "${THREADS:-32}" "$2" reads.fq.gz 2> "HG002_$1_mm2.log" |
        samtools sort -@ 8 -m 4G -T "sort_$1" -o "HG002_$1_mm2.part.bam" -
    mv "HG002_$1_mm2.part.bam" "HG002_$1_mm2.bam"
    samtools index "HG002_$1_mm2.bam"
}
align grch38 "$grch38" & first=$!
align chm13 "$chm13" & second=$!
wait "$first"  # a bare `wait` returns 0 even when an alignment failed
wait "$second"
md5sum reads.fq.gz > reads.fq.gz.md5
