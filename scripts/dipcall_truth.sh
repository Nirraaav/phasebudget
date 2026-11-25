#!/usr/bin/env bash
# Build a phased truth (dipcall VCF and BED) for one sample from its paternal
# and maternal assemblies, against GRCh38 no-alt. Outputs go to
# data/hprc/<sample>/. Round 4, direction W1 (the October directions memo (development repository)).
#
# Usage: scripts/dipcall_truth.sh <sample> <pat.fa.gz URL> <mat.fa.gz URL>
# Environment: THREADS (default 12), MALE=1 for a male sample (PAR handling),
# DIPCALL_KIT (default ~/.local/opt/dipcall.kit, dipcall 0.3).
set -euo pipefail
sample=$1 pat_url=$2 mat_url=$3
root=$(cd "$(dirname "$0")/.." && pwd)
kit=${DIPCALL_KIT:-$HOME/.local/opt/dipcall.kit}
ref=$root/data/reference/GRCh38_no_alt_analysis_set.fna
out=$root/data/hprc/$sample

fetch() {  # url dest
    if [ ! -s "$2" ]; then
        curl -sfL --retry 5 -o "$2.part" "$1"
        mv "$2.part" "$2"
    fi
}

if [ ! -s "$ref.fai" ]; then
    fetch "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" "$ref.gz"
    gzip -dc "$ref.gz" > "$ref.part" && mv "$ref.part" "$ref"
    "$kit/samtools" faidx "$ref"
fi

mkdir -p "$out"
cd "$out"
fetch "$pat_url" pat.fa.gz
fetch "$mat_url" mat.fa.gz
par=()
[ "${MALE:-0}" = 1 ] && par=(-x "$kit/hs38.PAR.bed")
"$kit/run-dipcall" -t "${THREADS:-12}" "${par[@]}" "$sample" "$ref" pat.fa.gz mat.fa.gz > "$sample.mak"
make -j2 -f "$sample.mak"
md5sum pat.fa.gz mat.fa.gz "$sample.dip.vcf.gz" "$sample.dip.bed" > "$sample.md5"
