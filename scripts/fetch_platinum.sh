#!/usr/bin/env bash
# Fetch the HG001 (NA12878) truth of the error-controlled phasing test
# (docs/design/confidence-run-plan.md): Illumina Platinum Genomes 2017-1.0,
# GRCh38 small variants, phased by pedigree inheritance, and the confident
# regions; then check them against scripts/data_platinum.md5.
#
#   bash scripts/fetch_platinum.sh
set -euo pipefail
cd "$(dirname "$0")/.."
B=https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants
D=data/platinum
mkdir -p "$D"

fetch() {  # url target: download once, atomically
  if [ ! -s "$2" ]; then
    curl -sfL --retry 3 -o "$2.part" "$1"
    mv "$2.part" "$2"
  fi
}

fetch "$B/NA12878/NA12878.vcf.gz" "$D/NA12878.vcf.gz"
fetch "$B/NA12878/NA12878.vcf.gz.tbi" "$D/NA12878.vcf.gz.tbi"
fetch "$B/ConfidentRegions.bed.gz" "$D/ConfidentRegions.bed.gz"
fetch "$B/ConfidentRegions.bed.gz.tbi" "$D/ConfidentRegions.bed.gz.tbi"
md5sum -c scripts/data_platinum.md5

# Exploratory, after the registered test was scored (development notebook, 4 Oct): GIAB
# HG001 v4.2.1, unphased, used only for its benchmark regions and genotypes.
G=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38
mkdir -p data/giab_hg001
for f in HG001_GRCh38_1_22_v4.2.1_benchmark.bed HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
         HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi; do
  fetch "$G/$f" "data/giab_hg001/$f"
done
md5sum -c scripts/data_hg001_giab.md5
