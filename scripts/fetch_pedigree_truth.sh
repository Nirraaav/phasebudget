#!/usr/bin/env bash
# Fetch the Platinum Pedigree small-variant truthset (Kronenberg et al. 2025),
# GRCh38, phased by inheritance across CEPH-1463, and its high-quality regions;
# then check them against scripts/data_pedigree_truth.md5. Used by
# docs/design/truth-adjudication-plan.md and the extension's pedigree arm.
#
#   bash scripts/fetch_pedigree_truth.sh
set -euo pipefail
cd "$(dirname "$0")/.."
B=https://platinum-pedigree-data.s3.amazonaws.com/variants/small_variant_truthset/GRCh38
D=data/pedigree_truth
mkdir -p "$D"
for f in CEPH1463.GRCh38.family-truthset.ov.vcf.gz CEPH1463.GRCh38.family-truthset.ov.vcf.gz.tbi \
         hq_regions_final.bed.gz; do
  if [ ! -s "$D/$f" ]; then
    curl -sfL --retry 3 -o "$D/$f.part" "$B/$f"
    mv "$D/$f.part" "$D/$f"
  fi
done
md5sum -c scripts/data_pedigree_truth.md5
