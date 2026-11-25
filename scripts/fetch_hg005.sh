#!/usr/bin/env bash
# Fetch and derive the truth inputs of the HG005 replication
# (docs/design/hg005-run-plan.md): the Chinese trio's v4.2.1 calls and BEDs,
# the hifiasm phase transfer, per-chromosome splits, and the transmission
# truth, by the same rules as scripts/fetch_data.sh uses for HG002.
#
#   bash scripts/fetch_hg005.sh
set -euo pipefail
cd "$(dirname "$0")/.."
PY=.venv/bin/python
R=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio
D=data/giab_hg005
mkdir -p "$D/split" "$D/truth"

fetch() {  # url target: download once, atomically
  if [ ! -s "$2" ]; then
    curl -sfL --retry 3 -o "$2.part" "$1"
    mv "$2.part" "$2"
  fi
}

for s in HG005_NA24631_son:HG005 HG006_NA24694_father:HG006 HG007_NA24695_mother:HG007; do
  dir=${s%%:*}; id=${s##*:}
  base=$R/$dir/NISTv4.2.1/GRCh38/${id}_GRCh38_1_22_v4.2.1_benchmark
  fetch "$base.vcf.gz" "$D/${id}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
  fetch "$base.vcf.gz.tbi" "$D/${id}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
  fetch "$base.bed" "$D/${id}_GRCh38_1_22_v4.2.1_benchmark.bed"
done
hifiasm=HG005_GRCh38_1_22_v4.2.1_highconf_hifiasm_v11_phasetransfer.vcf.gz
fetch "$R/HG005_NA24631_son/NISTv4.2.1/GRCh38/SupplementaryFiles/$hifiasm" "$D/$hifiasm"
md5sum -c scripts/data_hg005.md5

for i in $(seq 1 22); do
  c=chr$i
  for id in HG005 HG006 HG007; do
    out=$D/split/${id}_${c}_v4.2.1.vcf.gz
    if [ ! -s "$out" ]; then
      "$PY" - "$D/${id}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" "$c" "${out%.gz}" <<'PY'
import sys, pysam
source_path, chrom, out = sys.argv[1:]
source = pysam.VariantFile(source_path)
with pysam.VariantFile(out, "w", header=source.header) as target:
    for record in source.fetch(chrom):
        target.write(record)
PY
      gzip -9 "${out%.gz}"
    fi
    bed=$D/split/${id}_${c}_v4.2.1.bed
    [ -s "$bed" ] || awk -v c="$c" '$1==c' "$D/${id}_GRCh38_1_22_v4.2.1_benchmark.bed" > "$bed"
  done
  out=$D/split/HG005_${c}_hifiasm_phasetransfer.vcf.gz
  [ -s "$out" ] || zcat "$D/$hifiasm" | awk -v c="$c" '/^#/ || $1==c' | gzip -9 > "$out"
  truth=$D/truth/HG005_${c}_trio_truth.vcf.gz
  [ -s "$truth" ] || PYTHONPATH=src "$PY" benchmarks/build_trio_truth.py --chromosome "$c" \
    --data "$D/split" --child "HG005_${c}_v4.2.1.vcf.gz" \
    --father "HG006_${c}_v4.2.1.vcf.gz" --mother "HG007_${c}_v4.2.1.vcf.gz" \
    --father-bed "HG006_${c}_v4.2.1.bed" --mother-bed "HG007_${c}_v4.2.1.bed" \
    --phased "HG005_${c}_hifiasm_phasetransfer.vcf.gz:PATMAT" --out "$truth" --sample HG005 \
    > "$D/truth/$c.log"
done
echo "HG005 truth inputs ready"
