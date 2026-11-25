#!/usr/bin/env bash
# Rebuild the untracked inputs of the genome-wide runs in data/ from their
# public sources, check them against scripts/data.md5, and derive the
# per-chromosome files and the transmission truth exactly as the runs used
# them. The reads are not fetched here: benchmarks/genome_run.py streams them
# itself (docs/design/genome-run-plan.md, docs/design/ont-run-plan.md).
#
#   bash scripts/fetch_data.sh
#
# Checked on 28 September 2026: the split rules below reproduce
# data/giab_wg/split byte for byte on chromosome 21, and the truth command
# reproduces data/giab_wg/truth on chromosome 21.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=.venv/bin/python
R=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio
UCSC=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes
mkdir -p data/giab_wg/split data/giab_wg/truth data/giab_v5 data/reference

fetch() {  # url target
  [ -s "$2" ] || curl -sfL --retry 3 -o "$2.part" "$1" && { [ -s "$2" ] || mv "$2.part" "$2"; }
}

for s in HG002_NA24385_son:HG002 HG003_NA24149_father:HG003 HG004_NA24143_mother:HG004; do
  dir=${s%%:*}; id=${s##*:}
  base=$R/$dir/NISTv4.2.1/GRCh38/${id}_GRCh38_1_22_v4.2.1_benchmark
  fetch "$base.vcf.gz" "data/giab_wg/${id}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
  fetch "$base.vcf.gz.tbi" "data/giab_wg/${id}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
  fetch "${base}_noinconsistent.bed" "data/giab_wg/${id}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
done
hifiasm=HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.vcf.gz
fetch "$R/HG002_NA24385_son/NISTv4.2.1/GRCh38/SupplementaryFiles/$hifiasm" "data/giab_wg/$hifiasm"
for f in HG002_GRCh38_v5.0q_smvar.vcf.gz HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi \
         HG002_GRCh38_v5.0q_smvar.benchmark.bed; do
  fetch "$R/HG002_NA24385_son/v5.0q/$f" "data/giab_v5/$f"
done
for i in $(seq 1 22); do
  if [ ! -s "data/reference/chr$i.fa" ]; then
    curl -sfL --retry 3 "$UCSC/chr$i.fa.gz" | gunzip > "data/reference/chr$i.fa.part"
    mv "data/reference/chr$i.fa.part" "data/reference/chr$i.fa"
  fi
  [ -s "data/reference/chr$i.fa.fai" ] || "$PY" -c "import pysam, sys; pysam.faidx(sys.argv[1])" \
    "data/reference/chr$i.fa"
done

md5sum -c scripts/data.md5

for i in $(seq 1 22); do
  c=chr$i
  for id in HG002 HG003 HG004; do
    out=data/giab_wg/split/${id}_${c}_v4.2.1.vcf.gz
    [ -s "$out" ] || "$PY" - "$id" "$c" "$out" <<'PY'
import gzip, sys, pysam
sample, chrom, out = sys.argv[1:]
source = pysam.VariantFile(f"data/giab_wg/{sample}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")
with pysam.VariantFile(out[:-3], "w", header=source.header) as target:
    for record in source.fetch(chrom):
        target.write(record)
with open(out[:-3], "rb") as plain, gzip.open(out, "wb", compresslevel=9) as packed:
    packed.write(plain.read())
PY
    rm -f "${out%.gz}"
    bed=data/giab_wg/split/${id}_${c}_v4.2.1.bed
    [ -s "$bed" ] || awk -v c="$c" '$1==c' \
      "data/giab_wg/${id}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" > "$bed"
  done
  out=data/giab_wg/split/HG002_${c}_hifiasm_phasetransfer.vcf.gz
  [ -s "$out" ] || zcat "data/giab_wg/$hifiasm" | awk -v c="$c" '/^#/ || $1==c' | gzip -9 > "$out"
  truth=data/giab_wg/truth/HG002_${c}_trio_truth.vcf.gz
  [ -s "$truth" ] || PYTHONPATH=src "$PY" benchmarks/build_trio_truth.py --chromosome "$c" \
    --data data/giab_wg/split --child "HG002_${c}_v4.2.1.vcf.gz" \
    --father "HG003_${c}_v4.2.1.vcf.gz" --mother "HG004_${c}_v4.2.1.vcf.gz" \
    --father-bed "HG003_${c}_v4.2.1.bed" --mother-bed "HG004_${c}_v4.2.1.bed" \
    --phased "HG002_${c}_hifiasm_phasetransfer.vcf.gz:PATMAT" --out "$truth" \
    > "data/giab_wg/truth/$c.log"
done
echo "data/ rebuilt and checked"
