#!/usr/bin/env bash
# Regenerate the run outputs behind Part A of docs/submission/paper.tex
# from public data, in order, then build the paper.
# Each step is idempotent: finished jobs are skipped. The runs take days of
# machine time and about 60 GB of free disk; run them one at a time.
#
#   bash scripts/reproduce_paper.sh            everything
#   bash scripts/reproduce_paper.sh score      only the scoring and figures
set -euo pipefail
cd "$(dirname "$0")/.."
PY=.venv/bin/python
export PYTHONPATH=src:benchmarks
R=benchmarks/results/runs
mkdir -p "$R"

wait_for() {  # unit: block until a systemd user unit is no longer active
  while systemctl --user is-active --quiet "$1"; do sleep 300; done
}

if [ "${1:-}" != score ]; then
  bash scripts/fetch_data.sh            # HG002 truths and reference, checked
  bash scripts/fetch_hg005.sh           # HG005 truths, checked
  bash scripts/genome_run.sh            # HiFi run (docs/design/genome-run-plan.md)
  wait_for phasebudget-genome
  bash scripts/genome_run.sh ont chr20 chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
    chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr21 chr22
  wait_for phasebudget-genome-ont
  "$PY" benchmarks/trio_run.py --chromosomes chr20 chr16 chr18 chr21 chr22
  "$PY" benchmarks/tool_run.py --jobs 3   # tool-level test (docs/design/tool-run-plan.md)
  PHASEBUDGET_PLATFORM=hg005 "$PY" benchmarks/genome_run.py --jobs 4 \
    --chromosomes chr20 $(printf 'chr%s ' 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 21 22)
fi

# HiFi: registered analysis, the v5.0q re-score, mechanism, selector fix, calibration.
"$PY" benchmarks/score_genome.py --jobs 6 --tsv "$R/hifi_per_chromosome.tsv" > "$R/hifi_score.txt"
"$PY" benchmarks/score_genome.py --jobs 6 --truth v5 --tsv "$R/hifi_per_chromosome_v5.tsv" \
  > "$R/hifi_score_v5.txt"
"$PY" benchmarks/benchmark_spanning.py --count informative --depths full 10 --caps 3 5 \
  > "$R/hifi_spanning_informative.txt"
"$PY" benchmarks/benchmark_selector_fix.py > "$R/hifi_selector_fix.txt"
"$PY" benchmarks/calibrate_posterior.py --fit chr1 chr3 chr5 chr7 chr9 chr11 chr13 chr15 chr17 \
  chr19 chr21 --test chr2 chr4 chr6 chr8 chr10 chr12 chr14 chr16 chr18 chr22 --depths full 10 \
  > "$R/hifi_calibration.txt"

# ONT: e on chromosome 20, the registered prediction (reads no truth), then the scoring.
PHASEBUDGET_PLATFORM=ont "$PY" benchmarks/measure_allele_error.py chr20 > "$R/ont_allele_error.txt"
PHASEBUDGET_PLATFORM=ont "$PY" benchmarks/predict_regime.py --error 0.00503 \
  --tsv "$R/ont_prediction.tsv" > "$R/ont_prediction.txt"
PHASEBUDGET_PLATFORM=ont "$PY" benchmarks/score_genome.py --jobs 6 \
  --tsv "$R/ont_per_chromosome.tsv" > "$R/ont_score.txt"
python3 benchmarks/score_prediction.py benchmarks/results/runs/2026-10-01_ont_prediction_baselines.tsv \
  "$R/ont_score.txt" "$R/ont_prediction.tsv" > "$R/ont_prediction_score.txt"
PHASEBUDGET_PLATFORM=ont "$PY" benchmarks/benchmark_spanning.py --count informative \
  --depths full 10 --caps 3 5 > "$R/ont_spanning_informative.txt"
PHASEBUDGET_PLATFORM=ont "$PY" benchmarks/calibrate_posterior.py --T 0.7264 --eps 0 \
  --test $(printf 'chr%s ' 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 21 22) \
  --depths full 30 20 10 > "$R/ont_calibration.txt"

# Trio, the tool-level test, and HG005 (e and the prediction before the scoring).
"$PY" benchmarks/trio_run.py --score > "$R/trio_score.txt"
"$PY" benchmarks/tool_run.py --score > "$R/tool_score.txt"
PHASEBUDGET_PLATFORM=hg005 "$PY" benchmarks/measure_allele_error.py chr20 > "$R/hg005_allele_error.txt"
PHASEBUDGET_PLATFORM=hg005 "$PY" benchmarks/predict_regime.py --error 0.00011 \
  --tsv "$R/hg005_prediction.tsv" > "$R/hg005_prediction.txt"
PHASEBUDGET_PLATFORM=hg005 "$PY" benchmarks/score_genome.py --jobs 6 \
  --tsv "$R/hg005_per_chromosome.tsv" > "$R/hg005_score.txt"
python3 benchmarks/score_prediction.py --hg005 - "$R/hg005_score.txt" "$R/hg005_prediction.tsv" \
  > "$R/hg005_prediction_score.txt"

(cd docs/submission && latexmk -pdf paper.tex)
echo "paper regenerated"
