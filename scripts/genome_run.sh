#!/usr/bin/env bash
# Start, or resume, the registered genome-wide HG002 run (docs/design/genome-run-plan.md)
# and the trio budget benchmark, detached from any terminal.
#
# Both run as systemd user services, so closing the terminal or logging out
# does not stop them. Finished jobs are remembered in data/genome/done/, so running this
# script again after an interruption continues where it stopped.
#
#   bash scripts/genome_run.sh        start or resume
#   bash scripts/genome_status.sh     progress
#   bash scripts/genome_stop.sh       stop everything cleanly
#
# `bash scripts/genome_run.sh ont [chromosomes...]` starts the registered ONT run
# (docs/design/ont-run-plan.md) instead, as its own service; the other scripts
# take `ont` the same way.
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT="$PWD"
mkdir -p tmp data/genome

inhibit=()
# Keep the laptop from suspending while a run is going, if the system allows it.
if systemd-inhibit --what=sleep:idle --who=phasebudget --why=test true 2>/dev/null; then
  inhibit=(systemd-inhibit --what=sleep:idle --who=phasebudget --why="genome-wide phasing run")
fi

if [ "${1:-}" = ont ]; then
  shift
  chromosomes=()
  [ $# -gt 0 ] && chromosomes=(--chromosomes "$@")
  if systemctl --user is-active --quiet phasebudget-genome-ont; then
    echo "ONT run: already running"
  else
    systemd-run --user --unit=phasebudget-genome-ont --collect \
      -p WorkingDirectory="$ROOT" --setenv=PYTHONPATH="$ROOT/src:$ROOT/benchmarks" \
      --setenv=PHASEBUDGET_PLATFORM=ont \
      -p StandardOutput=append:"$ROOT/tmp/genome_ont_run.log" \
      -p StandardError=append:"$ROOT/tmp/genome_ont_run.log" \
      "${inhibit[@]}" nice -n 10 "$ROOT/.venv/bin/python" "$ROOT/benchmarks/genome_run.py" --jobs 4 \
      "${chromosomes[@]}" >/dev/null
    echo "ONT run: started (log tmp/genome_ont_run.log)"
  fi
  exit 0
fi


if systemctl --user is-active --quiet phasebudget-genome; then
  echo "genome run: already running"
else
  systemd-run --user --unit=phasebudget-genome --collect \
    -p WorkingDirectory="$ROOT" --setenv=PYTHONPATH="$ROOT/src:$ROOT/benchmarks" \
    -p StandardOutput=append:"$ROOT/tmp/genome_run.log" -p StandardError=append:"$ROOT/tmp/genome_run.log" \
    "${inhibit[@]}" nice -n 10 "$ROOT/.venv/bin/python" "$ROOT/benchmarks/genome_run.py" --jobs 5 \
    >/dev/null
  echo "genome run: started (log tmp/genome_run.log)"
fi

if grep -q T3DONE tmp/t3_trio_budget.txt 2>/dev/null; then
  echo "trio budget benchmark: already finished"
elif systemctl --user is-active --quiet phasebudget-trio; then
  echo "trio budget benchmark: already running"
else
  systemd-run --user --unit=phasebudget-trio --collect \
    -p WorkingDirectory="$ROOT" --setenv=PYTHONPATH="$ROOT/src" \
    bash -c "nice -n 19 ionice -c 3 '$ROOT/.venv/bin/python' benchmarks/benchmark_trio_budget.py --processes 4 > tmp/t3_trio_budget.txt 2>&1 && echo T3DONE >> tmp/t3_trio_budget.txt"
  echo "trio budget benchmark: started (output tmp/t3_trio_budget.txt)"
fi
