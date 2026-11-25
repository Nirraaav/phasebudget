#!/usr/bin/env bash
# Stop the genome-wide run and the trio benchmark. Finished jobs stay finished;
# running jobs are stopped and will be redone on the next start.
systemctl --user stop phasebudget-genome phasebudget-genome-ont phasebudget-tools phasebudget-trio 2>/dev/null
# Jobs run in their own memory-limited scopes, so stop them explicitly too.
pkill -f 'benchmarks/genome_ru[n].py --task' 2>/dev/null
pkill -f 'benchmarks/hg002_pipeline[s].py' 2>/dev/null
pkill -f 'benchmarks/extract_fragment[s].py' 2>/dev/null
pkill -f 'whatshap phas[e]' 2>/dev/null
pkill -f 'extractHAIR[S]|HAPCUT[2] ' 2>/dev/null
pkill -f 'longphas[e] phase' 2>/dev/null
pkill -f 'hiphas[e] --bam' 2>/dev/null
pkill -f 'whatshap_exac[t].py' 2>/dev/null
pkill -f 'benchmark_trio_budge[t].py' 2>/dev/null
echo "stopped"
