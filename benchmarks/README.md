# Benchmarks

Each script prints its parameters in the header and takes them on the command
line. Run outputs are in `benchmarks/results/runs/`; each starts with its
command, parameters and commit. Figures are in `benchmarks/results/`.

The sections below cover the read-selection scripts. The certificate tests
are one script per registration: each `<id>_*.py` in the `q`, `r`, `w`,
`k`, `l` and `z` series, and `stratified_*`, `third_score`,
`confidence_run` and `confidence_score`, implements the plan in
`docs/design/<id>-*.md`, and its usage is in its docstring.
`phase_confidence.py` scores any phaser's VCF from a BAM.
`scripts/reproduce_paper.sh` runs the read-selection scripts and the
`plot_*` scripts in order.

## `benchmark_heuristic_gap.py`

The headline measurement: every selector's additive objective divided by the
certified optimum on the same instance. Four levels are separate columns,
objective ratio, selection time apart from solve time, wMEC cost, and switch
error, and each selector's own objective gets a column so no method is scored
only on another's.

```bash
PYTHONPATH=src python3 benchmarks/benchmark_heuristic_gap.py \
    --variants 60 --coverage 20 --span 6 --replicates 25 \
    --caps 3 5 8 --error-rates 0.01 0.05 0.10 0.20 --quality-spread 2.0
```

Two columns need their labels read. `MEC kept` is the cost on the reads a
selector retained and is not comparable across selectors, since keeping fewer
reads lowers it; `MEC full` scores the same haplotype against every read and is
the comparable one. Run it with `--quality-spread 0` and the selection problem
stops depending on the error rate entirely, because every base then carries the
same weight and a read's score is its called-allele count. That control is in
The development notebook and it is the reason the flag exists.

Two Fischer and Marschall rows appear. `fischer_marschall` follows the
paper's pseudocode and `whatshap` is `select_whatshap`, the port of the code
WhatsHap 2.8 ships; they choose different reads, and `whatshap` is the
deployed baseline. `flips`, `N50` and `QAN50` follow `docs/design/metrics.md`.

## `benchmark_selection.py`

Every selector scored on additive weight, minimum coverage, covered variants,
connected components, read count, and runtime. The Fischer and Marschall
selector, the maximum-minimum coverage selector, and the additive min-cost-flow
selector have different objectives. This script does not call any difference
between them an optimality gap. It says nothing about phasing accuracy, which
`benchmark_phasing.py` measures.

```bash
PYTHONPATH=src python3 benchmarks/benchmark_selection.py --replicates 8
```

## `benchmark_phasing.py`

End-to-end accuracy under random, greedy, Fischer--Marschall,
maximum-minimum-coverage, additive exact, and adaptive selection at a fixed
coverage cap. Reports switch error, wMEC cost against the full read set, the
additive selection objective under the base weights, and wall time.

```bash
PYTHONPATH=src python3 benchmarks/benchmark_phasing.py --replicates 12
```

The adaptive strategy reweights reads, so its own objective value is not
comparable with the others. The script rescores every selection under the base
weights before printing.

## `benchmark_objective_truth.py`

Whether the wMEC optimum coincides with the true phasing, and whether MEC
margins identify the boundaries that carry switch errors. This is the project's
central hypothesis and the two columns are meant to be read against each other.

```bash
PYTHONPATH=src python3 benchmarks/benchmark_objective_truth.py \
    --variants 30 --coverage 6 --span 5 --replicates 25
```

Runs at full coverage, with no selection, so the dynamic program sees every
read. Coverage above about 12 makes it slow; that is the constraint the rest of
the project exists to work around.

`--coverages` and `--spans` accept grids. The legacy singular flags retain a
one-cell command line. The simulator starts with heterozygous variant sites,
so this benchmark varies coverage, span, and error rate only. It does not
pretend that changing the density of those sites is a heterozygosity sweep.

## `benchmark_margin_confound.py`

Whether the adaptive loop's margin signal is the fragility of the phasing or
the thinness of the subset selection just chose, which is OQ-10. Compares
subset margins against a full-matrix oracle on agreement, informativeness, and
switch error, and reports a matched pair at one reweight each so the round
budget cannot explain a difference between the two sources.

```bash
PYTHONPATH=src python3 benchmarks/benchmark_margin_confound.py --replicates 60
```

The oracle needs the whole matrix under the dynamic program's state limit, so
every instance it runs on is one that did not need selection. Nothing it
reports transfers to the coverage-constrained regime on its own.

## `benchmark_heterozygosity.py`

The objective-truth gap and margin informativeness along the variant-density
axis, using `simulate_genomic`. This is the axis
`benchmark_objective_truth.py` records as impossible for it to vary.

```bash
PYTHONPATH=src python3 benchmarks/benchmark_heterozygosity.py \
    --variants 100 --read-length 5000 --replicates 25
```

Connectivity is an outcome here, not a control: lowering the density is what
breaks a phase block. Blocks, phased fraction, and comparable pairs are
reported as columns, and the switch rate across density rows is a comparison
over different denominators. Keep `read_length` well below the genome span
`n_variants / heterozygosity`, or reads overhanging the ends dominate the
variants-per-read mean.

## `benchmark_sensitivity.py`

Whether the adaptive gain is a result or one setting of boost, quantile, and
iteration cap, which is OQ-6. Reports the fraction of the grid that beats
static selection and the worst setting in it. The best setting of a swept
parameter is a selection effect and the summary does not lead with it.

```bash
PYTHONPATH=src python3 benchmarks/benchmark_sensitivity.py --replicates 12 --full-grid
```

The static baseline is computed once per instance and shared by every cell, so
a difference between two cells is never a difference in instances.

## `benchmark_genotype_errors.py`

What the complementary-haplotype restriction costs when some sites called
heterozygous are not, which is OQ-7. Scores the complementary dynamic program
and the general model on read-bipartition error and on phase relations between
truly heterozygous sites. The general model is solved by `phasebudget.general`,
so it runs at the scale of the other benchmarks; `--solver bruteforce` uses the
`2^reads` oracle instead.

```bash
PYTHONPATH=src python3 benchmarks/benchmark_genotype_errors.py \
    --variants 60 --coverage 8 --span 5 --replicates 40 \
    --homozygous-rates 0.0 --error-rates 0.02 0.10     # one run per rate: 0.0 0.05 0.1 0.2 0.3
```

The general model may give both parts the same consensus at a variant, which is
it declining to phase that variant. Two switch columns are reported for it, one
charging every declined pair and one scoring only the pairs it phased, next to
the fraction it phased. Reading either alone gives the wrong answer.

## `benchmark_majidian_bound.py`

Checks Theorem 2 of Majidian, Kahaei and de Ridder 2020 against exact
unweighted MEC. The header prints a fixed counterexample in which no column is
half wrong and the truth still loses, then one row per coverage and error rate:
their formula, the event it computes, the fraction of instances whose truth is
an optimum, and the fraction where the two disagree.

```bash
PYTHONPATH=src python3 benchmarks/benchmark_majidian_bound.py \
    --variants 30 --span 5 --coverages 4 6 8 \
    --error-rates 0.05 0.10 0.20 --replicates 100 --seed-base 9000
```

## `benchmark_selection_at_scale.py`

Milestone 2: the heuristic gap on chromosome-scale simulated stretches from the
genomic simulator, with read profiles `hifi` (15 kb, error 0.01), `ont` (10 kb,
0.05) and `noisy` (10 kb, 0.10). Includes HapCHAT's merge-then-select and
`merge_optimal`, which merges and then selects by min-cost flow; those two have
no objective ratio because their reads are merged.

`--systematic-rate`, `--systematic-strength`, `--burst-rate`, `--burst-length`
and `--burst-error` switch on errors the base qualities do not report, and the
`mixed` column gives the share of merged reads holding both haplotypes.

```bash
PYTHONPATH=src python3 benchmarks/benchmark_selection_at_scale.py \
    --variants 2000 --replicates 10 --processes 16
PYTHONPATH=src python3 benchmarks/benchmark_selection_at_scale.py \
    --variants 2000 --replicates 10 --caps 3 5 10 --profiles ont noisy --processes 16
```

## `benchmark_phase_diagram.py` and `plot_phase_diagram.py`

Milestone 1: one grid over depth, error rate, variant density and base-quality
spread, every read to the exact solver, writing
`benchmarks/results/phase_diagram.csv` with its provenance in comment lines.
Instances above the state limit are skipped, redrawn and counted. The plot
script draws `phase_diagram.png` and `oq5_informativeness.png` and prints the
boundary table and the OQ-5 correlations; it needs matplotlib, which the
library does not.

```bash
PYTHONPATH=src python3 benchmarks/benchmark_phase_diagram.py --replicates 30
python3 benchmarks/plot_phase_diagram.py
```

## `benchmark_adaptive_signal.py`

Exact static selection, the loop steered by MEC margins and the loop steered by
posterior margins, paired per instance on 60-variant instances.

```bash
PYTHONPATH=src python3 benchmarks/benchmark_adaptive_signal.py --replicates 100
```

## `benchmark_adaptive_at_scale.py`

The adaptive loop paired against exact static selection on chromosome-scale
instances. `--signal posterior --error-probability p` uses the absolute
fragile rule, which marks a boundary fragile when the posterior says its call
is wrong with probability above `p`.

```bash
PYTHONPATH=src python3 benchmarks/benchmark_adaptive_at_scale.py --processes 16
PYTHONPATH=src python3 benchmarks/benchmark_adaptive_at_scale.py --processes 16 \
    --signal posterior --error-probability 0.05
PYTHONPATH=src python3 benchmarks/benchmark_adaptive_at_scale.py --processes 16 \
    --signal posterior --error-probability 0.05 --acceptance stitched --window 50 \
    --profiles ont noisy --caps 3 5 10 --replicates 32 --control
```

`--seed-base` moves the instances, for held-out checks. `--acceptance
stitched` accepts rounds window by window, default window 25
(`docs/design/adaptive-acceptance.md`); `--control` adds a third arm steered
by the posterior margins shuffled across boundaries and prints the signal
paired against it, which is how OQ-13 was answered. `--iterations`, `--boost`
and `--window` are the milestone-4 sensitivity axes, and `--merge` runs both
arms on reads merged by HapCHAT's rule, so its static arm is `merge_optimal`.

## `benchmark_convention_cost.py`

`OPT_span / OPT_called`: the weight the span convention gives up against the
called convention, whose optimum comes from a 0-1 program under CBC, and the
span coverage the called optimum would impose on the dynamic program.

```bash
PYTHONPATH=src python3 benchmarks/benchmark_convention_cost.py
```

## `build_trio_truth.py`

Phases HG002's heterozygous SNVs on one chromosome by transmission from the
GIAB v4.2.1 calls of HG003 and HG004 (`docs/design/truth-sets.md`), prints
where every candidate site went, scores GIAB's two phased files against it,
and writes the truth as a VCF. The inputs are chromosome-20 extracts in
`data/giab/`, which is not tracked; the development notebook records the download.

```bash
PYTHONPATH=src python3 benchmarks/build_trio_truth.py \
    --out data/giab/HG002_chr20_trio_truth.vcf.gz
```

## `whatshap_bridge.py`

Not a benchmark: converts a `FragmentMatrix` into a WhatsHap `ReadSet` and
runs WhatsHap's own dynamic program and read selection on it, for
`tests/test_whatshap_port.py` and for regenerating that test's fixture.
WhatsHap is an optional extra, pinned to the release the port reproduces;
the `whatshap` CI job installs it so the live comparisons run.

```bash
python3 -m venv .venv
.venv/bin/pip install -e ".[dev,ilp,whatshap]"
.venv/bin/python -m pytest tests/test_whatshap_port.py
```

## `diagnose_whatshap_score_rule.py`

Runs the WhatsHap port as shipped and with the paper's score update
substituted, nothing else changed, and reports objective ratio, block count,
switch error, N50 and QAN50 for each. It isolates the one line of WhatsHap
that fragments phase blocks at small caps.

```bash
PYTHONPATH=src python3 benchmarks/diagnose_whatshap_score_rule.py \
    --profiles ont noisy --caps 3 5 10 15 --replicates 10
```

## `benchmark_resources.py`

Runtime and peak memory of `select_optimal`, `select_whatshap`, the native DP
with and without checkpointing, with and without the backward pass, and the
Python DP, against chromosome length and cap. Each step runs in a fresh
process and peak memory is the rise in `VmHWM`, so it needs Linux. The
unchecked native pass is skipped when its table would exceed `--memory-gb`.

```bash
PYTHONPATH=src python3 benchmarks/benchmark_resources.py \
    --variants 500 2000 8000 --caps 5 10 15 20 --memory-gb 6
```

## `extract_fragments.py` and `describe_fragments.py`

The fragment step of the real-data plan: heterozygous SNVs of one sample
inside a BED, alleles read by WhatsHap's `PhasedInputReader` as `whatshap
phase` reads them, written as gzipped JSON. `phred_to_log_odds` converts
WhatsHap's phred qualities to the library's weights (amendment 3 of the
plan). `describe_fragments.py` prints coverage, read content, qualities and
whether each cap fits, which is the OQ-9 check.

```bash
.venv/bin/python benchmarks/extract_fragments.py --bam data/hg002/HG002_chr20_ccs11kb.bam \
    --vcf data/giab/HG002_chr20_v4.2.1.vcf.gz --bed data/giab/HG002_chr20_v4.2.1.bed \
    --reference data/hg002/chr20.fa --out data/hg002/HG002_chr20_fragments.json.gz
PYTHONPATH=src .venv/bin/python benchmarks/describe_fragments.py data/hg002/HG002_chr20_fragments.json.gz
```

## `hg002_pipelines.py`, `score_hg002.py`, `subsample_bam.py`, `compare_whatshap_reads.py`

The registered real-data run. `subsample_bam.py` draws the seed-0 half-depth
read set. `hg002_pipelines.py` runs one registered pipeline at one cap on a
fragment file in its own process and writes the phase call, selection,
objective, wMEC, time and peak memory. `score_hg002.py` scores every run and
the external tools' outputs against the transmission truth, applies the
registered decision rules with the 1 Mb window bootstrap, and prints the
hifiasm table and the switch sites most phasings share; `--exclude` is a
labelled post-hoc truth mask. `compare_whatshap_reads.py` checks the port's
selection against `whatshap phase --output-read-list`.

```bash
PYTHONPATH=src .venv/bin/python benchmarks/hg002_pipelines.py \
    --fragments data/hg002/HG002_chr20_fragments.json.gz --pipeline optimal --cap 5 \
    --weights logodds --signals --out data/hg002/runs/full_optimal_c5_logodds.json.gz
PYTHONPATH=src .venv/bin/python benchmarks/score_hg002.py
```

## `benchmark_pedigree.py`

A simulated trio phased by WhatsHap's `PedigreeDPTable`, each member's reads
selected at the per-member cap a trio gets (`cap // 3`), against the child
phased alone at the same cap. Needs the `whatshap` extra.

```bash
PYTHONPATH=src .venv/bin/python benchmarks/benchmark_pedigree.py \
    --profiles ont noisy --caps 3 5 --sites 2000 --replicates 16
```

`benchmark_adaptive_at_scale.py` takes `--merge-p` and `--merge-threshold`
for HapCHAT's merge rule and prints block N50 and QAN50 for both arms.

## `benchmark_trio_budget.py`

A simulated trio phased by WhatsHap's pedigree DP under three ways of
spending a total coverage budget: WhatsHap's per-member thirds, one shared
exact selection, and a shared selection weighted towards
transmission-ambiguous sites. Needs the `whatshap` extra.

```bash
PYTHONPATH=src .venv/bin/python benchmarks/benchmark_trio_budget.py --processes 4
```

## `genome_run.py`, `genome_common.py`, `recalibrate.py` and `scripts/`

The registered genome-wide run (`docs/design/genome-run-plan.md`).
`genome_run.py` schedules every job:
- streaming each chromosome's reads, subsampling, both fragment extractions,
  `whatshap phase`, HapCUT2, deleting the reads, the library pipelines,
  the recalibration fit and the base-quality pipelines;
- each job under a memory ceiling, at low priority, with memory and disk
  guards;
- resuming from `data/genome/done/`.

`recalibrate.py` fits the base-quality table on the odd chromosomes. Start,
follow and stop the run with the scripts, which run it as systemd user
services independent of any terminal:

```bash
bash scripts/genome_run.sh
bash scripts/genome_status.sh
bash scripts/genome_stop.sh
```

## Recording a result

A benchmark result must record the commit, the command line, the machine, the
Python version, and the seed policy. Simulated results are not evidence about
real data and are labelled as simulated wherever they appear. Every run file
in `benchmarks/results/runs/` starts with the command that produced it.
