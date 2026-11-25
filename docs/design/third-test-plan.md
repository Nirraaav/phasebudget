# Pre-registered third test: per-group error control on an unrelated genome, one basecaller

Written on 13 October 2026. The commit that adds this file is the
registration; its push time is the third-party record. Creating this file is
what lets `benchmarks/genome_common.py` open the test reads (`hg001_v5`).

## Why

The second test (`docs/design/stratified-run-plan.md`, scored 13 October)
supported its primary hypotheses, with two gaps a reader will see. Its six
Pedigree test genomes are relatives of the calibration genome, from the same
project; its R10 test genome, HG005, has no phased indel truth and was
sequenced with a newer basecaller than its calibration genome. No test put an
unrelated genome, sequenced with the same basecaller as its calibration
genome, against a truth that phases indels. This test does.

## Material seen before writing

Everything in the development notebook up to 13 October, in particular the second test's
results. HG001 (NA12878) is not new to this project:

- its 2023 ONT sup reads (`giab_2023.05`, flowcells PAO99309 and PAO98424)
  were the first registered test's HG001 genome (scored 4 October), and were
  used after it in development: track C, gate A within HG001, and the
  cross-genome comparison of 8 October (`stratified_gate.py --cross
  --calibration cal_r10`). There, thresholds calibrated on `cal_r10`'s SNV
  pairs held on HG001's SNV pairs against the Platinum Pedigree truth in
  every stratum, cell and beta; per-group thresholds trusted 0.830 to 0.952
  of pairs at `beta = 1e-3`, and the pooled threshold trusted none for
  whatshap and 0.617 to 0.967 for LongPhase (0.967 at 30x, above per-group's
  0.952). Indel pairs were not part of it;
- the Platinum Pedigree truth for NA12878 is the second test's Pedigree
  calibration truth, and NA12878's Platinum Pedigree ONT reads were its
  calibration genome.

What is new is the reads: `giab_2025.01` sup flowcells PAW79146 and PAW81754
for HG001, and PAW70337 and PAW71238 for HG002, none of them used before. On
13 October the headers of these four BAMs (and of HG005's) were read: each has
one read group, basecall model `dna_r10.4.1_e8.2_400bps_sup@v5.0.0`, and
GRCh38 contigs named `chr1` onwards. The flowcell lists of both releases show
that the 2025 flowcells are different sequencing runs from the 2023 ones, not
re-basecalled copies. No HG001 2025 read has been read.

The calibration run on `cal_v5` started at 12:44 on 13 October, before this
commit; calibration genomes are not test data. Its results had not been
looked at when this file was written.

## What is tested

The procedure is unchanged from the second test: Learn then Test calibrated
within each group, Vovk's conditional construction applied to Learn then
Test's risk, not claimed as new. The question is whether its guarantee
carries over from HG002 to an unrelated genome sequenced the same way.

## Frozen design

Everything in the second test's "Frozen design" holds unchanged: the score,
the groups (`indel`, `segdup`, `lowmap`, `tandem`, `difficult`, `easy`), the
100 kb unit, the grid, `delta = 0.05`, `beta` in {3e-4, 1e-3, 3e-3}, the
population (consecutive pairs of phased variants that the truth phases, with
the truth's own heterozygous calls, phase removed, as phaser input), the
depths 30, 20 and 10 drawn as there, both phasers with the same commands
(LongPhase with `--indels`), and the 21 confirmatory autosomes.

| | calibration (`cal_v5`) | test (`hg001_v5`) |
|---|---|---|
| genome | HG002 | HG001 (NA12878), unrelated to HG002 |
| reads | `giab_2025.01/basecalling/sup/HG002/{PAW70337,PAW71238}/calls.sorted.bam` | `giab_2025.01/basecalling/sup/HG001/{PAW79146,PAW81754}/calls.sorted.bam` |
| basecaller | R10.4.1 sup v5.0.0 | the same |
| phaser input and truth | GIAB v5.0q heterozygous SNVs and indels in its BED | Platinum Pedigree truthset, NA12878, inside `hq_regions_final.bed.gz` |

Both flowcells of a genome are drawn from with the same fraction and merged,
as for HG005 in the second test.

- **Shift stated now.** Calibration and test differ in person and in how the
  truth was built (assembly-based v5.0q against inheritance-based Platinum
  Pedigree), and they share a basecaller, a release and a lab.

## Hypotheses

A cell is (depth, phaser, beta, group). In a cell, `n` test units hold a
trusted pair and `k` hold a wrong one; it violates when the exact binomial
`P(X >= k | n, beta)` is below the hypothesis's level. A group that certifies
nothing cannot violate and is reported as abstaining. **Levels are fixed by
the cell counts registered here**, whatever number of cells turns out to have
test pairs (the second test's S3 was registered with a count its data could
not reach).

- **V1, the easy group carries over (primary).** 3 depths x 2 phasers x 3
  betas = 18 cells, level 0.05/18. *Supported* if no cell violates.
- **V2, the difficult and indel groups carry over (primary).** Groups
  `segdup`, `lowmap`, `tandem`, `difficult`, `indel`: 90 cells, level 0.05/90.
  *Supported* if no cell violates. The scorer reports how many of the 90
  certify something, since a group that certifies nothing passes by
  abstaining.
- **V3, usefulness (primary).** Group `easy`, `beta = 1e-3`, depths 30 and 20,
  both phasers: 4 cells, with the second test's S4 rule (at least 99 per cent
  of the group's pairs trusted, and wrong trusted pairs at most half the
  wrong pairs when there are at least 20, else under 1 per 10,000 trusted).
  *Supported* if 4 succeed, *partly* if 3.
- **V4, per-group against pooled (primary).** Depths 30 and 20, both phasers,
  `beta = 1e-3` (4 cells): per-group thresholds trust a larger share of all
  test pairs than one pooled threshold calibrated on the same rows.
  *Supported* if 4, *partly* if 3.
- **V5 (descriptive).** Per depth and phaser, the share of pairs in groups
  that certify nothing, at `beta = 1e-3`.
- **V6, whatshap's indel errors (secondary).** At 30x, on indel pairs both
  phasers phase, pairs only whatshap gets wrong outnumber pairs only
  LongPhase gets wrong (one-sided binomial, level 0.05).
- **V7 (descriptive).** The cells of V1 and V2 with thresholds from
  `cal_r10` (HG002, the 2023 sup reads, the second test's R10 calibration),
  to show what a basecaller change between calibration and test does on the
  same test reads. No decision is drawn from it.

## Software

- `benchmarks/confidence_run.py` and `benchmarks/genome_common.py` at the
  registration commit: platforms `cal_v5` and `hg001_v5`, three chromosomes
  at a time per genome, disk reserved per chromosome, reads deleted after
  phasing (`tmp/t3_runs.sh`, systemd user units `phasebudget-t3-cal_v5` and
  `phasebudget-t3-hg001_v5`).
- `benchmarks/third_score.py` at the registration commit decides V1 to V4
  and V6 and reports V5 and V7, reading only stored fragments, phaser VCFs
  and truths. It imports the second test's scoring functions from
  `benchmarks/stratified_score.py` unchanged. Its `--dry-run` was run on
  `cal_r10` (odd chromosomes calibrate, even ones stand in for the test) before
  this commit, as a check of the code
  (`benchmarks/results/runs/2026-10-13_third_score_dryrun.txt`, 4 min,
  2.5 GB): V1 18 of 18 cells certifying, 0 violating; V2 31 of 90
  certifying, 0 violating, smallest p 0.123; V3 and V4 4 of 4; abstention
  at `beta = 1e-3` 0.40 to 0.41 (whatshap) and 0.24 to 0.34 (LongPhase) on
  half of HG002; V6 2,045 against 136.
- Before scoring, `scripts/check_stratified_inputs.py`, extended to the two
  new platforms, checks every input exists and decodes, with no truth
  compared with any phasing.
- Every threshold is scored at full precision from the run file, never from
  printed text. Every deviation is a dated amendment below, stating whether
  results had been seen.

## Limits stated now

- One test genome. A supported V1 or V2 says the bound carried over to HG001
  in these conditions; it does not estimate how often it carries over.
- HG001 has been studied in this project, on other reads, and the 8 October
  development comparison already found the 2023-read analogue of V1 and of
  V2 without the indel group holding. This test adds the indel group, new
  reads and a newer shared basecaller under a registered rule; it is not an
  unseen person.
- The 18 and 90 cells are correlated (nested depths, shared reads between the
  phasers); the Bonferroni levels stay valid under that dependence and are
  conservative.

## Amendments

**1. The runs move to a KCDH server (14 October 2026; no results seen).**
The laptop runs stop whenever the laptop shuts down. On 13 October both
were paused at 13:45 with chromosomes 21 and 22 finished and 17 to 19 partly
done; no phasing had been compared with a truth. The runs restart from the
beginning on the KCDH compute server (Ubuntu 22.04, 96 cores, 1 TB memory), with
the code at this amendment's commit, Python 3.12.3, the package set of
`requirements-lock.txt` (checked identical with `pip freeze`; WhatsHap 2.8,
pysam 0.24.1) and the same LongPhase 2.0.2 binary. Its `data/` is on the
lab's network share, holding copies of the reference, the v5.0q files and
the Platinum Pedigree truth. The run script is unchanged and runs without
systemd. Scoring uses the server's outputs. The
laptop's finished chromosomes (21 and 22 of both genomes) are kept and
compared with the server's before scoring: hetsnv VCFs, fragments, and the
phaser VCFs decompressed and without `##` header lines, which hold paths and
dates. Any difference is reported with the scores.

**2. More chromosomes at once (14 October 2026; no results seen).** At the
author's request the server runs use more of the machine: `--workers 6`
per genome instead of 3, and at most 8 chromosomes phasing at once instead
of 4, a slot past the second only while 64 GB of memory is available.
`confidence_run.py` reads the last two from `PHASEBUDGET_PHASE_SLOTS` and
`PHASEBUDGET_SLOT_FREE_GB`, with its old values as defaults. Neither setting
changes an output: each draw has a fixed seed and does not depend on
which process makes it, and both phasers run single-threaded. The run is
stopped and restarted from its step markers; a step cut off is redone.

**3. Ten chromosomes at once (14 October 2026; no results seen).** One
remote draw ran at 1.6 to 1.9 MB/s beside the run's 12 others, which
together read 24 to 26 MB/s, so the limit is per draw; at the author's
request `--workers 10` per genome replaces 6; the 8 phasing slots
and the 64 GB rule stay. The reasons of amendment 2 for why outputs do not
change apply unchanged. The run is stopped and restarted from its markers.

*Correction, 14 October, 15:55; no results seen.* The reason given above
was wrong. With the run's 20 draws reading 25.7 MB/s in total, an extra
`curl` stream from S3 got 1.07 MB/s, the same share as each draw: the limit
is the server's link to S3 (about 25 MB/s at that hour), not each draw. Ten
workers do not raise the total. They stay, since another restart would lose
the draws in progress, and no output depends on them.
