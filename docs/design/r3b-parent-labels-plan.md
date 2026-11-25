# Pre-registered R3b: certificates calibrated on labels from HG005's parents

Written on 27 October 2026. The commit that adds this file is the
registration; its push time is the third-party record. The author asked for
it on 27 October.

## Why

R3 (`docs/design/r3-label-eligible.md`, development, 16 October) calibrated
per-group thresholds on HG002's own odd chromosomes, with phase labels from
ONT Clair3 calls on both flowcell replicates of HG003 and HG004, and scored
them on HG002's even chromosomes. On the pairs the labels can reach it had 0
violating cells of 108, and `easy` at 1e-3 certified in all six depth and
phaser cells, trusting 0.774 of all `easy` pairs. A procedure that needs only
the parents' reads from the same sequencing batch, and no benchmark truth, is
the answer to R2's floor: a new batch needs its own labels at genome scale.
R3b repeats R3 on a second trio.

## Material seen before writing

Everything in the development notebook up to 27 October. HG005's ONT reads (`hg005_ont`,
giab_2025.01 sup, both flowcells) were phased and scored against GIAB's trio
truth in the second test (X6, 13 October), so the phasers' errors on HG005
have been seen in aggregate. The parents' calls (HG006, HG007) have never
been fetched. No label, eligibility or label-calibrated threshold of HG005
has been computed.

## Design

Unchanged from R3 unless listed: the `agree` filter (a site is labelled
when both replicates give it the same paternal allele), eligibility (both
sites labelled), the calibration (eligible odd-chromosome pairs, wrong
against the labels), the groups, score, 100 kb unit, grid,
`learn_then_test_by_group`, `delta = 0.05`, `beta` in {3e-4, 1e-3, 3e-3},
depths 30, 20 and 10, whatshap and LongPhase as stored.

- **Child.** HG005, `hg005_ont`, the X6 phasings.
- **Labels.** ONT `wf-human-variation` sup calls (`SAMPLE.wf_snp.vcf.gz`)
  from `giab_2025.01/analysis/wf-human-variation/sup/`: replicate A is
  HG006 `PAY77227` with HG007 `PAY12990`, replicate B is HG006 `PBA16846`
  with HG007 `PBA20413` (flowcells in listing order). Fetched after this
  commit, md5 recorded.
- **Truth.** GIAB v4.2.1 HG005 trio transmission
  (`genome_common.primary_truth`), used only to score.
- **Cells.** (depth, phaser, `beta`, group), scored on eligible
  even-chromosome pairs: 3 x 2 x 3 x 6 = 108 cells; level 0.05/108 =
  4.63e-4 whatever number have pairs.

## Hypotheses

- **R3b-1, validity (primary).** No cell violates: in no cell is the exact
  `P(X >= k | n, beta)` below 4.63e-4, `n` the eligible test units with a
  trusted pair and `k` those with a wrong one.
- **R3b-2, usefulness (primary).** `easy` at `beta = 1e-3` certifies in all
  six depth and phaser cells, and the trusted share of all even-chromosome
  `easy` pairs, abstaining on every ineligible pair, is at least 0.70.
- **Descriptive.** R3's other outputs: the `all` and `double` filters,
  scoring on all even pairs, and the share and error of eligible,
  uninformative and other pairs per group.

## Risks stated now

- The truth comes from the same family's transmission, called from
  Illumina reads. A parental genotype error that the ONT calls and GIAB's
  calls share would be wrong in label and truth alike and invisible here.
- The cached rows hold only pairs the truth phases, so calibration runs on
  truth-phased eligible pairs, as in R3: the truth selects the calibration
  population, though it never sets a label.
- HG005's truth phases SNVs only, so the `indel` group has no pairs.

## Scorer and order

`benchmarks/r3_label_eligible.py`, generalized to take the trio from the
platform and to print the R3b verdicts, committed as a dated amendment with
a rerun on `cal_v5` that must reproduce R3's archived output
(`2026-10-16_r3_label_eligible.txt`). Then `PHASEBUDGET_PLATFORM=hg005_ont`
once. Deviations are dated amendments below.

## Amendments

### Amendment 1, 27 October 2026: the scorer rerun, and the score

The generalized scorer (8e82cab) reproduced R3's archived output on
`cal_v5` in all 164 result lines (the archive alone has a `/usr/bin/time`
trailer). The registered scoring ran once (`2026-10-27_r3b_score.txt`):
R3b-1 and R3b-2 are supported. Seen after scoring: because HG005's truth is
trio transmission and the rows hold only truth-phased pairs, 99.2 to 99.6
per cent of scored pairs are label-eligible, so R3b-2's trusted share
describes transmission-informative pairs only. The development notebook, 27 October.
