# Pre-registered test of per-stratum error control on seven new genomes

Written on 8 October 2026. The commit that adds this file is the
registration; its push time is the third-party record. Creating this file is
what lets `benchmarks/genome_common.py` open the test genomes.

## Material seen before writing

Everything in the development notebook up to 8 October, in particular:
- the registered confidence test (E1 to E5) and track C;
- gate A within HG001 and across genomes (HG002 to HG001): calibrated on the
  registered HG002 rows (v4.2.1 SNVs), the easy stratum held in 18 of 18
  cells and the difficult strata failed in 10; calibrated on `cal_r10`
  (v5.0q SNVs, which reach difficult regions), no stratum violated in any
  cell, and per-group thresholds trusted 0.83 to 0.95 of HG001's pairs at
  `beta = 1e-3` where one pooled threshold trusted none (whatshap) or 0.62
  to 0.97 (LongPhase) (development notebook, 8 October);
- the calibration genomes of both arms (`cal_r10`, `ped_NA12878`): their
  per-group rates, a within-genome odd/even check, and a dry run of this
  plan's scorer on `ped_NA12878`;
- Y1 (calibration from a trio's Clair3 calls passed on development data), Y2
  (certified relations in clinical genes; the duplicated-gene association),
  Z1, Z2, Z6 and Z7 (stopped), and whatshap's indel errors on NA12878.

No read of a test genome has been used. Two test-genome BAM **indexes** were
fetched without reads: `NA12879.minimap2.bam.bai` when the Pedigree BAM
headers were checked on 4 October, and its first 20 MB in a network speed
test on 7 October (development notebook, 7 October). HG005's BAM headers and the
release's depth summaries were read on 8 October.

## What is tested

The procedure is Learn then Test calibrated separately within each group, the
conditional inductive conformal construction of Vovk 2012 (arXiv
1209.2673) applied to Learn then Test's risk. It is not claimed as new. The
test asks where its guarantee carries over to genomes not used in
calibration: to relatives sequenced by the same project (Pedigree arm), and
to an unrelated genome sequenced with a newer basecaller (R10 arm).

## Frozen design

- **Score.** `phasebudget.confidence.score`, unchanged since `6cb6a74`.
- **Group of a consecutive pair.** `indel` if either variant is an indel;
  else the first of `segdup`, `lowmap`, `tandem` that holds for either SNV,
  else `easy` if both SNVs are inside `notinalldifficultregions`, else
  `difficult` (GIAB v3.6 GRCh38, md5 sums in `benchmarks/stratify.py`).
- **Unit.** `(chromosome, 100 kb bin of the pair's first variant, group)`.
- **Procedure.** `phasebudget.confidence.learn_then_test_by_group`, grid from
  0.5 in steps of 0.0025, `delta = 0.05`, `beta` in {3e-4, 1e-3, 3e-3}. A
  group with no certified threshold trusts nothing.
- **Population.** Consecutive pairs of phased variants that the arm's truth
  phases. The input variants are the truth's own heterozygous calls with
  phase removed; genotype errors are outside the population.
- **Calibration.** All 21 confirmatory autosomes of the arm's calibration
  genome, at the same depth and with the same phaser as the test cell.

## Arms and data

| arm | calibration | test | truth |
|---|---|---|---|
| Pedigree | `ped_NA12878`: Platinum Pedigree ONT, NA12878 | NA12877, NA12879, NA12881, NA12882, NA12885, NA12886 (`platinum-pedigree-data/data/ont/mapped/GRCh38/<sample>.minimap2.bam`) | Platinum Pedigree small-variant truthset, phased by inheritance, inside `hq_regions_final.bed.gz` |
| R10 | `cal_r10`: HG002 ONT sup, `giab_2023.05` `hg002_sup_60x`; v5.0q heterozygous SNVs and indels in its BED | HG005 ONT, `giab_2025.01/basecalling/sup/HG005/{PAW87816,PAW88001}/calls.sorted.bam`, R10.4.1, `dna_r10.4.1_e8.2_400bps_sup@v5.0.0`, 26.7x and 28.3x; both flowcells are drawn from with the same fraction and merged | GIAB v4.2.1 HG005 heterozygous SNVs and indels in its BED as input; transmission phase from HG006 and HG007 (SNVs only, `data/giab_hg005/truth/`) |

- **Shift stated now.** The R10 test differs from its calibration in person
  and in basecaller version. HG005's truth phases SNVs only, so the R10 arm
  has no scored indel pairs.
- **Depths** 30, 20 and 10, drawn as in amendment 1 of the confidence plan
  (a seed-0 draw aiming at 36x, trimmed to 30, then 20 and 10 from it).
- **Phasers** `whatshap phase` 2.8 (`--reference`, `--ignore-read-groups`)
  and LongPhase 2.0.2 (`--ont --indels -t 1`).
- **Chromosomes.** The 21 confirmatory autosomes; chromosome 20 enters no
  test.

## Hypotheses

A cell is (test genome, depth, phaser, beta, group). In a cell, `n` test
units hold a trusted pair and `k` hold a wrong one; it violates when the
exact binomial `P(X >= k | n, beta)` is below the hypothesis's level. Cells
of groups that certify nothing cannot violate and are reported as abstaining.

- **S1, the easy stratum carries over (primary).** Group `easy`, all seven
  test genomes: 7 x 3 x 2 x 3 = 126 cells, level 0.05/126. *Supported* if no
  cell violates.
- **S2, difficult groups within one project (primary).** Pedigree arm,
  groups `segdup`, `lowmap`, `tandem`, `difficult`, `indel`: 540 cells, level
  0.05/540. *Supported* if no cell violates.
- **S3, difficult groups across a basecaller change (primary).** R10 arm,
  the same groups: 90 cells, level 0.05/90. *Supported* if no cell violates.
  (A first draft predicted failure here, from the registered calibration
  rows; the `cal_r10` development result reversed that before
  registration.)
- **S4, usefulness (primary).** Group `easy`, `beta = 1e-3`, depths 30 and 20,
  both phasers, all test genomes: 28 cells. A cell succeeds when at least 99
  per cent of the group's test pairs are trusted and either wrong trusted
  pairs are at most half of its wrong pairs (when it has at least 20) or
  under 1 per 10,000 trusted pairs are wrong (when it has fewer).
  *Supported* if all 28 succeed, *partly* if at least 25 do.
- **S5 (descriptive).** Per cell, the share of pairs in groups that certify
  nothing.
- **S6, duplicated clinical genes withhold more (secondary).** 30x, `beta =
  1e-3`, per-group thresholds, each test genome and phaser (14 cells, level
  0.05/14): among genes with at least 10 ClinVar pathogenic or likely
  pathogenic alleles (`gene_specific_summary.txt` of 4 October 2026, md5
  `9761c2d9b82a96d8601b6129d607b312`) and at least 100 relations between
  truth-phased variants (`docs/design/chain-relations.md`), the trusted
  shares of genes with an unprocessed pseudogene named after them in
  GENCODE v50 basic are lower than the others' (one-sided Mann-Whitney U).
  *Supported* if all 14 are significant, *partly* if at least 10 are.
- **S7, whatshap's indel errors (secondary).** Each Pedigree test genome at
  30x, indel pairs both phasers phase: pairs only whatshap gets wrong
  outnumber pairs only LongPhase gets wrong (one-sided binomial, level
  0.05/6). *Supported* if all 6 are significant.

- **S8, per-group against pooled (primary).** Each test genome, both
  phasers, 30x and 20x, `beta = 1e-3` (28 cells): the share of all test
  pairs trusted by the per-group thresholds exceeds the share trusted by
  one pooled threshold calibrated on the same rows, all groups including
  `indel` pooled. *Supported* if all 28 cells, *partly* if at least 21. In
  the dry run on the calibration genomes the pooled threshold, pulled by
  indel pairs, certified nothing in 6 of 8 cells, so S8 mainly tests whether
  grouping lets indel-heavy data be certified at all.

## Software

- `benchmarks/confidence_run.py` and `benchmarks/genome_common.py` at the
  registration commit: platforms `hg005_ont` and `ped_<sample>`, three
  chromosomes at a time per genome, disk reserved per chromosome
  (`data/reservations/`), two phasing slots, reads deleted after phasing.
- `benchmarks/stratified_score.py` at the registration commit decides S1 to
  S4, S6 and S7 and reports S5, reading only stored fragments, phaser VCFs
  and truths. Its `--dry-run` was run on the calibration genomes before this
  commit, on both arms (`tmp/stratified_score_dryrun_both.txt`, summarized
  in the development notebook of 8 October); S6 had low power there (1 of 4 significant on
  half of each genome), which is a stated risk. Every threshold is scored at
  full precision from the run file, never from printed text.
- Every deviation is a dated amendment below, stating whether results had
  been seen.

## Amendments

### Amendment 1 (8 October, before any test-genome read)

Phasing slots: up to four instead of two, the third and fourth taken only
while at least 6 GB of memory is available. Scheduling only; no step,
parameter or output changes.

### Amendment 2 (8 October, before any test-genome read)

The phaser input keeps the GT field alone. HG005's GIAB v4.2.1 records also
carry IGT, IPS and PS (the original callsets' genotypes and phase), which the
first build of HG005's chromosome 21 and 22 inputs kept beside an unphased
GT. Those two inputs were deleted before any read was drawn and are rebuilt.
The Pedigree inputs held GT alone and are unchanged; `cal_r10`'s held GT and
AD, which carries no phase.
