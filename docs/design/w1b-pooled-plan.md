# Pre-registered W1b: per-group thresholds calibrated on many genomes at once

Written on 15 October 2026. The commit that adds this file is the
registration; its push time is the third-party record. Creating this file is
what lets `benchmarks/genome_common.py` open the Platinum Pedigree HiFi reads
(`phifi_<sample>`).

## Why

Learn then Test at `delta = 0.05` certifies `beta` in a group only if the
group has at least `ln(0.05) / ln(1 - beta)` units, about 2,995 at
`beta = 1e-3`, and then only if none of them is wrong. One genome's
segmental-duplication group holds about 2,900 units with a pair, so it can
never certify `beta = 1e-3`: in every registered test it abstained there at
that beta (development notebook, 14 October, W3). Calibrating on the units of many genomes at
once removes that floor. W1's gate also showed that a truth built from one
genome's assemblies is too noisy to calibrate on (about 0.6 per cent of units
disagree with trio transmission); W1 amendment 1 replaces it with the
consensus of dipcall and 1000 Genomes transmission, which this test uses for
its HPRC genomes.

## Material seen before writing

Everything in the development notebook up to 15 October, including W1's gate and the
third test's scores. The same five Platinum Pedigree members were test
genomes of the second test on ONT reads (scored 13 October) and W5 part 1's
children; their HiFi reads have not been read. No HPRC test genome's
phasing has been compared with any truth, and no comparison of dipcall with
the 1000 Genomes panel has been run.

## Design

Everything in W1's frozen design holds (score, groups `indel`, `segdup`,
`lowmap`, `tandem`, `difficult`, `easy`, the 100 kb unit, the grid,
`delta = 0.05`, `beta` in {3e-4, 1e-3, 3e-3}, depths 30, 20 and 10 as drawn
there, both phasers with LongPhase `--pb`, the 21 confirmatory autosomes),
except as listed.

- **Calibration pool.** HG002 (HPRC HiFi, truth v5.0q, W1's calibration
  genome) and the first 20 HPRC genomes with a consensus truth in the order
  of `w1-candidates.tsv`'s `rank` column: HG00735, HG01243, HG03516, HG00438,
  HG02622, HG02145, HG02723, HG01358, HG02148, NA20129, HG01071, HG03453,
  HG01978, HG03579, HG02818, HG01109, HG00741, HG02572, HG01175, HG00621
  (9 AFR, 9 AMR, 2 EAS). Their truth is W1 amendment 1's consensus truth.
  For each depth and phaser, `learn_then_test_by_group` runs on the pooled
  rows, with a unit identified by (genome, chromosome, 100 kb index).
- **Test genomes, accurate truth.** Five Platinum Pedigree members, NA12879,
  NA12881, NA12882, NA12885 and NA12886, from their HiFi reads aligned to
  GRCh38 (`data/hifi/mapped/GRCh38/<sample>.GRCh38.haplotagged.bam` in the
  public `platinum-pedigree-data` bucket), scored against the Platinum
  Pedigree truthset inside its high-quality regions, with its heterozygous
  calls (phase removed) as phaser input, as the pedigree arm of the second
  test did. These are the five smallest of the seven members with a truth
  (54 to 116 GB; NA12877 and NA12878 hold 236 and 212 GB).
- **Test genomes, consensus truth.** The other 19 HPRC genomes with a
  consensus truth: HG02055, HG01106, HG00673, HG02630, HG02080, HG01928,
  HG01361, NA18906, HG01891, HG03540, HG01258, HG02717, HG02257, HG01952,
  HG02886, HG00733, NA19240, HG03486, HG03098 (11 AFR, 6 AMR, 2 EAS).
- **Shift stated now.** The Pedigree genomes differ from the pool in person,
  ancestry (European), sequencing project and alignment pipeline; they share
  the platform. The held-out HPRC genomes share pipeline and truth
  construction with the pool and differ in person.

## Hypotheses

A cell is (genome, depth, phaser, beta, group), decided as in W1: `n` test
units hold a trusted pair, `k` hold a wrong one, and the cell violates when
the exact binomial `P(X >= k | n, beta)` is below the level. A group that
certifies nothing abstains and cannot violate. Levels are fixed by the
counts here.

- **B1, validity on accurate truths (primary).** The five Pedigree genomes,
  all six groups: 5 x 3 x 2 x 3 x 6 = 540 cells, level 0.05/540. *Supported*
  if no cell violates. The scorer reports how many certify.
- **B2, validity on held-out HPRC genomes (primary).** The 19 genomes, all
  six groups, counting only the (group, beta) cells that W1 amendment 1's
  gate makes testable: 2,052 registered cells, level 0.05/2,052. *Supported*
  if no testable cell violates.
- **B3, the floor (primary).** At `beta = 1e-3`, the pooled calibration
  certifies the `segdup` group in all 4 cells of depths 30 and 20 and both
  phasers. Decided on the calibration alone.
- **B4, usefulness (primary).** Pedigree genomes, group `easy`,
  `beta = 1e-3`, depths 30 and 20, both phasers: 20 cells with the second
  test's S4 rule. *Supported* if all 20 succeed, *partly* if 18 do.
- **B5 (descriptive).** On the Pedigree cells, the trusted share per group
  with W1's thresholds (HG002 alone, v5.0q) against the pooled ones, and the
  share of pairs in groups that certify nothing, at `beta = 1e-3`.
- **B6 (descriptive).** A locus-merged calibration: the same pooled rows
  with a unit identified by (chromosome, 100 kb index) only, wrong when any
  genome is wrong there. It does not treat one locus in different genomes
  as independent and is conservative for a new genome exchangeable with
  the pool. Which (group, beta) it certifies, against B3.

## Risk stated now

Units of different genomes at one locus are not independent: a collapsed
duplication in GRCh38 can mislead the phasers in every genome. Pooled Learn
then Test treats them as independent, so its bound can be optimistic in
exactly the groups the pool is meant to help. B1 and B2 test whether this
matters on new genomes; B6 shows what a calibration without that assumption
certifies. A failure of B1 or B2 confined to `segdup` would point here.

## Scorer and order

- `benchmarks/w1b_score.py` decides B1 to B4 and reports B5 and B6, reading
  only stored fragments, phaser VCFs and truths, with the second test's cell
  functions (`stratified_score.cells`, `fewer_errors`) unchanged. It is
  committed, as a dated amendment, with a dry run that reads no genome but
  HG002: the pool is HG002's odd chromosomes split into two pseudo-genomes,
  and HG002's even chromosomes stand in for a test genome.
- W1 and W1b are scored in one session after their runs end (W1 amendment
  1): the input check, `w1_score.py`, then `w1b_score.py`, each run once.
- If a pool genome does not finish all 21 chromosomes, the pool is the
  genomes that did, and the count is reported. Test genomes are scored when
  all 21 chromosomes finished.
- Every threshold is compared at full precision from the run files. Every
  deviation is a dated amendment below, stating whether results had been
  seen.

## Software and machine

`benchmarks/confidence_run.py` with `PHASEBUDGET_PLATFORM=phifi_<sample>` for
the Pedigree genomes (LongPhase `--pb`), W1's runs for the HPRC genomes, on
the KCDH compute server with the third test's Python environment, WhatsHap 2.8
and LongPhase 2.0.2.

## Amendments
