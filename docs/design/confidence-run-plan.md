# Pre-registered test of error-controlled phasing on new data

Written on 2 October 2026, before any HG001 or HG002 hac read was
downloaded. The commit that adds this file fixes that, and its GitHub push
time is the third-party record.

Material seen before writing:
- the follow-up plan (`docs/followup/error-controlled-phasing.md`) and its
  literature gate;
- the development gates G-A to G-D, the bin-dependence checks and the
  100 kb checks on HG002 ONT sup, HG002 HiFi and ONT chromosome 20
  (development notebook, 2 October);
- every PhaseBudget result.

No HG001 data, and no HG002 hac read, has been used.

## What is tested

A phaser's relation between two consecutive phased SNVs is *trusted* when
its score is at or above a threshold. The threshold is chosen on
calibration data so that the probability that a genome unit holds a wrong
trusted relation is at most `beta`. The test asks two things:
1. whether that bound holds on a different genome and a different
   basecaller, which is validity;
2. whether it is reached while keeping most relations and removing most
   errors, which is usefulness.

## Frozen design

- **Score.** `phasebudget.confidence.score` at commit `6cb6a74`: the Wilson
  lower bound (z = 2) on the share of reads calling both SNVs that agree
  with the called relation, with ties broken by the number of such reads.
  Reads are extracted as `whatshap phase` extracts them: realignment, MAPQ
  20, no supplementary alignments, read groups ignored.
- **Unit.** A 100 kb bin, `(chromosome, floor(pos1 of the pair's first
  SNV / 100,000))`.
- **Loss.** A unit is wrong at a threshold when at least one of its trusted
  pairs is wrong against the truth.
- **Grid.** `phasebudget.confidence.coverage_grid` on the calibration
  scores, from half the pairs trusted to all of them in steps of 0.25 per
  cent, ordered strict to loose.
- **Procedure.** `phasebudget.confidence.learn_then_test`, with exact
  binomial p-values, at `delta = 0.05`.
- **Targets.** `beta` in {3e-4, 1e-3, 3e-3}.
- **Calibration.** HG002 ONT R10.4.1 sup 60x, the 21 confirmatory autosomes
  of the ONT run, truth GIAB v5.0q (genome-run amendment 2). It uses the
  run's stored `realign` fragments and its `whatshap phase` and LongPhase
  outputs, at the same depth and with the same phaser as the test cell. If
  no calibration threshold is certified, the cell is reported as "none
  certified" and counts as neither a violation nor a success.

## Test data

| dataset | reads | variants | truth |
|---|---|---|---|
| HG001 ONT | `ont-open-data` `giab_2023.05/analysis/variant_calling/hg001_sup_60x/hg001.haplotagged.bam` | biallelic het SNVs of Platinum Genomes 2017-1.0 NA12878 (GRCh38) inside its `ConfidentRegions.bed.gz` | the same calls' phased genotypes, pedigree-phased; files checked by `scripts/data_platinum.md5` |
| HG002 ONT hac | `.../hg002_hac_60x/hg002.haplotagged.bam` | the ONT run's `hetsnv.vcf.gz` (GIAB v4.2.1 het SNVs in its BED) | GIAB v5.0q |

- **Chromosomes.** The 21 confirmatory autosomes. Chromosome 20 enters no
  test.
- **Depths.** 30, 20 and 10.
  - The 30-fold slice is drawn from the remote BAM with `samtools view
    --subsample`, seed 0. Its fraction is 30 over the depth estimated from
    non-secondary aligned bases in the middle 5 Mb of the chromosome.
  - The 20- and 10-fold slices are drawn from the 30-fold slice with
    fractions 2/3 and 1/3, seed 0.
  - Realized depths are measured on each slice and reported.
- **Phasers.** `whatshap phase` 2.8 (`--reference`, `--ignore-read-groups`,
  default `--internal-downsampling 15`) and LongPhase 2.0.2 (`--ont`,
  `-t 1`), the versions of the calibration runs.

## Hypotheses

There are 36 cells: two datasets, three depths, two phasers and three
values of `beta`.

- **E1, validity (primary).** In each cell, `n` test units hold a trusted
  pair and `k` hold a wrong one. A cell violates the bound when the exact
  binomial `P(X >= k | n, beta)` is below 0.05/36.
  - *Supported* if no cell violates.
  - The realized rate `k / n` and the number of cells with `k / n > beta`
    are reported for every cell.
- **E2, usefulness (primary).** Eight cells: `beta = 1e-3`, depths 20 and
  10, both phasers, both datasets. A cell succeeds when both hold:
  - at least 90 per cent of test pairs are trusted;
  - wrong trusted pairs are at most half the wrong pairs over all pairs.

  *Supported* if all 8 cells succeed, and *partly supported* if 6 or 7 do.
- **E3 (descriptive).** Risk-coverage curves for each phaser on each test
  genome: the pair error rate when 95, 98, 99 and 100 per cent of pairs are
  trusted, ranked by the score. The two phasers are compared at matched
  coverage.
- **E4 (descriptive).** The base unit rates, with every pair trusted, and
  the share of wrong pairs in GIAB v3.6 difficult-region strata among the
  pairs left untrusted.

## Software and safety

- `benchmarks/confidence_run.py` runs the data steps, one chromosome at a
  time:
  - slices and subsamples;
  - extracts fragments at each depth;
  - runs both phasers;
  - deletes the BAMs.

  A download starts only with at least 10 GB more free disk than the
  estimated 30-fold slice.
- `benchmarks/confidence_score.py` reads only stored fragments, phaser
  VCFs and truths. It calibrates on HG002 ONT sup and decides E1 and E2.
- Every deviation from this plan is a dated amendment below, stating
  whether results were seen.

## Amendments

1. **2 October 2026, about 21:30 IST, before any test output was scored or
   compared with a truth.**
   - *Defect.* The first run's subsampling was nested: the 20- and 10-fold
     slices were drawn from the 30-fold slice with the same seed (0) as the
     30-fold draw. Samtools keeps a read when a hash of its name and the
     seed falls below the fraction. So the 2/3 draw kept every read, and the
     1/3 draw kept 68 per cent of them. On HG001 chromosome 21, the only
     chromosome finished, the slices measured 37.4, 37.4 and 25.5 (mean
     aligned depth over the chromosome).
   - *Second defect.* The 30-fold slice was 37.4, not 30, because a single
     middle 5 Mb window misestimates the chromosome's depth.
   - *Seen.* Only `depth.json`. No phasing was scored, and no truth was
     read.
   - *Change:*
     - depth is the median, over 1 Mb windows with any aligned base, of
       non-secondary aligned depth;
     - the remote estimate uses ten 500 kb windows spaced evenly over the
       chromosome;
     - the 30-fold fraction is 30 over that estimate;
     - the 20- and 10-fold slices are drawn from the 30-fold slice with
       fractions 20 and 10 over its measured depth, using seeds 1 and 2,
       which are independent of the first draw;
     - realized depths of every slice are reported.
   - *Discarded.* The chromosome 21 and 22 outputs of HG001 are deleted and
     redone.
   - *Addendum, same evening, still before any scoring.* A test of the
     changed steps on chromosome 21 found two things:
     - the remote draw with seed 0 keeps more reads than its fraction. At
       0.4986 it kept 425 of 666 read names in one window, and the slice
       measured 39.9. Samtools' seeded hash is not uniform over these ONT
       read names at seed 0;
     - local draws with seeds 1 and 2 met their targets (19.9 and 10.1).

     The 30-fold slice is therefore drawn in two stages:
     - from the remote BAM at fraction `min(1, 36 / estimate)`, seed 0;
     - then, when its measured depth exceeds 30, locally to 30 over that
       depth with seed 3.

     The 20- and 10-fold slices follow from it, as above. Only `depth.json`
     was read.

2. **3 October 2026, 06:57 IST (the commit adding this amendment), before any test output was scored or
   compared with a truth. HG001 had 20 chromosomes finished, and HG002 hac
   had not started.**
   - *What prompted it.* A development comparison
     (`benchmarks/compare_hapcut2_quality.py`,
     `2026-10-03_hapcut2_quality_ont.txt`) ranked HapCUT2's own phasings of
     HG002 ONT sup two ways: by its per-variant quality (blocks column 11,
     the VCF's PQ) and by the frozen score. At 99 per cent of pairs trusted,
     the wrong pairs left were:

     | depth | frozen score | HapCUT2 quality |
     |---|---|---|
     | 30 | 1 | 134 |
     | 20 | 2 | 128 |
     | 10 | 24 | 169 |

     The AUCs were 1.000, 1.000 and 0.997 against 0.839, 0.864 and 0.925.
     Those data were seen.
   - *Addition, E5 (secondary).* On HG002 ONT hac only, at depths 30, 20
     and 10, the runner also phases with HapCUT2. It uses the calibration
     runs' commands: `extractHAIRS --ont 1 --ref` and then `HAPCUT2` with
     defaults, and pruned variants count as unphased. HapCUT2's pairs are
     ranked by the frozen score and, separately, by the smaller of the two
     variants' HapCUT2 qualities.
     - *Supported* if the frozen score leaves fewer wrong pairs than
       HapCUT2's quality at 99 per cent trusted, at all three depths.
     - The 99.5 per cent results are reported beside them.
     - Learn then Test on HapCUT2's phasings, calibrated on HG002 ONT sup
       HapCUT2, is reported descriptively.
   - *Not added for HG001.* Its reads were deleted after phasing, and
     fetching them again only for this baseline is outside the plan.
   - *Unchanged.* E1 to E4.

