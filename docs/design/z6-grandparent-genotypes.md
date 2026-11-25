# Z6 data: NA12878's parents genotyped at her heterozygous sites

Written on 7 October 2026, before any NA12891 or NA12892 read was read.
Direction Z6 (the October directions memo (development repository)) calibrates the
Pedigree arm on transmission labels instead of the Pedigree truth. Y1 step 3
found that labels from called genotypes, filtered for quality, calibrate as
safely as the truth on the HG002 trio (development notebook, 7 October). This run
produces the parents' genotypes the labels need. NA12891 and NA12892 are not
test genomes of X6.

## Sites

NA12878's biallelic heterozygous SNVs in the Platinum Pedigree truthset,
inside its high-quality regions, on the 21 confirmatory autosomes: the SNVs
of the `ped_NA12878` input (`extract_fragments.het_snvs`, phase removed).
Only the child's genotype is used, never its phase. That the child's sites
come from the truthset matches how X6 defines its population; the parents'
genotypes come from reads alone.

## Reads and counts

`NA12891.minimap2.bam` and `NA12892.minimap2.bam` from
`platinum-pedigree-data/data/ont/mapped/GRCh38/`, streamed per 10 Mb region
with `samtools mpileup -B -q 20 -Q 0 --ff UNMAP,SECONDARY,QCFAIL,DUP
--positions`, restricted to the sites. Per site: reads showing the reference
base and reads showing the alternative base. Nothing else is kept; no read is
stored.

## Genotype rule, fixed now

With `d` the reference plus alternative reads and `f` the alternative share:
- no call if `d < 8`;
- homozygous reference if `f <= 0.10`;
- homozygous alternative if `f >= 0.90`;
- heterozygous if `0.25 <= f <= 0.75`;
- no call otherwise.

A parent's no call makes the site unusable for that parent, as a low-quality
call did in Y1 step 3. Which parent is the father does not change any
relation between two labelled sites.

## Use

Labels by `phasebudget.truth.phase_by_transmission`, then calibration as in
Y1 step 2, with the Pedigree truth as the reference loss. How the labels
enter X6 (as a second calibration of the Pedigree arm, scored on the six test
relatives) is fixed by X6's registration, not here.

## Amendment 1 (7 October, before the run, after a 1 Mb test)

The local `samtools` 1.21 has no libcurl and cannot read https, and pysam's
`mpileup` wrapper exits silently on a remote file. Counts are therefore
taken with `pysam.AlignmentFile.count_coverage` in 1 Mb windows, with the same
filters (mapping quality at least 20; no unmapped, secondary, QC-failed or
duplicate reads; any base quality). The genotype rule is unchanged. A test
on NA12891 `chr22:20,000,001-21,000,000` (1,129 sites, 33 s) read reads
before this amendment; it showed only the genotype mix (323 `0/0`, 621
`0/1`, 167 `1/1`, 18 no calls) and was deleted, to be redone by the run.

## Development comparison on NA12878 (stated 8 October, before running)

Before X6 is registered, the labels are checked on NA12878's own data, as Y1
step 3 checked the HG002 trio's. No test genome is read.

- **Labels.** `phase_by_transmission` with NA12878's heterozygous SNVs as the
  child (genotype only) and NA12891 and NA12892 genotyped by the rule above;
  a parent's no call makes the site unusable for that parent.
- **Pairs.** `ped_NA12878` consecutive SNV pairs (pairs with an indel are
  left out: transmission labels here exist for SNVs only), at 30, 20 and 10x,
  both phasers, whose two SNVs both carry a label.
- **Losses.** Against the Platinum Pedigree truth (`truth`) and against the
  labels (`labels-called`); the phaser's relation is recovered from the
  row's loss against the truth.
- **Calibration and test.** Calibrate on the odd chromosomes, pooled and per
  stratum (`segdup`, `lowmap`, `tandem`, `difficult`, `easy`), `beta` in
  {3e-4, 1e-3, 3e-3}; score on all SNV pairs of the even chromosomes against
  the truth.
- **Rule.** The labels pass if (a) in no cell do their thresholds violate
  (p < 0.05) in more groups than the truth's, and (b) at `beta = 1e-3` and 30x
  their trusted share on the test is within 0.01 of the truth's for both
  phasers in both modes.
- **Also reported.** Labelled sites, label disagreements with the truth
  among consecutive labelled sites, and per cell the label errors that hide
  a phaser error and that invent one.

## Repair: a stricter genotype rule chosen on the odd chromosomes (stated 8 October, before running)

The development comparison failed (development notebook, 8 October): the pileup genotypes
hid 4 to 6 phaser errors per cell. Three rules are considered, the stored
counts reused:

| rule | min depth | 0/0 if alt share <= | 1/1 if alt share >= | 0/1 if alt share in |
|---|---:|---:|---:|---|
| R0 (fixed above) | 8 | 0.10 | 0.90 | [0.25, 0.75] |
| R1 | 15 | 0.05 | 0.95 | [0.30, 0.70] |
| R2 | 20 | 0.03 | 0.97 | [0.35, 0.65] |

**Choice, odd chromosomes only.** Among rules keeping at least 80 per cent of
R0's labelled SNVs on the odd chromosomes, the one with the fewest
disagreements with the Pedigree truth among consecutive labelled sites there;
ties go to the earlier rule.

**Evaluation.** The development comparison above, unchanged (calibrate on
the odd chromosomes, score on the even ones, the same (a) and (b)), with the
chosen rule. Labels on the even chromosomes are not used. The odd
chromosomes serve both the choice and the calibration, as training; the
even chromosomes are the test.
