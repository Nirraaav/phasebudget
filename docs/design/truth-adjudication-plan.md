# Pre-registered check: is the HG001 failure the truth's?

Written on 4 October 2026, before the Platinum Pedigree small-variant
truthset was downloaded. The commit that adds this file fixes that, and its
GitHub push time is the third-party record. This is track C of
the extension plan (development repository).

Material seen before writing:
- every result of `docs/design/confidence-run-plan.md`, scored on 3 October;
- the exploratory HG001 analysis inside GIAB regions
  (`2026-10-04_confidence_hg001_giab.txt`);
- of the Platinum Pedigree truthset
  (`s3://platinum-pedigree-data/variants/small_variant_truthset/GRCh38/CEPH1463.GRCh38.family-truthset.ov.vcf.gz`):
  its header, which lists the samples NA12877, NA12878, NA12885, NA12881,
  NA12879, NA12886 and NA12882 and only a `GT` field; and the genotype counts
  of about 3,000 records at the start of chromosome 1, read to check that
  genotypes are written with `|`. No record was compared with any phasing.

## Question

At the HG001 pairs that the registered run scored as wrong against the 2017
Platinum Genomes phase, does a second, independent-in-method truth agree with
the phaser or with the 2017 truth? The Platinum Pedigree truth (Kronenberg et
al. 2025) phases by inheritance across the CEPH-1463 family with HiFi, ONT
and Illumina data, and its high-quality regions reach difficult sequence.

## Data

- **Pairs.** The scorer's cached rows of the registered run
  (`data/confidence_rows/hg001_ont_{30,20,10}_{whatshap,longphase}_*.json`):
  each consecutive phased pair, whether it is wrong against Platinum 2017,
  and its score.
- **Relations.** For a pair, the 2017 relation (same or opposite haplotype)
  is read from `data/platinum/NA12878.vcf.gz`. The phaser's relation is the
  2017 relation, flipped when the pair is wrong. The Pedigree relation is read
  from the NA12878 column of the truthset.
- **Adjudicable pair.** Both SNVs are heterozygous, phased with `|`,
  biallelic with the same REF and ALT as the 2017 call, and inside
  `hq_regions_final.bed.gz`.
- **Outside GIAB.** At least one SNV of the pair is outside the GIAB HG001
  v4.2.1 benchmark BED (`scripts/data_hg001_giab.md5`).

## Hypotheses

- **C1 (primary).** Cell: `whatshap phase`, 30-fold. Among adjudicable wrong
  pairs outside GIAB, `x` of `m` have a Pedigree relation equal to the
  phaser's. With the exact (Clopper-Pearson) 95 per cent interval for
  `x / m`:
  - *truth errors dominate* if its lower end is above 0.5;
  - *phaser errors dominate* if its upper end is below 0.5;
  - *undecided* otherwise;
  - *not decidable* if fewer than 30 per cent of the wrong pairs outside
    GIAB are adjudicable.

  The same numbers are reported for the other five cells.
- **C2 (secondary).** The 18 registered HG001 cells are re-scored with the
  Pedigree truth: only adjudicable pairs, wrong meaning wrong against the
  Pedigree relation, and the registered calibration, grid, `delta` and level
  0.05/36. The number of violating cells is reported. If it is 0, the
  registered E1 failure is attributed to the 2017 truth outside GIAB; if it
  is 18, to the score; otherwise neither.
- **C3 (descriptive).** On the adjudicable pairs of the `whatshap phase`
  30-fold cell, the share where the two truths disagree, inside and outside
  the GIAB BED, regardless of the phaser.

## Limits stated now

- The two truths come from the same family and partly the same pedigree
  logic. The Pedigree truth also uses long reads, which can share a mapping
  artifact with the ONT reads being judged. Agreement of the Pedigree truth
  with a phaser is therefore evidence, not proof, of a 2017 truth error.
- Pairs that the Pedigree truth does not phase, or that lie outside its
  high-quality regions, are not adjudicated, and the adjudicable share is
  reported.

## Software

`benchmarks/adjudicate_truth.py`, written after this file is committed. It
reads only the cached rows, the two truths and the BEDs, and prints C1 to C3
with every parameter in its header. Deviations are dated amendments below.

## Amendments

None. After the registered output was read, an exploratory breakdown
(`adjudicate_truth.py --breakdown`) was added; it changes no registered
number.
