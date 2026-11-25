# W4: does a complete reference remove confidently wrong relations?

Written 15 October 2026, before any W4 phasing ran. Development gate, not
registered. Direction W4 of the October directions memo (development repository) (Z4 there).

## Claim

Relations that the phaser gets wrong while every read agrees with it sit in
regions where GRCh38 is wrong (false duplications and collapsed
duplications, Aganezov et al. 2022), and aligning the same reads to
T2T-CHM13 turns those relations right.

## Data

HG002 HPRC HiFi, the primary records of HPRC's winnowmap GRCh38 BAM,
realigned with one aligner, minimap2 2.31 `map-hifi`, to GRCh38 no-alt and to
CHM13 v2.0 (`scripts/realign_hifi.sh`, `data/w4/`). Two platforms run the
usual per-chromosome pipeline (`benchmarks/confidence_run.py`) at 30-fold
only, on the 21 confirmatory chromosomes:

- `w4_grch38`: variants and truth GIAB v5.0q GRCh38 small variants inside
  its benchmark BED, as `hprc_HG002`;
- `w4_chm13`: the same from v5.0q's CHM13 v2.0 files.

Both phasers: `whatshap phase` and LongPhase `--pb`.

## Definitions

A **relation** is two heterozygous variants that the phaser puts in one
block and the truth phases. It is **wrong** when the phaser's relative phase
differs from the truth's. Its **support** is the reads of each consecutive
phased pair between the two variants (for a consecutive pair, that pair
alone). A relation is **confident** when every one of those pairs has at
least 10 spanning reads and no read disagrees with the phaser. This needs no
calibration, so it means the same thing on both references. At 10 reads the
frozen score's Wilson bound is 0.714.

**Region list**, fixed now: GIAB stratifications v3.6, GRCh38
`OtherDifficult`, the union of `GRCh38_false_duplications_correct_copy`,
`GRCh38_false_duplications_incorrect_copy` and
`GRCh38_collapsed_duplication_FP_regions` (md5 89e87651…, 24112e2c…,
483cd58b…). The v5.0q GRCh38 benchmark excludes the incorrect copies
entirely and keeps 0.61 of 8.0 Mb of collapsed regions, so on GRCh38 the
list is testable mainly through the correct copies (0.89 Mb in the BED). A relation is in
the list when either of its variants lies inside it.

**Matching across references.** A GRCh38 variant maps to CHM13 by the
T2T chain `grch38-chm13v2.chain` (md5 913bb9a3…, point liftover, both
strands). It matches a CHM13 input variant at the lifted position with the
same pair of alleles (reverse-complemented on a `-` chain). Only SNVs match;
indels are left out of the matched comparison because their positions
depend on normalization.

## Measurements

Per phaser:

1. Per reference: consecutive relations, wrong ones, confident wrong ones.
2. GRCh38: confident wrong relations by GIAB group (`segdup`, `lowmap`,
   `tandem`, `difficult`, `easy`), and inside the region list against the
   list's share of all relations.
3. Matched: GRCh38 consecutive SNV relations whose two variants both match.
   For each, the CHM13 relation between the two matched variants (not
   necessarily consecutive there): right, wrong, confident, or not in one
   block.

## Gate

Fixed before looking. The gate passes when, for both phasers:

- **removal:** on matched relations, CHM13 has at most half as many
  confident wrong relations as GRCh38, and at least half of GRCh38's
  matched confident wrong relations are right on CHM13;

and, for at least one phaser,

- **location:** at least 5 GRCh38 confident wrong relations fall in the
  region list, and the list's share of confident wrong relations is at
  least 5 times its share of all relations.

Removal is the claim that matters for phasing; location says why. Removal
without location is reported as a different mechanism.

## Code

`benchmarks/genome_common.py` and `benchmarks/confidence_run.py` define the
platforms `w4_grch38` and `w4_chm13`; `benchmarks/w4_gate.py` computes the
measurements and the verdict; `tests/test_w4_gate.py` checks the liftover,
SNV matching and relations on fixed examples.

## Limits

One genome, one read type. The truth on each reference comes from the same
Q100 assembly through dipcall, so a relation wrong on one reference and
right on the other is a reference effect only where both truths agree;
matched relations rest on both truths, and where they disagree that is
counted separately. GRCh38's worst regions are outside its benchmark, so
the location half sees the correct copies, where reads are lost to the
false copy, and little of the collapsed regions.
