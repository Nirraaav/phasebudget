# Q3: how much phasing error lies outside the benchmark

Written 17 October 2026, after the phasing run started and before any relation
was scored against a truth. Development gate, not registered. Direction Q3 of
the October directions memo (development repository), round 6.

## Question

Every registered test scored relations inside a truth's benchmark BED. For
HG002, the GIAB v5.0q small-variant BED covers 2,521 Mb of the 21
confirmatory autosomes; the Q100 dipcall BED
(`GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.bed`) covers 2,607 Mb. Q3 asks how
many phasing errors fall in the part a v5.0q-based benchmark cannot see, how
the frozen-score certificates behave there, and what causes the errors.

Facts checked before scoring (counts of input variants, no phase read): the
Q100 dipcall VCF has 2,585,316 biallelic heterozygous SNVs and indels inside
its BED and inside v5.0q's BED on those chromosomes, and 235,997 inside its
BED but outside v5.0q's. v5.0q's VCF holds 1,212,575 records outside its own
BED, so the two call sets can be compared by position there.

## Data

HG002 HPRC HiFi, realigned with minimap2 2.31 `map-hifi` to GRCh38
(`data/w4/HG002_grch38_mm2.bam`, W4). Platform `q3_q100` in
`benchmarks/genome_common.py` runs the usual per-chromosome pipeline
(`benchmarks/confidence_run.py`) at 30-fold only: the same estimate, seed-0
draw and seed-3 trim as `w4_grch38`, so the 30-fold reads are W4's. The
variants are Q100's biallelic heterozygous SNVs and indels of at most 50
bases inside the Q100 dip BED, with phase removed (`het_snvs`, `indels=True`,
`unphase=True`, sample `syndip`). Phasers as in W4: `whatshap phase` 2.8 at
defaults with `--reference` and `--ignore-read-groups`, and LongPhase 2.0.2
with `--pb --indels -t 1`. The 21 confirmatory chromosomes decide; chromosome
20 runs afterwards and is reported as descriptive.

## Definitions

- **Relation**: two consecutive variants that the phaser phases in one block,
  both phased `0|1` or `1|0` by Q100 at the same alleles inside the Q100 dip
  BED (`genome_common.dipcall_truth` on the Q100 files). It is **wrong** when
  the phaser's relative phase differs from Q100's.
- **Inside**: both variants inside v5.0q's benchmark BED. **Outside**: any
  other relation, reported also as one variant inside and none inside.
- **Group**: `indel` when either variant is an indel, otherwise the stratum
  of `stratified_gate.strata` (`segdup`, `lowmap`, `tandem`, `difficult`,
  `easy`, first that holds), as in the stratified run.
- **Error per pair**: wrong relations over relations. **Unit**: 100 kb,
  `(chromosome, first position // 100,000)`, counted separately for inside
  and outside; its wrong share is units holding a wrong relation over units
  holding a relation.
- **Share outside**: wrong relations outside over all wrong relations.

## Gate

Q3 **finds a gap** when, for both phasers at 30-fold on the 21 chromosomes:

- the error per pair outside is at least 3 times the error per pair inside,
  and
- at least 20% of all wrong relations are outside.

The suggestion for this gate used the per-unit wrong share. The gate uses the
per-pair error instead because a unit's wrong share depends on how many
relations the unit holds, and outside-BED relations are spread thinly over
units that also hold inside relations: the same per-pair error gives a lower
unit share outside. In W4, whatshap had 12,601 wrong relations on the roughly
25,000 units of 100 kb these chromosomes hold, so a large part of the inside
units may hold a wrong relation, and a unit share cannot exceed 1: the unit
ratio can stay below 3 whatever the outside error is.
The per-unit ratio is printed beside the gate and does not decide it.

Prediction, stated now: the gate passes for whatshap; for LongPhase the share
outside passes and the ratio is uncertain.

## Certificates outside the benchmark

Thresholds come from `learn_then_test_by_group` on the cached `w4_grch38`
30-fold rows (the same reads, v5.0q's variants, scored against v5.0q inside
its BED, 21 chromosomes, grouped as above), betas 3e-4, 1e-3 and 3e-3, delta
0.05, grid from 0.5 in steps of 0.0025, 100 kb units per group. They are
applied to the Q100-scored relations inside and outside separately, as
`stratified_score.cells`. Inside cells are close to in-sample, since the reads
and most relations are the calibration's. An outside cell **violates** when
`P(X >= k | n, beta)` is below 0.05 over the number of outside cells with a
threshold and a unit; cells below 0.05 are listed too. Prediction: at least
one outside cell violates at beta 1e-3 for whatshap. This is reported and is
not part of the gate.

## Cause split

For every wrong relation, inside and outside:

- **spanning reads** are the reads of the `realign` fragments (WhatsHap's
  allele detection, all reads, before any selection) that call both variants;
  `agree` back the phaser's relation, `disagree` back the truth's. The
  majority is `wrong` (agree > disagree), `tie`, or `right`.
- **window optimum**: the fragment columns from 10 variants before the first
  variant to 10 after the second, without the columns strictly between them
  (so the pair is adjacent), with each read cut to its calls there; reads
  chosen by `phasebudget.selection.select_whatshap` (the library's port of
  WhatsHap 2.8's read selection) at coverage 15, the `whatshap phase` default
  `--internal-downsampling 15`; then the exact constrained optima
  `CoverageDP.boundary_costs` with the relation forced each way. The optimum
  is `right` when the true relation costs strictly less, `wrong` when it costs
  strictly more, `tie` otherwise. `tests/test_q3_outside_benchmark.py` checks
  the window and the constrained optima against exhaustive enumeration.

Causes, a partition:

1. **none** (i): no read calls both variants;
2. **data** (ii): the majority is `wrong` or `tie` and the window optimum is
   not `right`;
3. **search** (iii): the window optimum is `right` (the reads, under the
   phaser's own model and cap, prefer the truth; the phaser did not find it);
4. **other**: the majority is `right` but the window optimum is not, which
   happens when other reads in the window or the selection outweigh the
   spanning reads.

The full majority by optimum table is printed. LongPhase does not export its
allele calls, so its wrong relations are examined on the same WhatsHap
fragments; for LongPhase, `search` means the WhatsHap-model optimum gets the
relation right, not that LongPhase's own search failed on its own evidence.
Prediction: outside the BED, `data` and `none` together are at least half of
the wrong relations for both phasers. Reported regardless of the gate.

## Truth check

Q100 is an assembly truth and v5.0q was built from it, so their agreement is
not independent evidence. For each outside wrong relation, both variants are
looked up in v5.0q's VCF (all records, BED ignored) at the same position and
alleles: **genotype differs** when either is missing or not heterozygous in
v5.0q, else the relation is compared with v5.0q's phase (same, other, or
v5.0q unphased). A wrong relation where the genotype differs sits where the
two GIAB call sets disagree, and may be a truth error rather than a phasing
error.

## Code

`benchmarks/q3_outside_benchmark.py`; the platform in
`benchmarks/genome_common.py` and `benchmarks/confidence_run.py`; tests in
`tests/test_q3_outside_benchmark.py`. Output:
`benchmarks/results/runs/2026-10-17_q3_outside_benchmark.txt`.

## Limits

One genome, one read type, one depth. Q100 is the truth everywhere, and v5.0q
derives from the same assembly, so an error shared by both is invisible. The
outside regions are where v5.0q's authors found the assembly, the reference
or the alignments less reliable, so some outside wrong relations may be truth
errors; the truth check counts the ones where the two call sets differ, not
the ones where both are wrong. The cause split uses a window of 21 columns
and WhatsHap's model for both phasers.
