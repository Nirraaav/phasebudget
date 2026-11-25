# Q1: siblings without parents as calibration labels

Written 17 October 2026, before the gate computation. Development gate, not
registered. Direction Q1 of the October directions memo (development repository); it
follows R2 (a new batch needs about `ln(1/delta)/beta` labelled units of its
own) and R3 (Learn then Test on label-eligible pairs is valid for that
population).

## Statement

Children A and B of the same parents share 0, 1 or 2 haplotypes identical
by descent at each position (IBD0, IBD1, IBD2). In an IBD1 segment, a site
where A is heterozygous and B homozygous for allele `x` puts `x` on A's
shared haplotype. Two such sites `a`, `b` in one segment give A's relation
between them: the same first allele exactly when `x_a == x_b`, whichever of
A's haplotypes is the shared one. IBD0 and IBD2 segments give nothing: in
IBD0 B's alleles say nothing about A's phase, and in IBD2 A and B have the
same genotype.

Segments come from the two unphased genotypes alone. A three-state HMM runs
over the sites where A or B carries an alt allele. Its observation is the
genotype class: opposite homozygotes (IBS0), one heterozygous and one
homozygous (IBS1), both heterozygous, or both homozygous alt. Transitions
model the two parental IBD indicators as independent two-state chains, each
flipping at `lambda = 2e-8` per bp (a crossover in either of two meioses at
about 1 cM/Mb), so `P(stay) = (1 + exp(-2 lambda d)) / 2` per indicator over
`d` bp; the initial distribution is `(1/4, 1/2, 1/4)`. Emissions are fitted
per sibling pair by EM over all 21 chromosomes jointly (4 iterations,
transitions fixed). No allele frequencies are needed, and the genotype error
rate is part of the fitted emissions: IBS0 under IBD1 and IBD2, and IBS1
under IBD2, have no other source.

A **run** is a maximal stretch of consecutive HMM sites with posterior
`P(IBD1) >= 0.99`. A site is labelled by sibling B when A is heterozygous
and B homozygous there, it lies in a run, and it is at least 100 kb from
both ends of the run. A consecutive pair of A's phased sites gets B's
relation when B labels both sites in the same run. With several siblings a
pair is **eligible** when at least one sibling gives a relation and every
sibling that gives one agrees; a pair with disagreeing siblings is not
eligible. Eligibility uses the genotypes of A and its siblings only, so a
user knows it before any truth. As in R3, Learn then Test calibrated on the
eligible pairs of some units controls the wrong-unit share among the
eligible pairs of exchangeable units, and the procedure abstains on
ineligible pairs.

What the labels cannot do: they reach a pair only where some sibling is
IBD1 with A and homozygous at both sites, so a pair with both sites in a
region where every sibling is IBD0 or IBD2 with A, or where every sibling
is heterozygous, is never labelled. Whether such pairs are harder for
phasers is not assumed; it is reported per group as the eligible share.

## Threshold and margin, fixed on simulation

`tests/test_q1_relative_labels.py` and a one-off simulation (17 October,
server, not kept) fixed `0.99` and `100 kb` before any Pedigree data was
labelled. Two unrelated parents with independent sites (no linkage
disequilibrium), alt frequency `10^U(-2.3, 0)` capped at 0.95, one site per
600 bp, ten chromosomes of 150 Mb, Poisson crossovers at 1e-8 per bp per
meiosis, five sibling pairs. With exact genotypes one sibling labels 0.12 of
a child's consecutive heterozygous pairs. No label was wrong at posterior
thresholds 0.99 and 0.999 with any margin from 0 to 100 kb (about 282,000
labelled pairs). At 0.9 with no margin, 4 were wrong. A 100 kb margin costs
1% of the labels. With random genotype errors at rate 1e-3 and 1e-2 the
label error was 2.0e-3 and 1.9e-2 whatever the threshold or margin: a
wrong homozygous call inside a run is invisible to the HMM, and each pair
has two sites. On real data with linkage disequilibrium and gaps between
the high-quality regions the boundaries are less sharp, which is why the
margin is kept at 100 kb.

## Measurement

Children: the five Platinum Pedigree ONT runs `ped_NA12879`, `ped_NA12881`,
`ped_NA12882`, `ped_NA12885`, `ped_NA12886` (development data, scored in X6).
They are the five children in the truthset, so each child has four siblings
there. Sibling and child genotypes are the Pedigree truthset's
(`CEPH1463.GRCh38.family-truthset.ov.vcf.gz`) with phase removed: biallelic
records inside `hq_regions_final.bed`, one record per position (positions
with two records are dropped). Parents and grandparents are not used.

Versions of the genotypes:

- `exact`: the truthset genotypes (optimistic labels; gated);
- `e1e-3`, `e1e-2`: each genotype of each sample replaced by one of the
  other two with that probability, independently per site, seed 17, the
  same perturbed genotypes for a sample in every pair (descriptive).

Sibling sets: `all` four siblings (gated); `two` and `one`, the first two and
the first of the other children in the order NA12879, NA12881, NA12882,
NA12885, NA12886 (descriptive). Labelled-pair coverage at 30x whatshap is
also reported as the mean, minimum and maximum over every choice of 1, 2
and 4 siblings.

For each child, depth (30, 20, 10), phaser (whatshap, LongPhase), the
cached scored pairs of `stratified_score.grouped` form rows grouped as in
R3 (`indel`, else the stratum of `stratified_gate.strata`). Odd chromosomes
calibrate per group on eligible pairs with label-based errors (the phaser's
relation against the sibling relation); `learn_then_test_by_group` at
`beta` 3e-4, 1e-3, 3e-3, `delta` 0.05, 100 kb units. Even chromosomes are
scored against the Pedigree truth on eligible pairs (the statement) and on
all pairs (for the trusted share of all pairs). Each child is its own
calibration: no pooling across children.

Label error is the share of eligible pairs whose sibling relation differs
from the truthset's phased relation of the child.

## Gate

Fixed before looking, for `exact` genotypes and `all` siblings:

- **(a) labels are right:** in each of the six depth and phaser cells,
  pooled over the five children, at most 1e-4 of eligible pairs (odd and
  even chromosomes) have a wrong label. Pooled per cell rather than per
  child because each child is expected to have of the order of 10^6
  eligible pairs (an estimate made before the run), so a per-child
  bound would be decided by a handful of pairs, and per cell rather than
  once because the pairs depend on the phaser's blocks;
- **(b) valid where it claims:** scored on eligible even-chromosome pairs,
  no violating cell at the Bonferroni level over all cells (5 children, 3
  depths, 2 phasers, 3 betas, 6 groups: 540 cells, level 9.3e-5);
- **(c) worth having:** `easy` at `beta = 1e-3` certifies in all 30 child,
  depth and phaser cells.

Q1 passes when (a), (b) and (c) hold. The output also gives the trusted
share of all even-chromosome `easy` pairs (trusted eligible pairs over all
pairs, abstaining on ineligible ones) and the mean eligible share of pairs.
No threshold is set on either; R3's `easy` share was 0.774 with parents.

## What a pass would and would not show

A pass shows that sibling genotypes, if exact, label enough of a child's
pairs, correctly enough, to certify the `easy` group per child at
`beta = 1e-3`. It is optimistic in two ways. The genotypes are the
truthset's, and the truthset's phase was built from the same family's
inheritance, so label and truth share information. The perturbed versions
show how label error scales with genotype error; the simulation predicts
about twice the error rate per sibling, so labels at 1e-4 need genotype
error near 5e-5 at the labelled sites, and real genotype errors are not
independent of the phasers' errors (R3's post hoc: replicate disagreement
in the parents' calls marks the sites where phasers fail). Four siblings
sequenced together is also more than most batches have; the `one` and `two`
sets show what fewer give.

## Amendment 1, 18 October: segments from opposite homozygotes instead of the HMM

**Seen before writing this.** The first run (`2026-10-17_q1_relative_labels.txt`)
failed (a) with label error 2.7e-2 and certified no `easy` cell. A diagnostic
that used the parents' transmission only to explain the failure
(NA12879 and NA12881, 21 chromosomes) found that the fitted HMM called 61% of
the sites in truly IBD0 segments IBD1 with P(IBD1) >= 0.99; the EM had moved
its IBD0 state to a pattern with almost no sites where both are heterozygous.
In truly IBD0 segments 7.7% of informative sites are opposite homozygotes
(IBS0), never in IBD1 or IBD2. With IBD1 known from the parents, the sibling
rule labelled 502,761 consecutive pairs of NA12879 with 0 wrong. The rule is
sound; the segmentation failed.

**Change.** A site is labelled only when two rules agree: the first run's
HMM rule (P(IBD1) >= 0.99 over the run, 100 kb margins) and a window rule
with no fitted model. With `W` = 250 kb, the window rule labels a site where
A is heterozygous and B homozygous when the sites within `W` on both sides
(both genotypes called, not both 0/0) hold no IBS0 site and at least 25% of
them are IBS1 (one heterozygous, one homozygous). Two labelled sites relate
only when they lie in one HMM run and no IBS0 site lies between them.
Reasons, from the diagnostic's IBD0 rate and not tuned on any label: an
IBD0 stretch in these data has an IBS0 site about every 13 informative
sites (some 8 kb), so a 250 kb window inside one holds none with
probability far below 1e-10; IBD2 has no IBS1 sites with exact genotypes,
and IBD1 has about half. The window rule alone is not enough: on the
module's simulated sibships, which draw mostly rare alleles, it labelled
17 of 8,763 pairs wrongly (short IBD0 stretches without IBS0), where the
HMM labelled 0 of 8,645; requiring both gave 0 wrong of 8,645. This was
seen on simulation only, before the rule ran on any Pedigree data. The
`--labeller hmm` option reproduces the first output. The gate, its levels
and everything else are unchanged; the rerun writes
`benchmarks/results/runs/2026-10-18_q1_relative_labels.txt`.

## Code

`benchmarks/q1_relative_labels.py`; `tests/test_q1_relative_labels.py`
checks the genotype classes, the transitions, the forward-backward
posterior against enumeration of all state paths, the run and margin rule,
the agreement rule, and the labels on simulated sibships, with exact and
with noisy genotypes. Output:
`benchmarks/results/runs/2026-10-17_q1_relative_labels.txt`.
