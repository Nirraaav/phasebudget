# First gates of Y1 and Y2

Written on 7 October 2026, before either gate was run. The directions are in
the October directions memo (development repository). Each gate can stop its direction.
Neither uses a Platinum Pedigree member other than NA12878, so neither touches
the registered test genomes of X6.

## Y1: are transmission labels clean enough to calibrate on?

**Setting.** Y1 would calibrate the per-relation threshold on relations that
transmission phases in a family, instead of on a benchmark truth. Browning
and Browning 2022 (AJHG 109:1016-1025, `abstract`) show that genotype error
makes trio-based phase errors look 2.4 times more frequent than they are in
UK Biobank. A label error rate near the target `beta` would make calibration
meaningless, so the label error rate is measured first.

**Labels.** `phasebudget.truth.phase_by_transmission` on HG002 (child), HG003
(father) and HG004 (mother):
- child calls: GIAB v5.0q HG002, biallelic heterozygous SNVs inside the
  v5.0q BED;
- parent calls: GIAB v4.2.1 HG003 and HG004, with their benchmark BEDs (a
  parent site outside its BED is not used);
- the 21 confirmatory autosomes.

Both inputs are curated truths, so genotype error is far below a variant
caller's. The measured rate is therefore a lower bound on what labels from
called genotypes would carry.

**Measure.** Among consecutive labelled sites (consecutive in the list of
labelled sites), the relation (same or opposite parent) from the labels is
compared with the relation from v5.0q's phased genotypes, where v5.0q phases
both in one phase set. Reported: labelled relations, disagreements, rate with
a Clopper-Pearson 95 per cent interval, the share of v5.0q heterozygous SNVs
labelled, and both per stratum (`benchmarks/stratified_gate.py`'s strata).

**Rule.** Y1 continues if the disagreement rate's upper 95 per cent bound is
at most 1e-4 overall and in the easy stratum. Otherwise the labels cannot
support a `beta` of 1e-3, and Y1 stops.

**Known limit.** v5.0q is built from the Q100 assembly, whose haplotypes were
assigned to parents with parental data. Agreement between the two is then
partly by construction; disagreements are still errors of one of them.

## Y2: does the certified relation reach clinical genes?

**Genes.** Symbols in ClinVar's `gene_specific_summary.txt` dated 4 October
2026 (md5 `9761c2d9b82a96d8601b6129d607b312`, as published beside it) with
at least 10 alleles reported pathogenic or likely pathogenic, matched by
`gene_name` to protein-coding `gene` records of GENCODE v50 basic. 3,893
symbols pass the ClinVar filter before matching.

**Relations.** As in `docs/design/chain-relations.md`: pairs of truth-phased
heterozygous SNVs in one block and one listed gene. A link is trusted when
its score reaches its stratum's threshold from
`learn_then_test_by_group` on the registered HG002 ONT sup calibration rows
(as in `stratified_gate.py --cross`); a link in a stratum with no certified
threshold is untrusted. A relation is trusted when all its links are.

**Data.** 30x, `whatshap phase` and LongPhase, HG002 ONT hac (v5.0q) and
HG001 ONT sup (Platinum 2017, as in track D), `beta` 1e-3. Data seen before.

**Measure.** The share of relations trusted, wrong trusted relations and their
rate, the mean bound `beta * m`, and the 20 genes with the lowest trusted
share among genes with at least 100 relations.

**Rule.** Y2 continues if, on HG002 hac for both phasers, at least half of
the relations in listed genes are trusted and the wrong trusted rate is at
most the mean bound. Otherwise the certified relation does not reach enough
of the clinical population to carry a paper.

## Y1, step 2: calibrating on labels (stated 7 October, before running)

**Correction to step 1's framing.** Step 1's rule counted label errors per
relation. That is not what decides whether labels can calibrate. A label
error on a relation the phaser got right makes the relation look wrong, which
only makes the threshold stricter. Only a label error that hides a real
phaser error (label and phaser wrong together) makes calibration unsafe, and
its rate is about the label error times the phaser error. So the question is
answered directly: calibrate on labels and compare with calibrating on the
benchmark truth.

**Labels.** Three sources for the loss of each HG002 calibration pair:
- `truth`: v5.0q, as registered;
- `labels-truth`: transmission from the curated genotypes of step 1;
- `labels-called`: transmission from the Clair3 calls of the ONT
  `giab_2023.05` `wf_human_variation` runs, `hg002_sup_all`, `hg003_sup_all`
  and `hg004_sup_all` (`*.wf_snp.vcf.gz`), PASS records only. There is no
  callable-region BED, so a site absent from a parent's VCF counts as
  homozygous reference, as a lab without one would assume.

The phaser's relation comes from the row's loss against v5.0q, as in track C.

**Calibration.** HG002 ONT sup rows of the registered calibration, at 30, 20
and 10x, `whatshap phase` and LongPhase, restricted to pairs whose two SNVs
carry a label from both label sources, so the three calibrations see the
same pairs. Pooled and per-stratum Learn then Test, `beta` in {3e-4, 1e-3,
3e-3}, `delta` 0.05.

**Tests.** Each threshold is applied to all pairs of HG002 ONT hac (v5.0q)
and of HG001 ONT sup (Platinum Pedigree, track C's adjudicable pairs).
Reported per cell: violations at p < 0.05, trusted share, and wrong trusted
pairs.

**Rule.** Y1 continues if, with `labels-called`, (a) no test cell violates
that does not also violate with `truth`, and (b) on HG002 hac at `beta =
1e-3` and 30x the trusted share is within 0.01 of `truth`'s for both phasers.

## Y2, step 2: does withholding fall on independently known difficult genes? (stated 7 October, before running)

**Why a second step.** Step 1 printed the 20 clinical genes with the lowest
trusted share before any list of difficult genes was chosen, so any list
picked after that is not blind. The list used here is external and fixed:
GIAB's challenging medically relevant genes (CMRG v1.00,
`GRCh38_CMRG_benchmark_gene_coordinates.bed`, 273 genes, md5
`ab47ef868d48846c573b02932906d4d1`), matched by symbol.

**Hypothesis.** Among the clinical genes of step 1 with at least 100
relations, at `beta = 1e-3`, 30x, HG002 hac, per-stratum thresholds:
- the median trusted share of CMRG genes is below that of the other genes
  (one-sided Mann-Whitney U, level 0.05), for each phaser;
- CMRG genes are over-represented among genes with trusted share below 0.5
  (one-sided Fisher exact test, level 0.05), for each phaser.

**Limits stated now.** The per-stratum thresholds certify nothing in
segmental duplications, and many CMRG genes were chosen for overlapping
them, so part of any association is mechanical. The data were seen at gene
level in step 1. This step is exploratory; the confirmatory version belongs
in the registration of X6, on the test genomes.

## Y1, step 3: a genotype-quality filter on the labels (stated 7 October, before running)

**Why.** Step 2 found that labels from called genotypes hide no phaser errors
but invent 17 to 25 per cell, which blocks `beta = 3e-4` and costs
usefulness per stratum. Low genotype quality is the obvious source: about 1.5
million of HG003's 5.3 million PASS calls have GQ below 20.

**Filter `g`.** A child site is labelled only if its call has GQ >= `g`. A
parent call with GQ < `g` makes the site unusable for that parent (it is not
treated as homozygous reference, as an absent call is).

**Choosing `g` on the odd chromosomes only.** For `g` in {0, 10, 20}, count on
chromosomes 1, 3, ..., 21 the labelled sites and the disagreements with
v5.0q among consecutive labelled sites (step 1's measure, with called
labels). Choose the `g` with the fewest disagreements among those that keep
at least 80 per cent of the labelled sites of `g = 0`; ties go to the
smaller `g`.

**Evaluating on the even chromosomes only.** Step 2's procedure and rule,
with calibration rows restricted to chromosomes 2, 4, ..., 22, comparing
`labels-called` at the chosen `g` against `truth` on the same rows. Tests as
in step 2. Y1 continues if (a) and (b) of step 2 both hold, with (b) required
in both modes, as the step 2 script did.

## Y2, step 3: genes with a duplicated pseudogene copy (stated 7 October, before running)

**List.** A protein-coding gene `X` of GENCODE v50 basic is *duplicated* when
the same annotation has a `gene` record of type `unprocessed_pseudogene` or
`transcribed_unprocessed_pseudogene` named `X` followed by `P` and digits
(HGNC's convention for pseudogenes, e.g. `STRCP1`). Processed pseudogenes
are excluded: they lack introns, and reads spanning an intron place
unambiguously. The list depends only on the annotation.

**Not blind.** The genes with the lowest trusted share were printed in step 1
before this list was chosen, and the step 2 result was known. The choice of
an annotation-based list was informed by both. So this step is descriptive;
the test belongs in X6's registration, on genomes not yet seen.

**Measure.** As step 2, with `duplicated` in place of CMRG: medians,
Mann-Whitney U (one-sided), the share below 0.5, Fisher (one-sided), per
dataset and phaser, `beta = 1e-3`, 30x.

The preferred external list is Mandelker et al. 2016 (*Genet Med*
18:1282-1289, 193 genes, `recalled` until read); its supplement was not
obtained.
