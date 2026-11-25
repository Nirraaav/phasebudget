# R1: phase-truth error rates without a gold standard

Written 16 October 2026, before any fit ran. Development gate, not
registered. Direction R1 of the October directions memo (development repository).

## What it estimates

For a genome, a group `g` and a set of raters, each relation between two
consecutive heterozygous variants has a true relative phase `Z` (same or
opposite haplotype). Rater `r` reports `Z` wrong with probability `e[r, g]`,
independently of the other raters given `Z`. `pi[g]` is the share of
relations with `Z = same`. The model is fitted by expectation maximization
on the counts of agreement patterns, per group. It is identified when at
least three raters have `e < 0.5` (Dawid and Skene 1979, `recalled`).

The estimate of interest is `e[r, g]` for each truth `r`. A truth's error
here is an error of its relative phase between two consecutive variants
both raters phase; genotype errors enter only through the relations they
break.

## Raters

| genome | truths | read-based raters |
|---|---|---|
| NA12878 | Platinum Pedigree (`pedigree_truth`), Platinum Genomes 2017 (`platinum_truth`), 1000 Genomes transmission (`kg_transmission`) | whatshap 30x on ONT (`hg001_v5`); whatshap 30x on Platinum Pedigree HiFi (`phifi_NA12878`) |
| HG002 | v5.0q (`v5_truth`), GIAB trio transmission, hifiasm phase transfer | whatshap 30x on ONT (`cal_v5`); whatshap 30x on HPRC HiFi (`hprc_HG002`) |

Variants are the platform's phaser input (`hg001_v5`, `cal_v5`): SNVs and
indels up to 50 bases. A relation enters a fit when every rater of that fit
gives both variants a phase, and a read-based rater has them in one block.
Groups are those of `stratified_gate.strata`, with `indel` for a relation
touching an indel.

## Fits

- **Main fit (4 raters, no HiFi):** the three truths and ONT whatshap.
- **Check rater:** HiFi whatshap is left out of the main fit and used only
  to check the main fit's calls (gate part c).
- **Leave-one-out:** the main fit without each rater in turn, to show how
  much each estimate rests on one rater.
- **Fit diagnostic:** the likelihood-ratio statistic of the pattern counts
  against the fitted model, and per pair of raters the observed number of
  relations both get wrong (by the posterior) against the number expected
  under independence.

Uncertainty: 100 bootstrap resamples of 1 Mb blocks within each group.

## Gate

Fixed before looking. On NA12878:

- **(a) stable:** every truth with `e >= 1e-4` in a group has odd- and
  even-chromosome estimates within a factor 2 of each other;
- **(b) a measurable truth error:** in `lowmap`, `segdup` or `tandem`, at
  least one truth has `e >= 1e-3` with bootstrap 2.5% quantile `>= 1e-4`;
- **(c) the reads agree:** among relations in those groups where the main
  fit gives a truth's call a posterior of being wrong of at least 0.99, and
  that HiFi whatshap phases (at least 20 such relations), HiFi agrees with
  the fit's call in at least 75%.

R1 passes when (a), (b) and (c) hold. HG002 is reported, not gated.

## Limits

The NA12878 truths are all inheritance-based within CEPH 1463, and two of
them rest on Illumina reads, so their errors may be dependent; a dependence
that makes two truths err together makes the fit call the third one wrong.
The read-based raters fail together in repeats. Relations need every rater,
so the fitted population is the one all of them phase, which leaves out
the hardest relations.

## Code

`benchmarks/r1_truth_audit.py`; `tests/test_r1_truth_audit.py` checks the
fit against simulated raters with known error rates.
