# Pre-registered plan for the genome-wide HG002 run

Written on 27 September 2026, before any HG002 read outside chromosome 20
was downloaded or read. The commit that adds this file fixes that. The
chromosome-20 run (`real-data-plan.md`) is the pilot: it shaped this plan,
and chromosome 20 is excluded from every confirmatory test below. Other
material read before writing: the GIAB v4.2.1 benchmark VCFs and BEDs of
HG002, HG003 and HG004 and the hifiasm phase transfer (no reads), and the
simulation results in the development notebook.

Changes after the first non-pilot read is processed go in the amendments
section, dated and committed before the step they concern. The registered
analysis is always reported, next to any amended one.

## Amendments

1. **27 September 2026, after the run finished and before any phasing was
   scored.** The plan leaves some details of the scoring open. They are fixed
   here, before `benchmarks/score_genome.py` first runs.
   - *Pooled* means over the concatenation of the confirmatory chromosomes:
     switch error is total switches over total pairs; block N50 and QAN50
     are taken over the pooled list of blocks; the objective ratio is total
     objective over total certified optimum, with the worst chromosome
     reported beside it; runtime is summed and peak memory is the maximum.
   - The bootstrap windows are cut per chromosome, as in `bootstrap.py`, and
     the window lists are concatenated before resampling.
   - A sign test is *in X's favour* when X is wrong at fewer of the
     discordant pairs than the other pipeline. In G6, the repair *has fewer
     switch errors* in this sense, on the common pairs.
   - G7's premise is a sign test of `loop` against `optimal` at level 0.05/8
     (caps 3 and 5 at four depths) in `loop`'s favour. The share kept by
     `shuffled` is (rate of `optimal` minus rate of `shuffled`) over (rate of
     `optimal` minus rate of `loop`), each pipeline on its own pairs, as in
     the pilot's H4, with the window bootstrap.
   - G3 and G6 register no refutation rule; a result that is not supported
     is reported as such.
   - G8 and G9 use the boundaries between index-adjacent variants that the
     selected reads span and that the filtered truth phases on both sides.
     The expected calibration error uses ten bins of width 0.1 on the
     predicted probability, weighted by boundary count.
   - `whatshap phase` and HapCUT2 enter the metric tables against every truth
     and no test.
2. **28 September 2026, after the registered analysis was scored
   (`a49444e`) and before any phasing was compared with the truth below.**
   A post hoc re-score against GIAB v5.0q, prompted by the literature
   re-check (prior-art matrix, section H), not by the results.
   - *Truth.* `HG002_GRCh38_v5.0q_smvar.vcf.gz` inside
     `HG002_GRCh38_v5.0q_smvar.benchmark.bed`, built from the T2T-HG002
     Q100 v1.1 assembly with dipcall. A variant of the run's set (the
     v4.2.1 heterozygous SNVs) gets the first allele of a `0|1` or `1|0`
     call at the same position, reference and alternate allele. Every
     other variant is unphased: outside the v5.0q BED, absent from the
     v5.0q VCF, or not heterozygous there. The haplotype order is the
     assembly's and is one phase per chromosome; switch error does not
     depend on which haplotype is first.
   - *Analysis.* Exactly the registered one, with this truth in place of
     the filtered transmission truth: the same tables, tests, levels and
     decision rules. The recalibration table is not refitted. Every
     decision under this truth is reported as post hoc, beside the
     registered decision, which it does not replace.
   - *Also reported.* Per chromosome, how many v4.2.1 heterozygous SNVs the
     v5.0q BED covers, how many v5.0q calls differently, and how often the
     v5.0q and transmission relations agree between consecutive sites both
     phase.

3. **30 September 2026: a correction, and a re-run of every library
   pipeline.**
   - *The artifact.* The library reported an allele for every variant a
     selected read's span crosses, even when no selected read calls it.
     Such alleles are arbitrary. Against v5.0q they are wrong at 50 per
     cent: 140 of 278 boundaries under WhatsHap's selection at cap 3 and
     full depth, and 8 of 18 under the exact one.
   - *The fix.* `FragmentMatrix.blocks()` now leaves such variants unphased
     (`33b92bb`). Coverage, which the cap counts, is unchanged, and so is
     every selection.
   - *What was affected.* Every library pipeline's reported phasing, and so
     every switch-error, pair and N50 figure in the registered scoring and
     in the v5.0q re-score. The objective ratios of G1 are not affected. Nor
     are the external tools, which never phase an uncalled variant.
   - *The re-run.* Every `lib` and `bq` job is re-run with the fix, and the
     registered analysis and amendment 2 are scored again. The scores from
     before the fix stay in the repository, labelled as affected.

## Why a second run

The chromosome-20 run could not decide its accuracy hypotheses. The
35-fold HiFi reads left WhatsHap's selection 6 excess switch errors at cap 3
and 3 at cap 5, and exact selection none, among about 35,000 pairs, which a
1 Mb bootstrap on one chromosome cannot separate. Two more problems appeared:
- the ratio-of-rates test misbehaves when two pipelines share their errors;
- the transmission truth had an error, confirmed from the father's own reads.

This plan fixes all three:
- power, from 21 chromosomes and four depths;
- a paired test per boundary (`phasebudget.paired`);
- a truth filter fixed before any read is seen.

## Data

**Reads.** The HG002 PacBio Sequel II CCS 11 kb GRCh38 alignment chosen
under `real-data-plan.md` amendment 1, streamed one chromosome at a time
through its index with pysam's bundled samtools.

**Chromosomes.** Autosomes 1 to 22 except 20: 21 chromosomes, all
confirmatory. Chromosome 20 is reported as the pilot and enters no test.

**Depths.** The full depth of each chromosome's slice, and seeded read
subsamples targeting 20, 15 and 10-fold (`samtools view --subsample F
--subsample-seed 0`). `F` is the target divided by the slice's depth, which
is measured the same way as in amendment 1: mapped records times the mean
aligned length of a 2,000-read sample, over the chromosome length. A slice
already below a target skips that target, and the skip is reported.

**Variants.** HG002's biallelic heterozygous SNVs inside its v4.2.1 BED, on
each chromosome.

**Fragments.** Two extractions per depth, both with WhatsHap 2.8's reader at
mapping quality 20:
- `realign`: default realignment against the reference, quality 30 per
  allele, converted to log-odds (`real-data-plan.md` amendment 3);
- `basequality`: alleles from the alignment, each weighted by the read's
  base quality at the SNV, recalibrated as below.

## Truth

**Primary: filtered transmission truth.** The transmission truth
(`truth-sets.md`), built per chromosome from the v4.2.1 calls, minus every
site where the hifiasm phase transfer's PATMAT genotype disagrees with it.
The filter is fixed here because of the chromosome-20 finding: the two sites
where every phasing disagreed with the truth were sites where hifiasm
disagreed too, and HG003's own reads showed the truth's premise was false.
Genome-wide, 2,420 of 1,686,093 compared sites disagree (0.14 per cent).
The unfiltered truth is reported beside every result.

**Secondary.** The hifiasm phase transfer alone, reported and supporting no
hypothesis.

## Caps

3, 5, 10 and 15, as before.

## Pipelines

Library pipelines on the `realign` fragments, log-odds weights:

| name | selection | phasing | caps |
|---|---|---|---|
| `whatshap` | WhatsHap 2.8's compiled `readselection` | `CoverageDP` | all |
| `optimal` | min-cost flow | `CoverageDP` | all |
| `merge_optimal` | HapCHAT merge at 10^6, then `optimal` | `CoverageDP` | all |
| `merge_optimal_1e3` | merge at 10^3, then `optimal` | `CoverageDP` | all |
| `two_level` | `optimal` | `CoverageDP`, then block orientation by all reads (`two-level-phasing.md`) | all |
| `whatshap_two_level` | as `whatshap` | the same orientation | all |
| `merge_two_level_1e3` | as `merge_optimal_1e3` | the same orientation, by the original reads | all |
| `connected` | `select_connected`, k 1, penalty one median read weight | `CoverageDP` | all |
| `loop` | stitched loop as registered for chromosome 20 | its own | 3, 5 |
| `loop_oriented` | the same with `stitched_oriented` acceptance | its own | 3, 5 |
| `merge_loop_1e3` | merge at 10^3, then `loop` | its own | 3, 5 |
| `shuffled` | `loop` with the shuffled signal | its own | 3, 5 |

The loop-type pipelines run at caps 3 and 5 only. At 10 and 15 on the pilot
they stopped after one round and equalled static selection, at a large cost
in time.

On the `basequality` fragments, with recalibrated weights, `optimal` and
`whatshap` run at every cap, and `optimal` at cap 5 records the boundary
signals.

External tools, on the same BAM slices and VCF:
- `whatshap phase` 2.8, defaults (internal downsampling 15);
- HapCUT2 at commit `9a10aba`, `extractHAIRS --pacbio 1` with the reference,
  then `HAPCUT2`.

## Recalibration of base qualities

On the full-depth `basequality` fragments of the odd chromosomes (1, 3, ...,
21), each read is assigned to the truth haplotype it agrees with at the
majority of its truth-phased sites (reads with fewer than 3 such sites, or a
tie, are skipped). Its alleles at truth-phased sites are then counted as
right or wrong per base quality value.
- Values are pooled from the top down until every bin has at least 50
  wrong alleles.
- The error rate of a bin is (wrong + 1) / (alleles + 2).
- The weight of an allele is the log-odds of that rate.

The table is fitted on odd chromosomes only and applied to all. Calibration
is tested on the even chromosomes (2, 4, ..., 22 except 20), which never
touched the fit.

## Statistics

**Paired test.** For two pipelines at one depth and cap, the common pairs
are pairs of consecutive variants that the filtered truth and both
pipelines phase, inside one block of each (`phasebudget.paired`). Pairs
where only one pipeline is wrong are counted over all confirmatory
chromosomes and tested with the exact two-sided sign test. The effect size
is the difference in error rate on the common pairs, per 10,000 pairs, with
a 95 per cent percentile interval from 2,000 resamples of 1 Mb windows
pooled across chromosomes, seed 0.

**Multiple comparisons.** Within each hypothesis the significance level is
0.05 divided by the number of tests the hypothesis makes (Bonferroni).

**Per-pipeline metrics.** Switch error (count and denominator), flips,
phased fraction, block N50 and QAN50 in base pairs, the objective against
the certified optimum where reads are not merged, runtime and peak memory,
per chromosome and pooled.

## Power

From the pilot, WhatsHap's selection had 6 discordant pairs against
`optimal` at cap 3 and 3 at cap 5, at 35-fold, all in `optimal`'s favour.
The confirmatory chromosomes hold about 45 times chromosome 20's truth-phased
sites, so about 270 and 135 discordant pairs are expected at full depth, and
more at lower depth.

Let `n` be the number of discordant pairs, `pi` the share in `optimal`'s
favour, and `z = (pi - 1/2) / sqrt(1 / (4n))` the test statistic, compared
with 2.73 (8 tests at 0.05):

| `n` | `pi` | `z` | approximate power |
|---:|---:|---:|---:|
| 270 | 0.60 | 3.29 | 0.71 |
| 270 | 0.70 | 6.57 | above 0.99 |
| 135 | 0.70 | 4.65 | 0.97 |

The pilot's split was 6 to 0, so a share near 0.7 or above is the
expectation. The run has adequate power for effects of that size and
limited power for a 60 to 40 split at 35-fold, which the lower depths
should compensate for.

## Hypotheses and decision rules

**G1, deployed gap.** Pooled over confirmatory chromosomes at full depth and
cap 15, WhatsHap's selection reaches less than 0.95 of the certified
optimum. Descriptive; supported if below 0.95.

**G2, selection matters on real reads at small caps.** `whatshap` against
`optimal`, caps 3 and 5, each of the four depths: 8 sign tests, level
0.05/8. Supported if all tests at 10-fold and 15-fold are significant in
`optimal`'s favour. Refuted if no test at any depth is significant.

**G3, merging.** `merge_optimal_1e3` against `optimal`, caps 3 and 5, four
depths, level 0.05/8. Supported if significant in the merge's favour at all
tests at 10-fold and 15-fold.

**G4, two-level contiguity.** At cap 3 and each depth, `two_level` against
`optimal`. Supported if block N50 rises by at least 25 per cent at every
depth and no sign test is significant in `optimal`'s favour at level
0.05/4.

**G5, connected selection.** The same rule as G4, for `connected` against
`optimal`.

**G6, the loop and its repair.** `loop_oriented` against `loop` at cap 3,
four depths, level 0.05/4. Supported if the repair has fewer switch errors
and is significant at 10-fold and 15-fold. `loop` against `optimal` is
reported, two-sided, as the test of whether reselection helps real reads.

**G7, mechanism (H4 again).** Only if `loop` has significantly fewer errors
than `optimal` at some cap and depth: there, `shuffled` keeps at least 80
per cent of the reduction. Otherwise G7 is reported as without its premise.

**G8, calibration.** On the even confirmatory chromosomes, full depth, cap
5, `optimal` on recalibrated `basequality` weights. Calibrated if the
bootstrap interval of observed wrong calls over predicted ones contains 1,
and the expected calibration error over ten equal-width bins is below 0.01.

**G9, margins as a signal.** Same data as G8: the AUC of the posterior margin
for a wrong call at a spanned, truth-phased boundary is above 0.6.

**G10, depth.** Descriptive: the effect of G2 at each depth, to show where
along depth selection starts to matter.

## Safety and resources

- Each job runs at `nice 19` and idle I/O priority, under a
  `systemd-run --user --scope` memory ceiling (4 GB for library pipelines,
  6 GB for WhatsHap and HapCUT2).
- A scheduler starts a job only when at least 4 GB of memory is available,
  and runs at most 8 at once.
- It starts a download only when at least 20 GB of disk is free, and at
  most 2 at a time.
- BAM slices are deleted once every job that needs them has finished;
  fragments, phasings and logs are kept.

## Not in this plan

The ONT run and the trio run (research plan T5 and T3) will each be
registered separately, before their reads are touched.
