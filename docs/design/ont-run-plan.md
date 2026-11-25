# Pre-registered plan for the genome-wide HG002 ONT run

Written on 28 September 2026, before any ONT read was downloaded or read.
The commit that adds this file fixes that. Material read before writing:
- the HiFi genome-wide run and its scoring, registered (`a49444e`) and
  against GIAB v5.0q (amendment 2 of `genome-run-plan.md`);
- the D1 pilot (development notebook, 27 September);
- the file listing of the ONT data below, which contains no reads.

This is track T5's ONT arm, and D3 of the research plan (development repository) section 13. HiFi
could not say whether read selection matters on noisier reads. This run
asks that question with the same variants, truths, pipelines and
statistics, so the two platforms can be compared directly.

Changes after the first ONT read is processed go in an amendments section,
dated and committed before the step they concern.

## Amendments

1. **28 September 2026, after the chromosome-20 smoke test and before any
   confirmatory chromosome was downloaded.** The smoke test found two
   problems. Neither involved a phasing being compared with a truth.
   - *No read groups.* The BAM has no `@RG` lines, so WhatsHap's reader
     matched no read to sample HG002: 0 reads, 0 fragments. The fragment
     extraction and `whatshap phase` now run with read groups ignored,
     which WhatsHap documents as "assume all reads come from the same
     sample". HapCUT2 and LongPhase do not use read groups.
   - *Depth.* The registered depth counts every mapped record in the index,
     and this BAM holds many secondary alignments. On chromosome 20 that
     gave 87.2x for a release named 60x. Slice depth is now the summed
     aligned length of every record that is not secondary, over the
     chromosome length, in one pass over the slice. The index-based figure
     is still written beside it, and the subsample fractions use the new
     depth.

2. **30 September 2026, before any confirmatory chromosome was scored.**
   Genome-run amendment 3 applies here too: variants no selected read calls
   are left unphased (`33b92bb`). Every library job that finished before
   the fix is re-run. The external tools and the fragment extractions are
   not affected.

## Data

**Reads.** Oxford Nanopore's open-data release `giab_2023.05`: HG002,
R10.4.1, basecalled in super-accuracy mode, aligned to GRCh38 and
haplotagged by the release's own workflow:
`https://ont-open-data.s3.amazonaws.com/giab_2023.05/analysis/variant_calling/hg002_sup_60x/hg002.haplotagged.bam`
(205,484,249,144 bytes, with a `.bai`). The haplotag fields (`HP`, `PS`) are
never read. Reads are streamed one chromosome at a time through the index
with pysam's bundled samtools, as in the HiFi run.

**Chromosomes.** Autosomes 1 to 22 except 20: the 21 confirmatory
chromosomes of the HiFi run. Chromosome 20 runs first as a smoke test of
the pipeline and enters no test.

**Depths.** The full depth of each chromosome's slice, and seeded
subsamples targeting 30, 20 and 10-fold (`--subsample-seed 0`), with slice
depth measured as in `genome-run-plan.md`. A slice already below a target
skips it, and the skip is reported.

**Variants.** The HiFi run's set: HG002's v4.2.1 biallelic heterozygous
SNVs inside the v4.2.1 BED, per chromosome. The same set on both platforms
makes their pooled rates comparable.

**Fragments.** One extraction per depth, WhatsHap 2.8's reader at mapping
quality 20 with realignment against the reference, converted to log-odds
(`real-data-plan.md` amendment 3). No base-quality extraction: T4 is not
repeated on ONT.

## Truth

**Primary: GIAB v5.0q.** The phased calls of
`HG002_GRCh38_v5.0q_smvar.vcf.gz` inside its benchmark BED, assigned to the
variant set as in amendment 2 of `genome-run-plan.md`. On the HiFi variant
set, v5.0q and the transmission truth give the same relation between
consecutive sites at all but 2 of 1.63 million pairs. v5.0q phases 1.97
million sites against 1.63 million, so it is the primary truth here.

**Secondary.** The filtered transmission truth of the HiFi run, reported
beside every test and supporting no hypothesis.

## Caps and pipelines

Caps 3, 5, 10 and 15. Library pipelines on the log-odds fragments, each
phased with `CoverageDP`:

| name | selection | phasing |
|---|---|---|
| `whatshap` | WhatsHap 2.8's compiled `readselection` | `CoverageDP` |
| `optimal` | min-cost flow | `CoverageDP` |
| `merge_optimal_1e3` | HapCHAT merge at 10^3, then `optimal` | `CoverageDP` |
| `two_level` | `optimal` | `CoverageDP`, then block orientation by all reads |
| `connected` | `select_connected`, k 1, penalty one median read weight | `CoverageDP` |
| `error_optimal_1e2` | `select_error_optimal`, switch-model costs, allele error 0.00117, scale 1e2 (the D1 pilot's tuned arm) | `CoverageDP` |

The adaptive loops are not run. They were negative on HiFi (G6, G7) and
cost about 20 times the static pipelines.

**External tools**, on the same BAM slices and VCF:
- `whatshap phase` 2.8, defaults;
- HapCUT2 at commit `9a10aba`, `extractHAIRS --ont 1` with the reference,
  then `HAPCUT2`;
- LongPhase 2.0.2, `phase --ont`, defaults, SNVs only.

HiPhase is for PacBio HiFi only and is not run.

## Statistics

As in `genome-run-plan.md` and its amendment 1:
- the paired exact sign test on common pairs;
- the difference per 10,000 common pairs, with a 1 Mb window bootstrap
  (2,000 resamples, seed 0) pooled across chromosomes;
- Bonferroni within each hypothesis;
- the same per-pipeline metrics, pooled.

## Hypotheses and decision rules

**O1, deployed gap.** At full depth and cap 15, WhatsHap's selection reaches
less than 0.95 of the certified optimum, pooled. Supported if below 0.95.

**O2, selection matters beyond cap 3 on ONT.** `whatshap` against
`optimal` at caps 3, 5 and 10, four depths: 12 sign tests at level 0.05/12.
- Supported if the test at cap 5 is significant in `optimal`'s favour at
  every depth. On HiFi none of the cap-5 tests were significant.
- Refuted if no test at cap 5 or 10 is significant.
- Cap 3 is reported and expected to favour `optimal`, as on HiFi.

**O3, merging.** `merge_optimal_1e3` against `optimal`, caps 3 and 5, four
depths, level 0.05/8. Supported if significant in the merge's favour at
cap 5 at every depth.

**O4, two-level contiguity.** At cap 3 and each depth, `two_level`
against `optimal`: block N50 rises by at least 25 per cent at every depth.
The sign test on common pairs cannot see the joins, so the rule has a
second part. The switch rate on the pairs `two_level` adds over `optimal`
must be at most twice `optimal`'s own switch rate. That rate is
(`two_level` switches minus `optimal` switches) over (`two_level` pairs
minus `optimal` pairs), pooled per depth. Supported if both parts hold at
every depth.

**O5, the D1 negative.** `error_optimal_1e2` against `optimal` at cap 3,
four depths, level 0.05/4. The HiFi pilot found it worse. Replicated if
significant in `optimal`'s favour at two or more depths. Reversed if
significant in its own favour at two or more depths.

**O6, capped against uncapped (descriptive).** Pooled switch error, flips,
phased fraction, block N50 and QAN50 of every library pipeline at cap 15,
beside `whatshap phase`, HapCUT2 and LongPhase, at each depth.

**O7, platforms (descriptive).** The effect of O2 at each cap and depth
beside the HiFi run's G10 table. The HiFi rows are re-scored against v5.0q,
so both platforms use one truth.

## Safety and resources

The scheduler, limits and clean-up of `genome-run-plan.md`: `nice 19`, idle
I/O, a memory ceiling per job, at most 2 downloads and 3 chromosomes of
reads on disk, and at least 20 GB free before a download. Reads are
deleted once every job that reads them has finished. The slices are about
twice the HiFi run's size, and the run should take two to four days.

## Not in this plan

Base-quality extraction and recalibration (T4), the adaptive loops, the
trio (T3), and the difficult regions outside the v4.2.1 BED.
