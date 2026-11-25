# Pre-registered plan for the real-data run

Written on 26 September 2026, before any HG002 read was downloaded,
inspected or processed. The commit that adds this file fixes that date. The
only real data touched so far is the GIAB benchmark VCFs and BEDs of
`docs/design/truth-sets.md`, which contain no reads.

The simulation results that this run tests are in the development notebook. Each
hypothesis below names the simulation number it expects to reproduce, so a
reader can tell a confirmation from a result found after looking.

## Amendments

Any change made after the first read is processed is appended here with its
date and reason. The analysis as registered is still run and reported next to
the amended one.

- **1, 26 September 2026, the read set.** Chosen by the registered rule
  before any read was processed. GIAB's `giab_data_indexes` lists three
  GRCh38 HiFi alignments of HG002. Their depths are not in the index, so they
  were estimated on chromosome 20 from each BAM's remote index and a sample
  of 2,000 alignments (mapping quality 20 or more, primary, from 30 to 40 Mb):
  mapped records times mean aligned length over the chromosome length. The
  record count includes secondary and supplementary records, so the
  estimates are upper bounds.

  | data set | chr20 records | mean aligned length | depth |
  |---|---:|---:|---:|
  | PacBio Sequel II CCS 11 kb | 202,617 | 10,992 | 34.6 |
  | PacBio CCS 15 kb + 20 kb, chemistry 2 | 268,343 | 13,764 | 57.3 |
  | PacBio Sequel CCS 15 kb | 161,710 | 12,895 | 32.4 |

  The rule selects the 11 kb Sequel II set, which its README describes as
  about 32-fold: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_SequelII_CCS_11kb/HG002_GRCh38/HG002_GRCh38.haplotag.10x.bam`,
  69,133,873,504 bytes, last modified 30 August 2019, MD5
  `80d4a45c9272e352981c2899267e63c5` by the `.md5` file beside it. The
  `giab_data_indexes` row for this file gives `0719678b2c25a9bf9bfc7987b2ca6193`,
  which is the MD5 it also gives for the hs37d5 BAM, so the index row is
  wrong. The index file (`.bai`, MD5 `296f199e7e9721635d6f5085c94afa57`)
  matches. Only chromosome 20 was downloaded, through the index, so the whole
  file's MD5 could not be checked. The reads were aligned with pbmm2 to
  GRCh38_no_alt_analysis and carry WhatsHap `HP` tags from a 10x phasing,
  which nothing here reads. The reference is UCSC's hg38 `chr20.fa`, whose
  chromosome 20 has the BAM header's length, 64,444,167.

- **2, 26 September 2026, a second merge threshold, added as a secondary
  pipeline.** After registration, the milestone-4 sensitivity sweep
  (the development notebook) found HapCHAT's merge threshold of 10^6 too conservative on
  simulation: 10^3 lowered switch error and lengthened blocks on the tuning
  seeds, on held-out seeds and under unreported errors. `merge_optimal` and
  `merge_loop` stay as registered, with 10^6, and H3 is decided on them.
  `merge_optimal_1e3` and `merge_loop_1e3`, identical except for threshold
  10^3, are added and reported beside them, labelled as added after
  registration. This amendment was written before any fragment was
  extracted.

- **3, 26 September 2026, two corrections found while extracting fragments,
  before any pipeline was run.** First, the count of variants. The Variants
  paragraph says 46,344 sites inside the BED; 46,344 is the number of
  biallelic heterozygous SNVs in the whole chromosome-20 VCF, which is what
  `build_trio_truth.py` counted. Inside the BED, as the rule says, there are
  45,053, and those are the variants. Truth counts are recomputed on them
  when scoring. Second, the weight units. The Fragments paragraph says the
  WhatsHap qualities are used unscaled. WhatsHap's qualities are phred
  values, and with its default realignment every allele gets 30 (`realign`
  in `whatshap/variants.py` sets `base_qual_score = 30`); the library's
  weights are natural log-odds, the scale the posterior and the loop's rule
  `P(wrong) > 0.05` read. Read unscaled, 30 would state an error probability
  near `10^-13`, the loop would find no fragile boundary, and H5 would be
  decided by the units. The weights are therefore converted,
  `ln((1-e)/e)` with `e = 10^(-q/10)` (`extract_fragments.phred_to_log_odds`),
  which gives 6.91 for every allele. The minimum-cost pipelines choose the
  same reads and haplotypes either way, since every weight is scaled
  alike; the registered unscaled version is still run and reported for the
  pipelines it can change, `loop`, `merge_loop`, `shuffled` and H5 and H6.
  Because every weight is equal, the additive selection objective on HG002
  counts called alleles, the regime of the simulation control with
  `--quality-spread 0`.

## Data

**Reads.** One PacBio HiFi data set of HG002 aligned to GRCh38, chosen by
this rule: the HG002 HiFi alignment listed in GIAB's `giab_data_indexes`
repository with the highest mean depth that does not exceed 40, so that the
full-depth run needs selection and downsampling is not required. If no listed
HiFi alignment has depth 20 or more, the ONT alignment with the same rule is
used instead and the report says so. The chosen file's URL, size and MD5 are
recorded as the first amendment before fragment extraction starts.

**Region.** Chromosome 20, restricted to the HG002 v4.2.1 benchmark BED.

**Variants.** The biallelic heterozygous SNVs of the HG002 v4.2.1 benchmark
VCF inside that BED, 46,344 sites. Genotypes are taken from the benchmark and
trusted. Every site is heterozygous by construction, which is the
complementary model's assumption; indels and homozygous sites are excluded.

**Fragments.** Allele detection is WhatsHap 2.8's `ReadSetReader`, the reader
`whatshap phase` uses, with its default realignment and mapping-quality
threshold (20), on the GRCh38 reference GIAB used. Secondary and
supplementary alignments are dropped, as WhatsHap drops them. Each read's
alleles become a `Fragment`; its weights are the integer allele qualities the
reader reports, used unscaled. Reads calling one variant are kept in the
matrix for coverage accounting and are never selected, as on simulation.
`benchmarks/whatshap_bridge.py` already converts in the other direction; the
converter is tested on a hand-built BAM before it touches HG002.

**Depths.** Two runs: the full depth of the chosen file, and a 50 per cent
read subsample drawn with seed 0 before extraction.

## Truth

Primary: phasing by transmission from the parents' v4.2.1 calls,
`docs/design/truth-sets.md`, which phases 36,860 of the 46,344 sites and was
built without reads or a phasing tool. Switch and flip errors are counted
over consecutive variants of a predicted block that the truth phases,
`docs/design/metrics.md`.

Secondary, for the 9,484 sites the primary truth leaves unphased: the hifiasm
phase transfer, whose measured discordance with the primary truth is 13 of
36,834 sites. It is reported in a separate table with that figure beside it
and supports no hypothesis below.

## Caps

OQ-9 is decided by this rule. The DP's state limit is 22 active reads
(`dp.MAX_ACTIVE_READS`), and every selector keeps span coverage at or below
the cap, so any cap up to 22 fits. The caps run are 3, 5, 10 and 15: 15 is
WhatsHap's default, 5 is the per-member cap of a trio at that default, and 3
and 10 bracket them. If the full-depth coverage at the 95th percentile of
sites is below 15, cap 15 is still run and reported as unconstrained at most
sites.

## Pipelines

The same code and parameters as the simulation tables, with no tuning on
HG002:

| name | selection | phasing |
|---|---|---|
| `whatshap` | `select_whatshap`, the port that reproduces WhatsHap 2.8 | `CoverageDP` |
| `optimal` | `select_optimal`, min-cost flow | `CoverageDP` |
| `merge_optimal` | `hapchat.merge_reads`, then `select_optimal` | `CoverageDP` |
| `loop` | `adaptive_selection`, posterior signal, error probability 0.05, 4 rounds, boost 2, stitched acceptance, window 25 | its own |
| `merge_loop` | `merge_reads`, then `loop` | its own |
| `shuffled` | `loop` with the shuffled signal | its own |
| `whatshap phase` | WhatsHap 2.8 as released, default settings, same VCF and BAM | its own |

HapCUT2 is added as a second external tool if it installs without root
privileges; otherwise its absence is stated. This settles OQ-8 unless the
author overrules it before the run. Every selection is scored under
the base weights.

## Metrics

Switch error, flip errors, phased fraction, block N50 and QAN50 in base
pairs, wMEC of the reported haplotype against the full read matrix, the
additive selection objective over the certified optimum, runtime and peak
memory. Accuracy and objective value stay in separate columns.

**Uncertainty.** Chromosome 20 is one instance, so intervals come from a
block bootstrap over the variants: 1 Mb windows of the chromosome, resampled
with replacement 2,000 times with seed 0, recomputing each statistic on the
resample. A ratio of switch error rates is reported with its 95 per cent
percentile interval. The bootstrap is `phasebudget.bootstrap`, written and
tested on simulated instances before the real run (`tests/test_bootstrap.py`
checks its window counts against the metrics module). A pair belongs to the
window of its second variant. A resample on which a ratio is undefined, because
the denominator pipeline has no switches in it, is left out of the
percentiles, and the number of such resamples is reported.

## Hypotheses and decision rules

**H1, the deployed gap.** At cap 15 on the full-depth run, `whatshap` reaches
less than 0.95 of the certified additive optimum. Supported if its ratio is
below 0.95. Simulation, with the port, expects the value in the
chromosome-scale table of the development notebook.

**H2, where optimal selection matters.** At caps 3 and 5, `optimal` has lower
switch error than `whatshap`, with the bootstrap interval of the ratio
`whatshap / optimal` above 1. At cap 15 the interval contains 1. Supported if
both parts hold; refuted if the cap-3 and cap-5 intervals both contain 1.

**H3, the pipeline.** At caps 3 and 5, `merge_loop` has lower switch error
than `optimal`, with the interval of the ratio `optimal / merge_loop` above 1
at both caps. Simulation under unreported errors expects ratios of 3.6 to
4.4. Refuted if the interval contains 1 at both caps.

**H4, the mechanism.** `shuffled` keeps at least 80 per cent of the loop's
reduction, that is `(optimal - shuffled) / (optimal - loop)` is at least 0.8,
at caps 3 and 5. Simulation found almost all of the gain kept (OQ-13). If it
is below 0.5 at both caps, the fragility signal matters on real data and OQ-13
is reopened.

**H5, calibration.** The posterior's predicted switch probabilities
under-predict the observed switch rate on real reads: observed divided by
predicted, over spanned boundaries of the selected reads at cap 5, is above
1, with its bootstrap interval above 1. On simulation this happened only when
errors were unreported, so a value near 1 would say the reported qualities
are faithful.

**H6, margins as a signal.** The area under the ROC curve of the posterior
margin for predicting a switch at a spanned boundary, at cap 5, is above
0.6. Simulation's spanned-only value at the density of chromosome 20 is in
the corrected heterozygosity table of the development notebook.

H1 and H6 are descriptive. H2 and H3 test the conditional hypothesis of
`PROPOSAL.md`; H4 tests the OQ-13 finding.

## What this plan does not settle

Pedigree mode, other chromosomes, and other samples. A result on chromosome
20 of one sample is reported as that.
