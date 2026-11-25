# Pre-registered tool-level test: WhatsHap with exact read selection

Written on 30 September 2026, before any HG002 HiFi read of a confirmatory
chromosome was read for this test. The commit that adds this file fixes
that.

Material seen before writing:
- the HiFi and ONT library results, and the uncalled-variant correction;
- the mock review of 30 September;
- a smoke test on the chromosome-20 pilot BAM. It is not part of any test.
  At `--internal-downsampling 3`, stock WhatsHap had 4 switches in 43,368
  pairs and WhatsHap with exact selection 2 in 43,102; HiPhase and LongPhase
  had 0.

## Why

The library pipelines phase with this project's DP, not with `whatshap
phase`. A reader will ask whether exact selection changes the tool
people use. `benchmarks/whatshap_exact.py` runs WhatsHap 2.8 unchanged,
except that `readselection` is replaced by the maximum-quality selection
under the same cap. The weights are WhatsHap's own allele qualities; reads
calling fewer than two variants are never selected, as in WhatsHap.

## Data and tools

- **Reads.** The HG002 PacBio Sequel II 11 kb GRCh38 alignment of the
  genome-wide run, streamed one chromosome at a time, full depth only. Each
  slice is deleted after its tools finish.
- **Chromosomes.** The 21 confirmatory autosomes.
- **Variants.** The genome-wide run's `hetsnv.vcf.gz`.
- **Runs.**
  - `whatshap phase` 2.8 (`stock`) at `--internal-downsampling` 3, 5 and
    15;
  - `benchmarks/whatshap_exact.py phase` (`exact`) at 3, 5 and 15;
  - HiPhase 1.7.0 with `--ignore-read-groups`;
  - LongPhase 2.0.2 with `--pb`.
- **Truth.** GIAB v5.0q, assigned as in genome-run amendment 2.

## Hypotheses

Statistics are those of `genome-run-plan.md`, amendment 1: common pairs,
the exact sign test, and the 1 Mb window bootstrap.
- **W1.** At cap 3, `exact` has fewer discordant errors than `stock`,
  significant at 0.05/3.
- **W2.** At caps 5 and 15, the same comparison, at 0.05/3 each. These are
  reported, and no direction is predicted.
- **W3 (descriptive).** Switch error, phased fraction, block N50 and
  runtime of every run, including the uncapped HiPhase and LongPhase.

## Safety

One slice on disk at a time, a download only with at least 20 GB free, and
tools at `nice 19`, each under a memory ceiling of 6 GB.
