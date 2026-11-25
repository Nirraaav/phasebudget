# Pre-registered tool-level test on ONT: WhatsHap with exact read selection

Written on 1 October 2026, before any ONT read was phased by `whatshap
phase` with exact selection, and before the genome-wide tie-aware re-score
(`benchmarks/tie_rescore.py`) finished. The commit that adds this file
fixes that.

Material seen before writing:
- every registered result, including the HiFi tool-level test
  (`tool-run-plan.md`). There, W1 was not supported: +0.11 [-0.04, +0.25]
  per 10,000 at cap 3.
- the ONT library results (O1 to O7). In the library, WhatsHap's selection
  minus the exact one is +8.81 at cap 3 and +0.83 at cap 5, full depth, per
  10,000.
- the chromosome-20 finding that the library phases tied sites, which
  `whatshap phase` leaves unphased, together with the chromosome-20
  tie-aware counts (`2026-10-01_ties_ont_chr20_10x.txt`, and the development notebook, 1
  October, night).

## Why

The fourth mock review asks whether the ONT effects are properties of the
library or of the tool people use. `whatshap phase` abstains at tied sites,
so the library's effect may shrink inside the tool.

## Data and runs

- **Reads.** The HG002 ONT R10.4.1 sup alignment of the ONT run
  (`genome_common.BAM_URL` for `PHASEBUDGET_PLATFORM=ont`), streamed one
  chromosome at a time, full depth only. Each slice is deleted after its
  runs finish.
- **Chromosomes.** The 21 confirmatory autosomes.
- **Variants.** The ONT run's `hetsnv.vcf.gz`.
- **Runs.** All runs use `--ignore-read-groups`, as in the ONT run.
  - `whatshap phase` 2.8 (`stock`) at `--internal-downsampling` 3, 5 and
    15;
  - `benchmarks/whatshap_exact.py phase` (`exact`) at 3, 5 and 15.
- **Truth.** GIAB v5.0q.

## Hypotheses

The statistics are those of `tool-run-plan.md`: common pairs, the exact
sign test, and the 1 Mb window bootstrap.
- **WO1.** At cap 3, `exact` has fewer discordant errors than `stock`,
  significant at 0.05/3.
- **WO2.** At cap 5, the same comparison at 0.05/3.
- **WO3.** At cap 15, the same comparison at 0.05/3. It is reported, and no
  direction is predicted.
- **WO4 (descriptive).** For every run: switch error, phased fraction,
  block N50 and runtime.

## Safety

One slice on disk at a time. A download starts only with at least 45 GB
free, since the chromosome-1 slice is about 37 GB. Tools run at `nice 19`,
each under a memory ceiling of 6 GB.
