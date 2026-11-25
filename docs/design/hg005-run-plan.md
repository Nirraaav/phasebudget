# Pre-registered replication on a second individual: HG005 HiFi

Written on 30 September 2026, before any HG005, HG006 or HG007 file was
downloaded or read. The commit that adds this file fixes that. Material
seen before writing:
- the HG002 HiFi and ONT designs and the HiFi results, including the
  uncalled-variant correction (`genome-run-plan.md`, amendment 3);
- the GIAB file listings for the Chinese trio, which contain no calls or
  reads.

## Why

Every HG002 result is one individual. HG005 is a Han Chinese son with the
same kind of GIAB benchmark and HiFi reads from the same chemistry, so it
tests whether the HG002 HiFi findings are a property of the method or of one
genome.

## Data

- **Reads.** `.../ChineseTrio/HG005_NA24631_son/PacBio_SequelII_CCS_11kb/HG005_GRCh38/HG005_GRCh38.haplotag.10x.bam`
  (about 30-fold; README of 16 May 2019), streamed one chromosome at a time.
  Haplotag fields are never read.
- **Chromosomes.** Autosomes 1 to 22 except 20, all confirmatory.
  Chromosome 20 is the smoke test and the held-out chromosome for the
  vote model's allele error, and enters no test.
- **Depths.** The full slice and a seeded 10-fold subsample.
- **Variants.** HG005's biallelic heterozygous SNVs inside its v4.2.1
  benchmark BED.
- **Fragments.** WhatsHap 2.8's reader with realignment, mapping quality
  20, log-odds, and read groups ignored if the BAM has none.

## Truth

GIAB has no v5.0q-style assembly benchmark for HG005. The primary truth is
transmission from the parents' v4.2.1 calls (HG006, HG007), built by
`build_trio_truth.py` exactly as for HG002, without the hifiasm filter. On
HG002 that filter removed 2,354 correct sites and one error. The hifiasm
phase transfer (`HG005_GRCh38_1_22_v4.2.1_highconf_hifiasm_v11_phasetransfer.vcf.gz`)
is reported as a secondary truth.

## Pipelines and caps

`whatshap`, `optimal`, `two_level`, `merge_optimal_1e3` and `connected`, at
caps 3, 5, 10 and 15, phased with `CoverageDP`, with uncalled variants
unphased. `whatshap phase` 2.8 with defaults as the external tool.

## Hypotheses

Statistics are those of `genome-run-plan.md`, amendment 1.
- **H5-1, deployed gap.** WhatsHap's selection reaches less than 0.95 of the
  certified optimum at full depth and cap 15, pooled. Supported if below
  0.95.
- **H5-2, the HiFi regime replicates.** `whatshap` against `optimal`, sign
  tests at caps 3 and 5, full and 10-fold, level 0.05/4.
  - Supported if both cap-3 tests are significant in `optimal`'s favour and
    neither cap-5 test is significant.
  - Refuted if neither cap-3 test is significant.
- **H5-3, the vote model.** Its prediction, with `e` measured on HG005
  chromosome 20, is computed and committed before scoring. It is judged by
  rules R1 and R2 of `ont-prediction.md`.
- **Descriptive.** Block N50 of `two_level` and `connected` against
  `optimal` at cap 3.

## Safety

The ONT run's scheduler, limits and clean-up. This run starts only after
the ONT run has finished, so that no two runners share the disk.
