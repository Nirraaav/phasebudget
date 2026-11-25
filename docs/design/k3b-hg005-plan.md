# Pre-registered K3b: the benchmark gap on HG005 against an assembly truth

Written on 27 October 2026. The commit that adds this file is the
registration; its push time is the third-party record. Creating this file is
what lets `benchmarks/genome_common.py` open the platform `k3_hg005`. The
author asked for it on 27 October.

## Why

Q3 (HG002 HiFi, v5.0q) and K3 (HG001 ONT, Platinum Pedigree truthset) found
that pairs with a variant outside GIAB's benchmark BED hold most phasing
errors: 0.53 and 0.94 of whatshap's and LongPhase's wrong pairs on HG002,
0.799 and 0.887 on HG001. The directions entry for K3
(the October directions memo (development repository), round 7) asked for a genome with an
assembly truth independent of the phasers. HG005's HPRC release 2 dipcall
(`data/hprc/HG005/`, made for W1's gate on 14 and 15 October) is one.

## Material seen before writing

Everything in the development notebook up to 27 October. HG005 was read before: GIAB
11 kb HiFi (H5, 1 October), ONT (X6, 13 October), HPRC HiFi with GIAB inputs
(K1b, scored 26 October), and W1's gate compared the release 2 dipcall with
GIAB transmission (about 9e-5 disagreement per pair). No phasing of HG005 on
the dipcall variants exists, and no count of errors outside GIAB's HG005 BED
has been made.

## Design

- **Reads.** HG005 HPRC HiFi (HudsonAlpha), winnowmap to GRCh38, the BAM
  `k1_hg005` uses, drawn at 30x as `confidence_run.py` draws it.
- **Phaser input.** Heterozygous SNVs and indels of dipcall 0.3 on the HPRC
  release 2 HG005 assemblies, inside dipcall's BED, phase removed, as W1's
  test genomes (`confidence_run.hetsnv`).
- **Truth.** The same dipcall calls' phase (`genome_common.dipcall_truth`).
- **Phasers.** WhatsHap 2.8 and LongPhase 2.0.2, as in K1b.
- **Pairs.** Consecutive phased pairs whose two variants the truth phases,
  21 confirmatory autosomes. Inside: both variants inside
  `HG005_GRCh38_1_22_v4.2.1_benchmark.bed`. Outside: at least one is not.
- **Platform.** `k3_hg005` in `genome_common.py`; data in
  `data/ext_k3_hg005`.

## Hypothesis

- **K3b (primary).** For both whatshap and LongPhase, the share of wrong
  pairs that lie outside the GIAB BED is at least 0.3.
- **Descriptive.** Error per pair inside and outside and their ratio, per
  group (`indel`, then the GIAB v3.6 strata), and the K2 cause split on each
  side, as `benchmarks/k2_k3_error_causes.py` prints them for K3.

## Risks stated now

- The release 2 assemblies used HG005's HiFi reads, among other data. Where
  assembly and phaser make the same error, the pair counts as right, which
  lowers the measured error on both sides of the BED.
- Outside GIAB's BED the truth itself is less certain. An outside "error" may
  be a truth error. K3b measures where the phaser and an independent
  assembly disagree; it does not adjudicate which is right.

## Scorer and order

`benchmarks/k2_k3_error_causes.py` with `--platforms k3_hg005 --k3-platform
k3_hg005 --k3-bed data/giab_hg005/HG005_GRCh38_1_22_v4.2.1_benchmark.bed`,
once, after `confidence_run.py` has finished all 21 chromosomes. Output
`benchmarks/results/runs/<date>_k3b_hg005.txt`. Deviations are dated
amendments below.

## Amendments
