# L7: NA12878 HiFi, five predictions registered before scoring

Registered 22 November 2026, before any NA12878 HiFi phasing exists. The
`phifi_NA12878` run was started once for R1 part (c) on 16 October and
stopped after variant extraction on a few chromosomes; no phase was computed
and nothing was compared with a truth. NA12878's ONT data were scored in X8,
K3 and R1; her HiFi reads were not.

## Data

- Reads: Platinum Pedigree HiFi, `NA12878.GRCh38.haplotagged.bam`
  (platform `phifi_NA12878`, `benchmarks/genome_common.py`), depths 30x, 20x
  and 10x, whatshap and LongPhase, 21 confirmatory chromosomes, the same
  pipeline as W1b and K1b (`benchmarks/confidence_run.py`).
- Truth: the Platinum Pedigree truthset for NA12878. Its phase comes from
  transmission in the pedigree, not from an assembly.
- Assembly: single-sample hifiasm 0.25.0 and dipcall 0.3, as in K1b, with one
  stated deviation: the reads are a random subsample to about 30x (the K1b
  pool used every read). NA12877's full-depth assembly took 67 hours, which
  does not fit the schedule. The subsample fraction is fixed from the BAM's
  read count and mean read length before hifiasm runs, seed 78.
- Parents for the family key: NA12891 and NA12892 genotypes from the 1000
  Genomes high-coverage panel (`data/kg_panel`).

## Predictions and decision rules

Cells are (depth, phaser, beta, group) with betas 3e-4, 1e-3, 3e-3, delta
0.05, groups `easy`, `difficult`, `segdup`, `lowmap`, `tandem`, `indel`, 100
kb units, as in W1b and K1b.

- **P1, W1b's cutoffs break on all links.** W1b's pooled thresholds (HG002
  and the 20 consensus-truth HPRC pool genomes, `benchmarks/w1b_score.py`),
  tested on NA12878 against the Pedigree truth: at least one cell in
  `segdup`, `lowmap` or `indel` violates at 0.05/108. This is L1's
  prediction; NA12878 shares no family with the W1b pool.
- **P2, they hold on links the family key labels.** Same thresholds, test
  links restricted to those where NA12878's assembly relates the pair (Q6
  contig rule) with the relation transmission from NA12891 and NA12892 gives
  (L1 amendment 1's key): 0 violating cells at 0.05/108. Random control as
  in L1: links kept at random with the same share per cell; reported.
- **P3, K1b replicates.** K1b's pooled agreeing-pair thresholds
  (`benchmarks/k1b_score.py` pool), agreement with NA12878's assembly, 30x,
  both phasers: 0 violating cells at 0.05/36; `segdup` at 1e-3 certifies
  with trusted share at least 0.5 for one phaser.
- **P4, L2 replicates.** L2's end-variant rule and pool
  (`benchmarks/l2_end_calls.py`), assembly arm, r = 1e-2: L2-1 validity
  holds, and at least 0.7 of NA12878's within-gene pairs 10 to 100 kb apart
  are certified for one phaser.
- **P5, the benchmark gap with an assembly-free truth.** At 30x, of wrong
  pairs against the Pedigree truth, the share with a variant outside the GIAB
  HG001 v4.2.1 benchmark BED is at least 0.5 for both phasers, and the wrong
  rate per pair outside is at least 3 times inside.

P1 and P5 need no assembly. If the assembly is not done by 25 November, P1
and P5 are scored alone and P2 to P4 later, each still before its own first
comparison with the truth.

## Known limits

NA12878 is the mother of the five Pedigree children in the K1b and L2 pools,
so P3 and P4 are not independent of their calibration: sibling and parent
wrong units share loci (L1 step C). P1 and P2 use the W1b pool, which has no
member of this family. The 30x assembly may switch phase more often than the
pool's full-depth assemblies, which makes P2 to P4 harder, not easier.

## Amendment 1 (22 November 2026, before any NA12878 result)

aria2c failed to finish the NA12878 BAM (TLS handshake failures at 199 of
212 GB). The assembly reads are drawn instead through htslib with
`scripts/l7_reads.py`: 10 Mb chunks, a read kept when crc32 of
"78:name" falls below the fraction for about 30x, primary alignments only.
Unmapped reads are not drawn, which K1b's assemblies included. No result
of any kind has been computed for NA12878.

## Amendment 2: scorer and dry run (22 November 2026, before any NA12878 result)

The scorer is `benchmarks/l7_score.py`, with `tests/test_l7_score.py` for P5's
inside and outside split. It changes no rule above. It reuses W1b's pool and
`learn_then_test_by_group` (P1, P2), L1's `rows_with_pos`, `restrict` and
`matched_random` with a family key built from the test genome's K1b table and
its own parents' panel genotypes (P2), K1b's `pair_rows` and pool (P3), and
L2's extraction, tallies and pool (P4). Three details the registration left
open:

- P2's random control keeps links at random with the key's share per group,
  seed 0, as in L1. It is reported and decides nothing.
- P4's level is 0.05 over this test's cells: 2 phasers times 2 targets, so
  0.05/4, built the way L2 built 0.05/8 for two test genomes.
- P5 counts a pair inside when both variants lie in the BED, as K3 did.

`--no-assembly` scores P1 and P5. `--dry-run` puts NA12879 (parents NA12877
and NA12878) in NA12878's place and reads no NA12878 phase. NA12879 is in the
K1b and L2 pools, so the dry run's P3 and P4 numbers check the code only.

## Result (24 November 2026, `2026-11-24_l7_score.txt`, 9,325 s)

All five predictions are supported.

- **P1 holds.** On all links, 31 of 108 cells violate at 0.05/108: 13
  `segdup`, 11 `lowmap`, 6 `indel`, 1 `tandem`.
- **P2 holds.** On links the family key labels, 0 of 108 cells violate. The
  random control (key's share kept at random per group) has 20 violating
  cells, so the key's 0 is not lost power.
- **P3 holds.** With K1b's pool thresholds on agreeing pairs at 30x, no cell
  violates at 0.05/36 (smallest p 0.90). `segdup` at 1e-3 certifies with
  trusted share 0.779 (whatshap) and 0.779 (LongPhase). At 3e-4 `segdup`
  does not certify for either phaser.
- **P4 holds.** At r = 1e-2 no cell violates at 0.05/4 (1 and 2 wrong picks
  of about 2,460), and the certified share of within-gene pairs 10 to 100 kb
  apart is 0.904 (whatshap) and 0.893 (LongPhase). At r = 1e-3 that share is
  0 for both.
- **P5 holds.** At 30x the share of wrong pairs with a variant outside the
  GIAB HG001 v4.2.1 BED is 0.914 (whatshap) and 0.848 (LongPhase); the wrong
  rate per pair outside over inside is 40.9 and 23.0.

The known limits above still apply: NA12878 is the mother of the five
Pedigree children in the K1b and L2 pools, so P3 and P4 are not independent
of their calibration. P1 and P2 use the W1b pool, which has no member of
this family.
