# Pre-registered K4b: within-gene cis or trans calls with assembly agreement

Written on 27 October 2026. The commit that adds this file is the
registration; its push time is the third-party record. The author asked for
it on 27 October, after K1b was scored.

## Why

K4 (`docs/design/k1-assembly-agreement.md`, development, HG002 HiFi, two
folds by chromosome parity) certified a within-gene cis or trans call for
0.548 (whatshap) and 0.647 (LongPhase) of the pairs at a per-call rate of
1e-2 when every link of the chain had to agree with the assembly, against
0.038 and 0.252 without agreement, with 0 and 1 wrong picks over about 2,600
genes. K4b asks whether that holds when the thresholds come from K1b's
six-genome pool and the calls are made on K1b's two test genomes.

## Material seen before writing

Everything in the development notebook up to 27 October, including K1b's score
(`benchmarks/results/runs/2026-10-26_k1b_score.txt`). K1b read both test
genomes: per group, their pair counts, agreement shares, wrong pairs among
agreeing and other pairs, and every certificate cell. No within-gene pair,
chain or call of NA12877 or HG005 has been computed. This is the main limit of
K4b: its test genomes are not unseen, only its question is.

## Design

Unchanged from K1b unless listed: the pool (HG002 `hprc_HG002` and the five
Platinum Pedigree members `phifi_<sample>`, 21 chromosomes, 30x, agreeing
pairs only), the agreement rule, the groups, the score, the 100 kb unit
(genome, chromosome, index), `learn_then_test_by_group` and `delta = 0.05`.
Unchanged from Q4 (`docs/design/q4-clinical-cis-trans.md`) unless listed: the
link, the chain, the call and its certificate, `b(m, r)`, the gene list
(`chain_relations.clinical_symbols()`, protein-coding GENCODE v50 genes with
at least 10 ClinVar P/LP alleles), the pair population (both variants
truth-phased, inside one listed gene, at most 100 kb apart, 21 autosomes) and
the CRC32 pick.

- **Grid.** Q4's `beta` grid {1e-4, 3e-4, 1e-3, 3e-3, 1e-2}, thresholds
  checked for nesting as in Q4.
- **Thresholds.** One set per phaser, from the pool's agreeing pairs.
- **Links on the test genomes.** A link the test genome's own assembly
  (K1b's hifiasm 0.25.0 assembly and dipcall table) does not agree with gets
  a score below every threshold (`k1_assembly_agreement.rewrite_links`).
- **Test genomes.** NA12877 (`phifi_NA12877`, Pedigree truthset) and HG005
  (`k1_hg005`, GIAB v4.2.1 trio transmission, SNVs only; a chain may pass
  through indels the truth does not phase, as in Q4).
- **Targets.** `r` in {1e-3, 1e-2}.
- **Cells.** (test genome, phaser, `r`): 2 x 2 x 2 = 8. Level 0.05/8 =
  6.25e-3.

## Hypotheses

- **K4b-1, validity (primary).** In each cell, each gene with a certified
  call contributes the certified pair with the smallest CRC32 of
  `"gene:pos_i:pos_j"`; `k` of these `n` picks are wrong. No cell has
  `P(X >= k | n, r)` below 6.25e-3.
- **K4b-2, usefulness on NA12877 (primary).** At `r = 1e-2` at least one
  phaser certifies at least 0.5 of NA12877's within-gene pairs.
- **K4b-3, usefulness on HG005 (primary).** The same on HG005.
- **K4b-4 (descriptive).** The same calls with thresholds from all pooled
  pairs and no agreement rule (the frozen-score baseline), in the same
  cells; wrong certified pairs over all certified pairs; certified share by
  distance (under 1 kb, 1 to 10 kb, 10 to 100 kb); the named paralog genes
  of Q4.

## Scorer and order

`benchmarks/k4b_score.py` reads K1b's tables (`data/k1/tables/`) and cached
rows for the pool, and Q4's chain files (`data/r6_q4/chains/`, extracted if
missing) for the test genomes. It is committed as a dated amendment with a
dry run that reads no test genome: HG002's odd chromosomes (split in two
pseudo-genomes) and the five Pedigree genomes form the pool, and HG002's
even chromosomes stand in for a test genome. Then the scorer runs once.
Every threshold is compared at full precision. Deviations are dated
amendments below.

## Software and machine

the KCDH compute server, the third test's Python environment, the K1b outputs.

## Amendments

### Amendment 1, 27 October 2026: the scorer and its dry run

Written before any test genome's chain was read. `benchmarks/k4b_score.py`
(9b999f8) ran with `--dry-run`
(`benchmarks/results/runs/2026-10-27_k4b_dry_run.txt`, 1,121 s): pool
HG002's odd chromosomes in two pseudo-genomes and the five Pedigree genomes,
test HG002's even chromosomes. No cell violates at 6.25e-3 (0 wrong picks of
1,145 genes in every agreement cell). At `r = 1e-2` the certified share is
0.607 (whatshap) and 0.712 (LongPhase) with agreement, against 0.039 and
0.286 for the frozen-score baseline. The registered scorer runs next,
unchanged.

### Amendment 2, 27 October 2026: the score

The registered scorer ran once at beb824b
(`benchmarks/results/runs/2026-10-27_k4b_score.txt`): K4b-1, K4b-2 and K4b-3
are supported. Seen after scoring: in some duplicated clinical genes
(PMS2, OTOA, CYP11B1, CYP11B2, PKD1 on NA12877) the agreement rule certifies
fewer pairs than the frozen score, because the assembly does not relate them.
The development notebook, 27 October.
