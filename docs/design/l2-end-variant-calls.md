# L2: within-gene cis or trans calls checked at the two end variants

Written on 22 November 2026, before any computation described here. The
commit that adds this file states the rules; it is not a blind registration,
because K4b already read both test genomes (see "Material seen").

## Claim

A cis or trans call between variants `i` and `j` is the phaser's relation
between their two alleles. It is wrong exactly when that relation differs
from the truth's. The links in between matter only through that relation. Q4
and K4b certified a call only when every link of the chain from `i` to `j`
was trusted and agreed with the assembly. One link the assembly does not
relate, often an indel the assembly writes differently, then blocks the
call, and more links means more chances to block. In K4b, whatshap on
NA12877 certified 0.863 of pairs up to 1 kb apart but 0.553 of pairs 10 to
100 kb apart at `r = 1e-2` (LongPhase 0.857 and 0.599).

The end rule checks the thing reported: the call is trusted when the
person's assembly gives the same relation for the two end variants. The
error rate of agreeing calls is calibrated directly, per call, on the same
six-genome pool as K4b.

## Material seen

Everything in the development notebook up to 22 November 2026, including K4b's score
(`benchmarks/results/runs/2026-10-27_k4b_score.txt`): per test genome and
phaser, certified shares overall and by distance, wrong certified pairs,
wrong CRC picks, and per-gene shares of the named paralog genes. No end-rule
quantity (end agreement, its error, or read support at the two ends) has
been computed on any genome. L2 is a development analysis with rules stated
before computing, not a blind test.

## Definitions

Unchanged from K4b unless listed: the pool (`hprc_HG002` with v5.0q, and
`phifi_NA12879`, `NA12881`, `NA12882`, `NA12885`, `NA12886` with the
Pedigree truthset; 21 chromosomes; 30x), the test genomes (`phifi_NA12877`,
`k1_hg005`), both phasers, the assemblies and their values (K1b's tables in
`data/k1/tables/`, Q6 contig rule), the gene list
(`chain_relations.clinical_symbols()`), the pair population (truth-phased
input variants `i < j` inside one listed gene, at most 100 kb apart), the
targets `r` in {1e-3, 1e-2}, the distance bins (up to 1 kb, 1 to 10 kb, 10 to
100 kb) and the CRC32 pick of `"gene:pos_i:pos_j"`.

- **Phaser call.** `q6.relation` of the phaser's `[allele, block]` values at
  `i` and `j`: defined when both are phased in one block.
- **End agreement (primary arm).** The assembly's relation at `i` and `j`
  (`q6.relation` on its `[allele, block]` values) exists and equals the
  phaser call (`k1_assembly_agreement.agrees`).
- **Read agreement (descriptive arm).** Among the 30x realigned fragments
  that call an allele at both `i` and `j`, at least 2 exist and a strict
  majority give the phaser's relation. This arm needs no assembly. Read
  support is counted only for variants inside a listed gene.
- **Group.** (distance bin, class); class is `hard` if either end is an
  indel or the pair's stratum (`stratified_gate.strata`) is not `easy`, and
  `easy` otherwise. Six groups.
- **Wrong.** Phaser call differs from the truth's relation.

## Calibration

For each phaser, arm, group and `r`: over the pool genomes, every agreeing
pair in the group is a candidate, and each (genome, gene) contributes one
pick, the candidate with the smallest CRC32. With `n` picks and `k` wrong,
the group is certified at `r` when `P(X <= k | n, r) <= 0.05` (exact
binomial). With `k = 0` that needs `n >= 299` at `r = 1e-2` and `n >= 2,995`
at `r = 1e-3`.

At test time a pair is certified at `r` when it agrees under the arm and its
group is certified at `r`. Nothing else is trusted.

The guarantee is about the same random draw the validity test makes: a
gene's pick among certified calls. It assumes the test genome's agreeing
calls in a group are no more error-prone than the pool's, the assumption
K1b and K4b already make.

## Decision rules (primary arm)

Cells: (test genome, phaser, `r`), 2 x 2 x 2 = 8, level 0.05/8 = 6.25e-3,
as in K4b.

- **L2-1, validity.** In each cell each gene with a certified call
  contributes the certified pair with the smallest CRC32; `k` of `n` picks
  are wrong. No cell has `P(X >= k | n, r)` below 6.25e-3.
- **L2-2, long-range gain on NA12877.** At `r = 1e-2`, at least one phaser
  certifies at least 0.70 of NA12877's within-gene pairs 10 to 100 kb apart
  (K4b: whatshap 0.553, LongPhase 0.599).
- **L2-3 (reported, not decided).** HG005's shares, overall and by distance
  (K4b at 10 to 100 kb: 0.862 and 0.870).

L2 succeeds when L2-1 and L2-2 both hold.

## Descriptive outputs

Per cell: overall certified share against K4b's; shares by distance; wrong
certified pairs over all certified pairs; the union with K4b's chain rule is
not computed (its validity would need a separate bound). The read arm in the
same cells, with the same validity test printed but not decided. Pool picks
and wrong picks per group. The named paralog genes of Q4 (PMS2, OTOA,
CYP11B1, CYP11B2, PKD1 lost calls under K4b).

## Risk stated in advance

hifiasm switches phase inside contigs (Q6: 8,549 runs of agreement over its
blocks on HG002, median 47 variants). At far-apart ends a switch makes the
assembly disagree with a correct phaser, which abstains: safe, but it lowers
the gain. If a phaser switch and an assembly switch fall between the same
two ends, the pair agrees and is wrong; calibration must catch that rate.

## Order

1. `benchmarks/l2_end_calls.py` with a test of the end rule against a
   brute-force parity check.
2. `--dry-run`: HG002's odd chromosomes split into two pseudo-genomes plus
   the five Pedigree genomes as pool, HG002's even chromosomes as the test.
   Output `benchmarks/results/runs/2026-11-22_l2_dry_run.txt`.
3. The run on NA12877 and HG005, once, unchanged after the dry run unless an
   amendment below says why. Output
   `benchmarks/results/runs/2026-11-22_l2_end_calls.txt`.

## Amendments

### Amendment 1, 22 November 2026: the dry run

Written before NA12877 or HG005 was read by L2. `benchmarks/l2_end_calls.py`
(63fed98) ran with `--dry-run` on the server
(`benchmarks/results/runs/2026-11-22_l2_dry_run.txt`, 780 s). No cell
violates at 6.25e-3 (at most 1 wrong pick of 1,147 genes). At `r = 1e-2`
the assembly arm certifies 0.946 (whatshap) and 0.923 (LongPhase) of HG002's
even-chromosome pairs, and 0.940 and 0.917 of pairs 10 to 100 kb apart. The
risk stated above shows in the pool: agreeing pairs 10 to 100 kb apart are
wrong about 3 times in 10,000 (whatshap 18,909 of 62.3 million in `easy`),
against none up to 10 kb in `easy`, so that bin certifies at 1e-2 and not
at 1e-3. The read arm certifies few pairs beyond 10 kb (about 0.10), as
expected from read length. The run on the test genomes follows, unchanged.
