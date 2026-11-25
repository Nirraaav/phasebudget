# Q4: certified cis or trans calls in clinical genes

Written 17 October 2026, before any Q4 number was computed. Development
gate, not registered. Direction Q4 of the October directions memo (development repository)
(round 6), which revives Y2 with the per-group certificates of the second
and third tests.

## Material seen before writing

Everything in the development notebook up to 17 October. Relevant to Q4:

- Track D (5 October, `docs/design/chain-relations.md`): the chain relation
  and its union bound, measured with pooled thresholds on HG002 ONT hac and
  HG001 ONT sup at 30x, all distances, all protein-coding genes.
- Y2's gate (7 October, `chain_relations.py --clinical --stratified`):
  3,893 ClinVar symbols with at least 10 P/LP alleles; per-stratum
  thresholds from HG002 ONT sup trusted 0.786 to 0.789 of HG002 hac
  relations and 0.690 to 0.702 of HG001's at `beta = 1e-3`, all distances.
  Genes that withheld almost everything included GBA1, CYP11B1, CYP11B2,
  STRC, SDHA, CHRFAM7A, CLCNKB, KRT6A, CATSPER2 and FLG.
- Y2 step 2 (CMRG genes do not withhold more) and step 3 (genes with an
  unprocessed pseudogene withhold more, descriptive), and S6 of the second
  test (partly supported, 13 of 14).
- The third test (X8): per-group thresholds from `cal_v5` held on
  `hg001_v5` in `easy` and in 88 of 90 difficult cells; whatshap `lowmap` at
  `beta = 3e-3` violated at 30x and 10x.

None of these used a per-call target rate, the `hprc_HG002` HiFi run, a
100 kb window, or the call-level validity test below.

## What a call is

The objects are track D's. Within one chromosome the phaser input variants
(the positions the stored fragments keep) are ordered by position. A
**link** joins two consecutive phased variants of one block, as
`phasebudget.confidence.phased_pairs` defines them; its score is the frozen
score, its **group** is `indel` when either end is an indel and otherwise the
stratum of `stratified_gate.strata`, and its **cell** is `(group,
chromosome, first position // 100 kb)`, the unit of the per-group
calibration.

For variants `i < j`, both phased, the phaser's **cis or trans call** is the
parity of the chain of links from `i` to `j`: cis when `a_i = a_j`. The chain
exists when every consecutive pair of phased variants from `i` to `j` is a
link; otherwise the pair is a **block break**. The call is **wrong** when
`a_i XOR h_i != a_j XOR h_j`, with `h` the truth's first-haplotype allele.

## Certification at a per-call rate

Per group `g` and `beta` in the grid `B = {1e-4, 3e-4, 1e-3, 3e-3, 1e-2}`,
`learn_then_test_by_group` gives a threshold `t_g(beta)` or none. For one
group the thresholds are nested: the grid does not depend on `beta`, and the
binomial p-value `P(X <= k | n, beta)` falls as `beta` grows, so the
fixed-sequence scan goes at least as far at a larger `beta`. A link trusted at
`beta` is trusted at every larger `beta`. The code checks this.

For a target rate `r` and a chain whose links lie in `m` distinct cells, let
`b(m, r)` be the largest `beta` in `B` with `m * beta <= r`. The call is
**certified at `r`** when every link of the chain is trusted at its group's
`t_g(b(m, r))`. By nesting, this is the same as "certified at some `beta` in
`B` with `m * beta <= r`".

**Claim.** If the call is certified at `r`, the truth phases every variant of
the chain, the test cells are exchangeable with the calibration cells of
their group, and every group the chain crosses was calibrated correctly
(probability at least `1 - s * delta` for `s` groups), then
`P(call certified and wrong) <= m * b(m, r) <= r`.

Proof. A wrong call has an odd number of wrong links, so at least one
(track D, part 1). Every link of a certified call is trusted at
`b = b(m, r)`, a quantity fixed by the positions and groups alone, not by any
score. So the event lies in the union over the call's `m` cells of "the cell
holds a wrong link trusted at `t_g(b)`", each of probability at most `b`.

The read of `r` used below is the one Learn then Test itself uses for units:
among certified calls, the share that is wrong. The same limits as track D
apply: chains through variants without truth phase, and the selection of
cells by gene, are not covered by the claim. Calls are reported with and
without such variants.

## Population and settings

- **Gene list, fixed now.** Protein-coding genes of GENCODE v50 basic (md5
  11e77cf1d4b78b43f57d608e3f7f2202) whose symbol has at least 10
  pathogenic or likely pathogenic alleles in ClinVar's
  `gene_specific_summary.txt` of 4 October 2026 (md5
  9761c2d9b82a96d8601b6129d607b312): `chain_relations.clinical_symbols()`,
  Y2's list, which the second test's S6 also used. A threshold of 10 rather
  than 1 keeps Y2's list (11,256 symbols have at least one P/LP allele, most
  of the protein-coding genome). No list of recessive genes or of the ACMG
  secondary findings is on the server, so no such subset is used.
- **Paralog lists, fixed now, descriptive.** (i) `cmrg_withholding.duplicated()`:
  protein-coding genes with an unprocessed pseudogene named after them in
  GENCODE v50 basic. (ii) A named list of genes known for paralog problems
  (`recalled`, not checked against a source): PMS2, SMN1, SMN2, GBA1,
  CYP21A2, HBA1, HBA2, STRC, TNXB, NEB, CYP11B1, CYP11B2, NCF1, PKD1, CYP2D6,
  SDHA, CHEK2, OTOA, IKBKG, HYDIN. Genes on chromosome 20 or X are not in the
  21 confirmatory autosomes and drop out.
- **Pairs.** Every unordered pair of phaser input variants that the truth
  phases, both inside one listed gene (`gene` record span), at most 100 kb
  apart, on the 21 confirmatory autosomes. A pair inside two overlapping
  genes counts once in pair totals and once in each gene's row.
- **Development genomes.** `cal_v5` (HG002 ONT sup v5.0.0), `cal_r10`
  (HG002 ONT sup 2023), `hprc_HG002` (HG002 HPRC HiFi, winnowmap), all
  against v5.0q; `hg001_v5` (HG001 ONT sup v5.0.0) against the Platinum
  Pedigree truth. Depths 30, 20, 10; phasers whatshap and LongPhase. All
  were seen before: this is development.
- **Calibration.** Within each genome, two folds by chromosome parity:
  thresholds from the odd chromosomes certify calls on the even ones and the
  reverse, and the two folds are pooled into one cell. A pair's thresholds
  never saw its chromosome. Separately, the transfer setting of the third
  test: `hg001_v5` with thresholds from `cal_v5` on all 21 chromosomes.
  Groups, score, units, grid start and step, and `delta = 0.05` are the
  third test's.
- **Targets.** `r` in {1e-3, 1e-2}.
- **Cells.** (setting, genome, depth, phaser, `r`): 4 x 3 x 2 x 2 = 48
  cross-fitted cells and 3 x 2 x 2 = 12 transfer cells, 60 in all.

## Measures

Per cell: pairs; certified pairs and their share; wrong certified pairs and
their rate; the same for certified calls whose chain holds a variant the
truth does not phase; the mean `m`; abstentions by cause: an endpoint
unphased, block break, untrusted link (counted once for each group holding an
untrusted link at `b(m, r)`). Per gene: pairs, certified share. Genes with
every pair certified, and the median certified share over genes with at
least 10 pairs. The genes that abstain most, with causes; the named paralog
genes; listed duplicated genes against the rest (median share, one-sided
Mann-Whitney, descriptive).

## Gate

Fixed now, before any Q4 number:

- **(a) Validity.** Calls in one gene share links, so certified pairs are not
  independent. Per cell, each gene with at least one certified call
  contributes one: the certified pair with the smallest CRC32 of
  `"gene:pos_i:pos_j"`, a choice that does not look at phase or truth. Over
  the `n` genes, `k` of these calls are wrong. A cell violates when the exact
  `P(X >= k | n, r)` is below `0.05 / 60`. (a) holds when no cell violates.
  Rates over all certified pairs are reported beside it, with no test.
- **(b) Usefulness.** On `hprc_HG002`, 30x, `r = 1e-2`, cross-fitted, at
  least 0.80 of the pairs get a certified call for at least one phaser.

Q4 passes when (a) and (b) both hold. If it passes, a registration for the
44 W1 genomes is drafted (drafted for review, not committed) for the
author's review; it is not committed by the run that drafts it.

## Code

`benchmarks/q4_clinical_cis_trans.py`; `tests/test_q4_clinical_cis_trans.py`
checks `b(m, r)`, the chain parity, the union bound on small blocks by
enumeration, and the pair enumerator against a brute-force recount. Chains
are extracted per (genome, depth, phaser, chromosome) into `data/r6_q4/`.
Output `benchmarks/results/runs/2026-10-17_q4_clinical_cis_trans.txt`.
