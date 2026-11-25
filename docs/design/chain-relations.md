# Relations between any two variants of a block (track D)

Track D of the extension plan (development repository). The registered confidence test
scored only consecutive pairs. A compound heterozygous call joins two
variants that can be far apart in one gene. This document states the relation,
when it is trusted, and what the per-unit guarantee says about it.

## Objects

- A block holds phased variants `v_1 < ... < v_n`. Link `k` joins `v_k` and
  `v_{k+1}`; its score is the frozen score of the consecutive pair, and it is
  **trusted** at threshold `t` when its score is at least `t`.
- The **relation** of `v_i` and `v_j`, `i < j`, is cis or trans as the phaser
  calls it. It is **trusted** when every link `i, ..., j - 1` is trusted.
- Under the complementary diploid model (contract item 3), write
  `d_k = a_k XOR h_k` for the phaser's first-haplotype allele `a_k` and the
  truth's `h_k`. The relation is **wrong** when `d_i != d_j`. It is defined
  only when the truth phases both endpoints.
- A **unit** is as in the registered test: `(chromosome, 100 kb bin)` of a
  link's first variant, or with a stratum added in track A.

`phasebudget.confidence.chain_counts` counts relations, wrong relations,
trusted relations and wrong trusted relations of one block in linear time.
Its oracle is the enumeration of every pair in `tests/test_confidence.py`.

## Proposition

Let the relation of `v_i` and `v_j` be trusted, and let the truth phase every
variant `v_i, ..., v_j`. Let `U` be the set of units that hold links
`i, ..., j - 1`, with `m = |U|`.

1. **Deterministic part.** If no link of the chain is wrong, the relation is
   right. Proof: `d_j - d_i` is the sum, modulo 2, of `d_{k+1} - d_k` over the
   links, and a link is wrong exactly when its term is 1. The converse fails:
   two wrong links cancel, so the bound below is conservative.
2. **Probabilistic part.** The event "relation trusted and wrong" lies inside
   the union, over `u` in `U`, of "`u` holds a wrong trusted link". Assume
   the test units are exchangeable given the calibration data. Then each unit
   has the same probability of holding a wrong trusted link, the risk `R(t)`
   that Learn then Test bounds. On the event, of probability at least
   `1 - delta` over the calibration data, that `R(t) <= beta`, the union
   bound gives

   `P(relation trusted and wrong) <= m * beta`.

   With track A's per-stratum thresholds, the chain's units can belong to
   several strata. The bound then holds on the event that every stratum the
   chain crosses is certified correctly, of probability at least
   `1 - (number of those strata) * delta`.

## What it does not claim

- **Variants the truth does not phase.** The calibration only sees links
  whose two variants the truth phases. A chain through a variant without
  truth phase, which includes every false heterozygous call, has links the
  guarantee never covered. The proposition says nothing about such chains,
  and the development numbers report them separately.
- **Selection by gene.** A relation is chosen because both variants lie in
  one gene. The units of a gene are not a random draw of units; if genic
  units hold wrong links more or less often than the average unit, the
  per-relation bound does not follow from the per-unit guarantee. The
  development numbers measure it; they cannot prove it.
- **No bound on a whole gene.** The statement is per relation. A gene with
  many trusted relations can hold a wrong one with probability well above
  `m * beta`.

## Development measure

On HG002 ONT hac (truth v5.0q) and HG001 ONT sup (truth Platinum 2017), 30x,
`whatshap phase` and LongPhase, with the registered pooled thresholds
calibrated on HG002 ONT sup:

- the population is every pair of truth-phased heterozygous SNVs in one
  phaser block and one protein-coding gene of GENCODE v50 basic (`gene`
  records with `gene_type "protein_coding"`; nested or overlapping genes each
  count their own pairs);
- per beta: relations, the share trusted, wrong trusted relations, their
  rate, and the mean bound `beta * m` over trusted relations, with `m`
  counted as the bins from the first to the last endpoint, an upper bound on
  the units the chain holds;
- the same counts restricted to chains whose every variant the truth phases,
  the population of the proposition.

Development data were seen before. The registered test of this relation is
track B's, after `docs/design/stratified-run-plan.md` is committed.
