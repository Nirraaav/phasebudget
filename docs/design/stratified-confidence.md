# Per-stratum error control (track A)

Track A of the extension plan (development repository). The registered confidence test
calibrated one threshold on HG002 variants inside the GIAB v4.2.1 benchmark
regions. On HG001 it failed, and after scoring the failure was confined to
variants outside those regions. This document states the change precisely
enough to be wrong.

## Objects

- A **pair** is two consecutive phased variants of one block, with a score
  `s` (the frozen score, `phasebudget.confidence.score`) and a loss `w`,
  which is 1 when the phaser's relation disagrees with the truth.
- A **group** `g` is a label fixed for each pair before any threshold is
  chosen. Track A uses the region stratum of the pair (below); track D adds
  the variant type.
- A **unit** is `(chromosome, 100 kb bin, group)`. A unit is *active* at
  threshold `t` when it holds a pair with `s >= t`, and *wrong* when it holds
  such a pair with `w = 1`.

## Strata

The stratum of a pair is the first of these that holds for either of its two
variants, by GIAB v3.6 GRCh38 stratifications (`benchmarks/stratify.py`):

1. `segdup`: segmental duplications;
2. `lowmap`: low mappability, all;
3. `tandem`: tandem repeats and homopolymers, with 5 bp of slop;
4. `difficult`: outside `notinalldifficultregions`, and in none of the above;
5. `easy`: both variants inside `notinalldifficultregions`.

The order puts mapping problems before sequence-context problems, because a
read placed on the wrong copy makes every allele it reports suspect.

## Procedure

For each group `g` separately:
1. The grid is `coverage_grid` over the calibration scores of `g`, from half
   of its pairs trusted to all of them, in steps of 0.25 per cent.
2. The threshold `t_g` is `learn_then_test` over the units of `g` at `beta`
   and `delta`.
3. If no threshold is certified, `g` trusts nothing.

`phasebudget.confidence.learn_then_test_by_group` implements it. Its oracle
is `learn_then_test` applied to each group's rows alone, which it must
reproduce exactly.

## Guarantee

Assume that, within each group, the test units and the calibration units are
exchangeable. Then for each group with a certified threshold, with
probability at least `1 - delta` over the calibration data, a test unit of
that group is wrong with probability at most `beta`. This is the registered
guarantee, applied to each group's units; the groups are calibrated on
disjoint data, so each statement holds separately.

A 100 kb bin holds up to five units, one per stratum. By the union bound, the
probability that a bin holds any wrong trusted pair is at most `beta` times
the number of its strata with a certified threshold. Stating each group at
`beta / 5` gives a bin-level bound of `beta`.

## What it does not claim

- It does not make units of different genomes exchangeable. The assumption
  is per group, and a test genome whose units inside a group differ from the
  calibration genome's (more paralogous copies, another platform) can break
  it, as HG001 broke the single threshold.
- A group that certifies nothing is reported, not hidden. Its pairs are
  untrusted, and the share of pairs in such groups is part of every
  usefulness number.
- Variants that the truth does not call heterozygous are outside the
  population. Track C found that many HG001 test SNVs outside GIAB are not
  heterozygous by the Platinum Pedigree truth; a phased relation to such a
  site has no right answer, and the guarantee says nothing about it.

## Development rerun with track B's HG002 calibration (stated 8 October, before running)

`stratified_gate.py --cross --calibration cal_r10`: as `--cross`, but the
calibration rows come from `cal_r10` (HG002 ONT sup, v5.0q heterozygous SNVs
and indels inside the v5.0q BED), restricted to SNV pairs because the HG001
test rows are SNV pairs. The question is the extension plan's: does
calibrating on v5.0q's larger set of difficult-region SNVs change which
strata hold on HG001? Descriptive, no rule: the same tables as the `--cross`
run of 6 October, which used the registered calibration rows (v4.2.1 SNVs),
are reported for comparison.
