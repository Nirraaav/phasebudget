# R3: certificates for the pairs family labels can reach

Written 16 October 2026, before any run. Development gate, not registered.
Direction R3 of the October directions memo (development repository); it follows W5 part 3
(`docs/design/w5-label-noise.md`, the development notebook 15 October).

## Statement

A child's pair is **eligible** under a label filter when the filter labels
both of its sites from the parents' calls. Eligibility is known at test
time, before any truth, from the same calls. Learn then Test calibrated on
the eligible pairs of some units controls the wrong-unit share of the
eligible pairs of exchangeable units, counted on eligible pairs only: the
calibration and the test then see the same population. It says nothing
about ineligible pairs. A site is ineligible with probability one when both
parents are heterozygous there, so no reweighting of eligible pairs can
stand in for those pairs; for the rest (the replicates disagree, or a call
is missing) reweighting needs an assumption that the reason a label is
missing is unrelated to phase error, which W5 part 3 found false in
`lowmap`. The procedure therefore abstains on every ineligible pair.

## Measurement

HG002 `cal_v5`, labels from Clair3 calls on both flowcell replicates of
HG003 and HG004 (`benchmarks/w5_noise.py`), filters `agree` (both
replicates give the same label), `all` (replicate A) and `double` (both
parents homozygous, opposite). Odd chromosomes calibrate on eligible pairs
with label-based errors, as in W5; even chromosomes are scored against
v5.0q twice: on eligible pairs (the statement) and on all pairs (W5's
scoring, which the statement does not cover). Depths 30, 20, 10; whatshap
and LongPhase; `beta` 3e-4, 1e-3, 3e-3; groups of `stratified_gate.strata`
plus `indel`.

Ineligible sites are split by reason under `agree`: `uninformative` (in
both replicates both parents are heterozygous or have a non-PASS record, so
these calls cannot label the site) and `other` (the replicates differ). Per group it reports the share of
pairs in each class and whatshap's and LongPhase's error there.

## Gate

Fixed before looking, for the `agree` filter:

- **(a) valid where it claims:** no violating cell when scored on eligible
  pairs, at the Bonferroni level over those cells;
- **(b) worth having:** `easy` at `beta = 1e-3` certifies in all six depth
  and phaser cells, and the trusted share of all even-chromosome `easy`
  pairs (abstaining on ineligible ones) is at least 0.70.

R3 passes when both hold.

## Code

`benchmarks/r3_label_eligible.py`; `tests/test_r3_label_eligible.py` checks
the eligibility and reason rules on fixed sites.
