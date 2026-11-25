# R2: certifying a held-out sequencing batch from its own labelled units

Written 16 October 2026, before any run. Development gate, not registered.
Direction R2 of the October directions memo (development repository).

## The floor

Learn then Test certifies a group at `beta` when the share of 100 kb units
holding a wrong trusted pair is shown to be at most `beta` at level `delta`.
Suppose the new batch's labels are the only information about its error
rate given the score (no assumption ties the pool's error rate to the new
batch's). Compare two batches with the same reads and scores: in one the
risk at a threshold is `beta + epsilon`, in the other 0. With `n` labelled
units, the first gives all-correct labels with probability at least
`(1 - beta - epsilon)^n`. A procedure valid at `delta` that certifies the
second batch with probability `p` when all labels are correct certifies the
first with probability at least `p (1 - beta - epsilon)^n`, so

    n >= ln(p / delta) / -ln(1 - beta)   (as epsilon goes to 0)

which is about 2,995 units at `beta = 1e-3`, `delta = 0.05`, `p = 1`. The
pool's units cannot lower this, and neither can unlabelled units of the new
batch, which are the same in both. Prediction-powered inference corrects
the pool's estimate with the new batch's labels and meets the same bound.

## What the gate measures

Batch shift as W3 found it: HG002 (`cal_v5`, 2025 basecaller) is the new
batch, the other six genomes of `benchmarks/w3_pooled.py` (HG001 and five
Platinum Pedigree children, `hg001_v5`, `ped_*`) are the pool. Frozen score;
depths 30 and 20; whatshap and LongPhase; groups of `stratified_gate.strata`
plus `indel`; `beta` in 3e-4, 1e-3, 3e-3; `delta` 0.05.

From HG002's odd chromosomes, `k` units are drawn at random (`k` = 300,
1,000, 3,000, and all odd units), 20 seeds for each `k` below all. Methods,
each calibrated per group:

- **pool**: Learn then Test on the pool alone (W3 pooled's procedure, no
  HG002 labels);
- **target**: Learn then Test on the `k` HG002 units alone;
- **checked**: the target procedure with its grid cut at the pool's
  threshold, so it is never looser than the pool (and abstains where the
  pool does);
- **merged**: pool and `k` units in one calibration. Not valid under shift;
  it shows what merging does. Run for the first seed of each `k` only.

Each threshold is scored on all pairs of HG002's even chromosomes: whether
it certifies, the trusted share, and whether the wrong-unit share exceeds
`beta` at the Bonferroni level over all cells of valid methods.

## Gate

Fixed before looking:

- **(a) the floor holds:** at `k` = 300 and 1,000, `target` and `checked`
  certify no group at `beta = 1e-3` in any seed (the floor says they
  cannot);
- **(b) the floor is reachable:** with all odd units, `target` or `checked`
  certifies `easy` at `beta = 1e-3` in every depth and phaser;
- **(c) valid:** no violating cell for `target` or `checked` over all `k`.

R2 "passes" when (a), (b) and (c) hold: then the honest statement is that a
new batch needs about 3,000 of its own labelled units per group at
`beta = 1e-3`, which one trio on that batch provides for `easy` and
`difficult` and never for `segdup`. It does not give a few-label method.

## Code

`benchmarks/r2_batch_labels.py`; `tests/test_r2_batch_labels.py` checks the
floor against `learn_then_test`.
