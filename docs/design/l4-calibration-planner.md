# L4: how many calibration units a group needs to certify

Round 8, idea 4. Written 22 November 2026, before any subsampling was run.
Development analysis: both pools (K1b's six genomes, W1b's 21) have been
read and scored before, so nothing here is blind.

## What Learn then Test needs

`learn_then_test` (`src/phasebudget/confidence.py`) scans the coverage grid
from strict to loose and stops at the first threshold whose exact binomial
p-value `P(X <= k | n, beta)` exceeds `delta`. A group therefore certifies
something if and only if the first grid threshold passes. That threshold,
`t0`, is the score of the pair at rank `round(0.5 N)` among the group's `N`
pairs ranked by score: half the pairs are trusted there. At `t0` let `n` be
the units holding a trusted pair and `k` the units holding a trusted wrong
pair.

For fixed `n`, `beta` and `delta`, define

    k*(n) = max { k >= 0 : P(Bin(n, beta) <= k) <= delta },  or -1 if none.

The group certifies if and only if `k <= k*(n)`. If the `n` units are
independent and each holds a trusted wrong pair with probability `r`, then

    P(certify) = P(Bin(n, r) <= k*(n)).                               (1)

This is exact for the test as implemented, given `n` and `r`. It is not
monotone in `n` (`k*` steps up by one at a time), so "units needed for power
`pi`" is defined through the predicate "(1) is at least `pi` at `n` and at
50 evenly spaced points of `[n, 2n]`", and `units_needed` finds the `n` where
it turns true by bisection, to 0.1%, from the `r = 0` floor upward.

Special cases:

- `r = 0`: (1) is 1 when `k*(n) >= 0`, that is `(1 - beta)^n <= delta`, so
  `n >= ln(1/delta) / -ln(1 - beta)`, about `ln(1/delta)/beta`. At
  `delta = 0.05` this is 2,995 units at `beta = 1e-3`, 998 at `3e-3` and
  9,985 at `3e-4`. This is the floor met in W3 and K1a.
- `0 < r < beta`: the normal approximation gives

      n ~ (z_{1-delta} sqrt(beta (1-beta)) + z_pi sqrt(r (1-r)))^2 / (beta - r)^2,

  which grows like `beta / (beta - r)^2` as `r` approaches `beta`. Near the
  target, halving the gap `beta - r` multiplies the units needed by four.
- `r >= beta`: (1) tends to 0 as `n` grows. No pool size certifies.

## The planner

Inputs: a group's unit rate `r` at `t0` and the units per genome `u` at
`t0`. Output: units `n(pi)` from (1) and genomes `ceil(n(pi) / u)` for
`pi` in {0.8, 0.95}, at `delta = 0.05` and `beta` in {3e-4, 1e-3, 3e-3}.
The pooled `r` and `u` are measured on each full pool. The planner assumes
independent units with a common rate. It says nothing about whether the
certified threshold is valid on a new genome (W1b's B1 failure is about
that), only about whether calibration certifies at all.

## The check

Pools, at 30x, both phasers (whatshap, LongPhase), groups `easy`,
`difficult`, `segdup`, `lowmap`, `tandem`, `indel`:

- K1b: HG002 (`hprc_HG002`) and five Pedigree HiFi genomes, agreeing pairs
  only (`benchmarks/k1b_score.py`'s rows, mode `agree`).
- W1b: HG002 and the 20 HPRC pool genomes, all pairs, consensus truth
  (`benchmarks/w1b_score.py`'s rows).

For each pool size `m` in 1 to 6, draw `B = 400` bootstrap pools of `m`
genomes with replacement (a genome drawn twice counts as two genomes with
distinct units), run the group's first grid threshold exactly as
`learn_then_test_by_group` does, and record whether it certifies. The
observed frequency is compared with the predicted power: the mean over the
same draws of (1) with each draw's own `n` and the full pool's `r`.

The threshold is chosen from the data, and `r` is measured at the full
pool's `t0`, not each draw's. The prediction is therefore a statement that
`t0` is stable across draws and that genomes share one unit rate. Genome
to genome differences in rate (overdispersion) make the observed frequency
lower than (1) when `r < beta` is close to `beta` and higher when `r` is
above `beta`; that is the departure the check can find.

A cell is (pool, phaser, group, beta, `m`). It agrees when the predicted
power lies inside the 95% Wilson interval of the observed frequency over
the 400 draws. A cell is informative when the predicted power is in
[0.02, 0.98] or the observed frequency is strictly between 0 and 1; cells
where both are at the same end agree trivially and are counted separately.

**Success rule.** At least 80% of informative cells agree. All cells,
informative cells, agreeing cells and the worst disagreement (largest
absolute difference between observed and predicted) are printed. With 95%
intervals, about 5% of cells are expected to disagree even when (1) is
exact.

The implementation of (1) is checked against brute force over `k` in
`tests/test_calibration_planner.py`, and the fast first-threshold check in
the benchmark is compared with `learn_then_test_by_group` on a sample of
draws before any frequency is reported.
