# A registered prediction of the ONT regime map from a vote model

Written on 30 September 2026, while the registered ONT run
(`ont-run-plan.md`) was running and before any ONT phasing of a
confirmatory chromosome had been compared with a truth. The commit that adds
this file fixes the model, its one parameter, and the rule for judging it.

## Amendments

1. **30 September 2026, before any prediction was reported and before the
   ONT run was scored.**
   - *What prompted it.* An exploratory HiFi table
     (`2026-09-30_spanning_informative_hifi.txt`) finished after the model
     was registered. At cap 3 and full depth, WhatsHap's selection leaves
     642 boundaries that no selected read calls on both sides, yet the
     phasing joins them through reads with gaps. Those boundaries are
     wrong at 22 per cent, and they hold 142 of its 320 scored errors. The
     exact selection has 25 such boundaries.
   - *The problem.* The registered model counts only boundaries with
     `i >= 1`, so it cannot see this component.
   - *The amended model.* It scores every boundary that both phasings join
     in one block, read from the stored run files, not from a truth. It uses
     `f(i)` where `i >= 1`, and a fixed rate `g0 = 0.215` where `i = 0`.
     `g0` is the pooled HiFi rate at gap-linked boundaries at caps 3 and 5,
     full depth: 150 wrong in 667 at cap 3, 20 in 122 at cap 5, so 170 in
     789.
   - *Both are reported.* The registered model and the amended one are both
     judged by the rule below; neither replaces the other.

2. **30 September 2026, still before any ONT scoring.**
   - *The artifact.* Amendment 1's gap-linked boundaries were mostly an
     artifact of how the library reported blocks. A variant that no
     selected read calls was still given an allele. It now is not
     (`33b92bb`).
   - *The consequence.* Its 50 per cent error rate made `g0`. With the fix,
     such boundaries leave the phasings, so the amended model is withdrawn
     and not computed.
   - *What stands.* The registered model stands unchanged. It reads only
     the selections, not the phasings, so the artifact does not touch it.
   - *Its predictions on 9 confirmatory chromosomes*, computed before any
     ONT scoring, are in
     `benchmarks/results/runs/2026-09-30_ont_prediction_9chrom.txt`. The
     pooled prediction over all 21 chromosomes follows once their fragments
     exist.

3. **30 September 2026, before any ONT scoring: baselines and a stricter
   score.** R1 and R2 compare the model with nothing. Three baselines are
   added, each defined here and computed before the ONT scoring:
   - **B0, transfer.** The corrected HiFi observed difference (v5.0q, after
     genome-run amendment 3) at the same cap. Depths are matched by nominal
     depth: ONT full to HiFi full, ONT 30 to HiFi full (about 33-fold), ONT
     20 to HiFi 20, ONT 10 to HiFi 10.
   - **B1, scaled transfer.** B0 times `q_ONT / q_HiFi = 0.010009 / 0.000260
     = 38.5`.
   - **B2, null.** 0 in every cell.

   The score is computed over the 16 cells (4 caps, 4 depths):
   - the mean absolute error against the observed difference;
   - sign agreement, over the cells whose observed interval excludes 0.

   **R3.** The model beats the baselines if its mean absolute error is below
   all three, and its sign agreement is at least theirs. R1, R2 and R3 are
   all reported, together with every cell of every predictor.

## What is predicted

O7 of the ONT plan: WhatsHap's selection minus the exact selection, in
switch errors per 10,000 common pairs, at caps 3, 5, 10 and 15 and depths
full, 30, 20 and 10, pooled over the 21 confirmatory chromosomes, against
GIAB v5.0q.

## The model

At the boundary between adjacent variants `b` and `b + 1`, let `i_b` be the
number of selected reads that call both variants. The model takes the
boundary's phase relation to be decided by a majority vote of those reads:
- each read calls the relation wrong independently, with probability
  `q = 2e(1 - e)` for a per-allele error `e`;
- a tie is broken at random.

The boundary is then wrong with probability

```text
f(i) = P(Binomial(i, q) > i/2) + P(Binomial(i, q) = i/2) / 2.
```

For a pair of selections, the predicted difference is the mean of
`f(i_whatshap) - f(i_optimal)` over the boundaries that both selections
call, those with `i >= 1` in each, times 10,000. The selections are
recomputed from the fragments, and the computation reads no truth.

The model ignores three things, and each is a way for it to fail:
- the DP decides relations jointly;
- errors cluster at some sites;
- common pairs skip variants the truth does not phase.

## The one parameter

`e = 0.00503`, the per-allele error of the ONT `realign` alleles on
chromosome 20. It was measured against v5.0q, using reads assigned to the
haplotype they agree with at most of their sites, with at least 3 such
sites. Chromosome 20 is the ONT run's smoke test and enters no test. Only
5 of its 44,437 sites have an allele error above 0.2, holding 0.6 per cent
of the errors, so the independence assumption is not contradicted there.
On HiFi the same measurement gives `e = 0.00013`, with 16 per cent of the
errors at 3 sites.

## The rule

Let `P` be the predicted difference and `O` the observed one, with its 95
per cent window-bootstrap interval, at each cap and depth.
- **R1, level.** At caps 3 and 5, at every depth, `P` lies inside the
  interval of `O`, or within a factor of 2 of `O`.
- **R2, regime.** At every cap and depth, the model's call and the data's
  call agree. The model calls "selection matters" when `P >= 0.5` per
  10,000; the data calls it when the lower end of `O`'s interval is above 0.

The model is supported if R1 and R2 both hold, and partly supported if only
R2 holds. Every prediction is reported, with its observed value, whatever
the outcome.

## The HiFi check

The same computation on the HiFi run, with `e = 0.00013`, is reported
beside the HiFi G10 table. It is exploratory, since HiFi was seen, and it
supports nothing.
