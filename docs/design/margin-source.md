# Where the margin signal is measured

## The confound

The adaptive loop of `phasebudget.adaptive` computes margins on the subset that
selection just chose. That subset is the thing selection changed, so a small
margin has two possible causes and the loop cannot tell them apart:

1. the phase relation at that boundary is genuinely uncertain given the data,
   which is the signal the method wants; or
2. the selected subset happens to be thin at that boundary, so the objective is
   nearly indifferent there whatever the full read set says.

Under cause 2 the loop is a coverage-repair procedure with a margin-shaped
interface. It would still raise accuracy, and the accuracy gain reported in
The development notebook would still be real, but the stated mechanism would be wrong. This
is OQ-10.

## The two sources

`adaptive_selection` takes `margin_source`:

- `"subset"` computes margins on the selected reads. This is the deployable
  source and the only one whose cost the loop can afford, because the whole
  reason for selection is that the dynamic program is exponential in active
  coverage.
- `"full"` computes margins once on the entire fragment matrix and reuses them
  in every round. This is an oracle. It is available only when the full matrix
  is under the state limit, which restricts it to instances small enough that
  selection was not needed in the first place.

Boundary indices mean the same thing under both sources because
`FragmentMatrix.subset` keeps `n_variants`, which is correctness contract 1.
Without that a full-matrix margin could not be compared with a subset margin at
all.

`Margin.spanned` does differ between the sources, and the difference is the
confound made concrete: a boundary spanned in the full matrix and unspanned in
the subset has a structural zero margin under `"subset"` and a real margin
under `"full"`. `fragile_boundaries` drops unspanned boundaries, so such a
boundary is invisible to the subset loop and visible to the full one.

## What the full-margin loop does differently

A full-matrix margin does not depend on the selection, so it does not change
between rounds. The fragile set is therefore constant and the same reads are
boosted every round, their weights growing geometrically until the selection
stops moving, at which point the loop stops on its repeat check.

This was first written down here as the prediction that the full source is a
one-shot reweighting rather than an iterative one. That prediction is wrong.
Measured on the instances of `benchmark_margin_confound.py`, the full source
runs 2.85 to 3.20 rounds of a 4-round cap against 3.35 to 3.90 for the subset
source, and its selection moves in 15 to 18 of 20 instances. A constant
fragile set still admits several distinct selections, because the boost is
multiplicative and each round pushes a different read past the coverage cap.

The round counts do differ, though, and more rounds is a weak advantage on its
own: the loop keeps the best round it has seen. Any comparison between the two
sources therefore has to hold the round budget fixed, which is why the
benchmark reports a matched pair at one reweight each next to the full-cap
numbers.

## What settles OQ-10

Three measurements, in `benchmarks/benchmark_margin_confound.py`:

1. **Agreement.** Rank correlation between subset margins and full margins at
   the same boundary. If subset margins are a noisy copy of full margins, the
   confound is benign whatever else is true.
2. **Informativeness.** The rank-sum probability that a boundary carrying a
   switch error has a smaller margin than one that does not, computed for both
   sources against the same switch positions. This is the statistic of
   `benchmark_objective_truth.py`, applied to two signals on one instance.
3. **Accuracy.** Switch error of the adaptive loop under each source, against
   static selection on the same instances. If `"full"` beats `"subset"`, the
   subset margin is a degraded signal and the loop is leaving accuracy behind.
   If they tie, the subset margin carries what the loop needs. If `"subset"`
   beats `"full"`, the gain is not coming from the fragility signal as stated,
   and cause 2 above is the live explanation.

Outcome 3 is the one that would force a rewrite of the component-3 claim, so it
is stated here before the measurement rather than after.
