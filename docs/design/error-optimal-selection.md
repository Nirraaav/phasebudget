# Read selection against a model of switch error

Direction D1 of the research plan (development repository), section 13. Written on 27 September
2026, before the selector existed.

## The model

Boundary `b` lies between variants `b` and `b + 1`. A selected read spans
it when its span contains both, and each spanning read votes on the phase
relation across `b`. A read calls the relation wrong when exactly one of its
two alleles is wrong, so with independent allele errors at rate `e` its
relation error is `q = 2e(1 - e)`. If the relation were decided by a vote of
the `s` spanning reads alone, its error probability would be a function
`f(s)` of the count only.

The exact majority error is not convex in `s`: with ties split at random,
`f(1) = f(2) = q`. The model therefore uses the Bhattacharyya bound on a
vote of `s` independent readings,

```text
f(s) = rho^s,   rho = 2 sqrt(q (1 - q)),
```

which is convex and decreasing, bounds the exact error from above, and has
the same exponential rate in `s`. At `e = 0.00117`, the pooled allele error
of the recalibration fit on the odd chromosomes (`data/genome/recalibration.json`), `rho` is about
0.097.

The model ignores three things, and each is a way for it to be wrong:
- the phasing DP decides every relation jointly, and reads that do not span
  `b` constrain it through their neighbours;
- read errors are not independent or equal;
- the relation of a boundary no read spans is not decided at all. It ends
  the block, which is a different error from a switch.

## What the selector optimizes

For every boundary `b` that some positively weighted read spans, `s_b` is
the number of selected reads spanning it. Given a cost sequence
`g(0), ..., g(c)` that is non-increasing and convex,

```text
maximize    sum_r w_r x_r  -  sum_b g(s_b)
subject to  sum_{r covers j} x_r  <= c     for every variant j
            x_r in {0, 1}
```

The model's choice is `g(s) = lambda * rho^s`. `lambda` converts expected
switches into read weight; with log-odds weights, one unit is one nat of
evidence. `select_connected` is the case `g(s) = penalty * max(k - s, 0)`,
which is also non-increasing and convex.

## Why it is exact

Use the split network of `connectivity-selection.md`. The flow on boundary
`b`'s arcs is `u_b = c - s_b`, so a cost `g(s_b)` is a cost `h(u_b) =
g(c - u_b)` on that flow. It is convex in `u_b`, because it is `g` composed
with an affine map, and non-decreasing, because `g` is non-increasing.

Replace the boundary arc by `c` parallel arcs of capacity 1. Give the `u`-th
arc the cost `h(u) - h(u - 1)`. Convexity makes these costs non-decreasing
in `u`, so a min-cost flow fills the cheaper arcs first. The cost of any
integral flow is then `h(u_b) - h(0)`, and `h(0) = g(c)` is a constant.

The constraint matrix is the interval matrix of the connected selector, and
the capacities are integers, so the min-cost flow is an integral optimum of
the program. A cost sequence that is not convex breaks the argument, and
the selector rejects it.

## Oracle and tests

`tests/test_error_optimal_selection.py`:
- enumerates every subset of small random matrices and compares the
  program's optimum with the selector's, for random convex non-increasing
  cost sequences and for the model's `lambda * rho^s`;
- checks that `select_connected` equals the selector with the step cost;
- checks that a non-convex sequence is rejected.

## What would make it useful, and what would not

It is useful if, at a fixed cap, the error-optimal selection has fewer
switch errors than `optimal` on real reads it was not tuned on. It is not
useful if it only lengthens blocks, which `connected` and two-level
orientation already do, or if its gain disappears once blocks are oriented.
The pilot tunes `lambda` on the odd chromosomes and tests on the even ones.

## Pilot protocol, fixed before the pilot ran

- **Data.** The `realign` fragments of the genome-wide HiFi run at all four
  depths, with log-odds weights. HiFi has already been seen, so the pilot
  is exploratory and supports no claim beyond deciding whether D1 goes on.
- **Model.** `switch_model_costs` with allele error `e = 0.00117`, the
  pooled rate of `recalibration.json` (odd chromosomes). Scales `lambda` of
  `1e2`, `1e4` and `1e6` nats per expected switch, at caps 3 and 5.
- **Phasing.** As for `optimal`: `CoverageDP` on the selection.
- **Tuning.** On the odd confirmatory chromosomes, the scale with the lowest
  pooled switch error against the filtered truth, pooled over caps 3 and 5
  and all four depths, wins; a tie goes to the smaller scale.
- **Test.** On the even confirmatory chromosomes, that scale is compared
  with `optimal` by the paired sign test on common pairs, at each cap and
  depth. Pooled switch error, pairs and block N50 are reported beside
  `optimal`, `connected` and `whatshap`.
- **Gate.** D1 goes on if the error-optimal selection has fewer discordant
  errors than `optimal` at cap 3 at every depth on the even chromosomes,
  and is significant at 0.05 in total. The total is the sign test on the
  discordant pairs summed over the four depths.
