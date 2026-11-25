# Posterior margins

## Why a second signal

The MEC margin of a boundary is a max-product quantity: the optimum with the
phase relation reversed, minus the optimum. It compares the best solution with
the single best alternative. When reads are noisy there are many near-optimal
alternatives, and the margin sees one of them. `benchmark_objective_truth.py`
measured its informativeness falling from 0.977 to 0.613 as the error rate
rose, and `benchmark_heterozygosity.py` measured it at 0.72 at low density
and error 0.20 (a first measurement of 0.423 was a scoring error, corrected on
26 September). The sum-product version counts every solution. ProbHap
(Kuleshov 2014, section 4.4) computes it as a transition probability; this is
the same quantity on this project's model, as a log-odds.

## Definition

With log-odds weights, `w = log((1-e)/e)`, the cost of a solution is minus
the log-likelihood ratio of the reads against error-free reads, so
`exp(-cost(h, a))` is proportional to the likelihood of haplotype `h` and read
assignment `a` under a uniform prior on both. For a boundary `b` and a relation
`r` in `{0, 1}`, let

```text
Z_r(b) = sum over (h, a) with h[b] xor h[b+1] = r of exp(-cost(h, a)).
```

The posterior margin of `b` is `|log Z_0(b) - log Z_1(b)|`, and its relation is
the `r` with the larger `Z_r`. It is a log posterior odds, in the same units as
the MEC margin: when one solution dominates both sums the two margins
coincide, and they part when many solutions are close to optimal. A boundary
no read spans has `Z_0 = Z_1` by relabelling one side, so its posterior margin
is zero, like its MEC margin.

## Computation

The forward and backward passes of `phasebudget.dp` with `(min, +)` replaced
by `(logsumexp, +)` over `-cost`, on the same states. Every `(h, a)` is exactly
one state sequence, because every read is active at some variant and its side
is fixed wherever it is active. `Z_r(b)` combines the forward table at `b`
with the backward table at `b+1` by shared-read code and allele, as
`boundary_costs` does, with sums in place of minima.

## Oracle

Summing over assignments read by read,

```text
sum over a of exp(-cost(h, a)) = product over reads r of
    (exp(-cost_r(h, side 0)) + exp(-cost_r(h, side 1))),
```

so `Z_r(b)` is a sum over the `2^n` haplotypes of that product.
`exact.bruteforce_posterior_margins` enumerates it, and
`tests/test_posterior.py` requires agreement on random instances.

## What would make it the better signal

It is better if it separates boundaries carrying a switch error from those
that do not with a higher rank probability than the MEC margin, at the error
rates where the MEC margin fails, and if the adaptive loop steered by it
lowers switch error against the loop steered by margins at the same cap,
boost, quantile and iteration count. Either outcome is reported.
