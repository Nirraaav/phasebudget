# Block-wise margin extraction

## What a margin is

For adjacent variants `j` and `j+1`, the phase relation is
`d_j = h_j XOR h_(j+1)`. The margin of that decision is the cost the objective
would have to pay to reverse it:

```text
margin_j = |  best cost with d_j = 0  -  best cost with d_j = 1  |.
```

The absolute value is the whole definition, because the smaller of the two is
the unconstrained optimum. `phasebudget.exact.constrained_optimum` computes one
side by enumerating haplotypes, and `bruteforce_margins` calls it twice per
boundary. That is the definition; everything below is a way to compute it
faster.

## The state

`phasebudget.dp.CoverageDP` puts one state per pair

```text
( labelled bipartition of the reads active at j , allele of h at j ).
```

The bipartition is labelled rather than unordered: bit `i` says which haplotype
the `i`-th active read belongs to. Labelling doubles the state count and
introduces a global symmetry, since flipping every label and every allele at
once gives the same cost. The symmetry is harmless here and is what makes the
margin well defined without a symmetry-breaking constraint: under the flip,
`h_j` and `h_(j+1)` both change, so `d_j` does not.

Column cost is computed once per column for all `2^k` masks using

```text
cost(mask, a=0) = sum_i c0_i + sum_(i in mask) (c1_i - c0_i)
cost(mask, a=1) = sum_i c1_i - sum_(i in mask) (c1_i - c0_i),
```

where `c0_i` and `c1_i` are the weighted cost of giving read `i` the value `0`
or `1` at this variant. The inner sum is filled in by one pass over the masks
in increasing order, so the tables cost `O(2^k)` rather than `O(k 2^k)`.

## The transition

Reads active at both `j` and `j+1` must keep their labels. Reads that end are
forgotten and reads that start are free. Grouping the states of column `j` by
their labels on the shared reads turns the transition into a lookup: a state at
`j+1` reads the best predecessor with the matching projection. The allele is
unconstrained across columns, so the projection is over labels only.

## Reading the margins off

Let `F_j(mask, a)` be the best cost of variants `0..j` ending in that state and
`G_j(mask, a)` the best cost of variants `j..m-1` starting from it. Both
include the cost of column `j`, and the transition between `j` and `j+1` adds
no cost of its own, so for any compatible pair

```text
total = F_j(mask, a) + G_(j+1)(mask', a').
```

Grouping both tables by their projection onto the shared reads and minimizing
over the projection and over `a` gives

```text
best(delta) = min over p, a of  Fproj_j(p, a) + Gproj_(j+1)(p, a XOR delta),
```

which is exactly the constrained optimum at `d_j = delta`. Two passes give
every margin.

`tests/test_dp.py::test_block_wise_margins_match_constrained_re_solving` checks
this against `bruteforce_margins` on 200 random instances with ragged spans,
interior gaps, and uneven weights. `test_boundary_costs_bracket_the_optimum`
checks the identity `min(best(0), best(1)) = wMEC` separately, since a
consistent error in both directions would pass the first test alone.

## Two facts worth stating

A boundary no read spans has margin zero. Nothing in the code special-cases it:
with no shared reads the projection is trivial, the suffix can be relabelled
freely, and the two constrained optima coincide. The block structure falls out
of the arithmetic.

A spanned boundary can also have margin zero.
`tests/test_dp.py::test_a_spanned_boundary_can_still_have_zero_margin` carries
one: a variant called by exactly one read on each side leaves its allele free,
so both relations cost the same. A zero margin therefore means the objective is
indifferent, which is a weaker statement than "no read connects these
variants". Any use of margins as a confidence signal has to keep the two cases
apart, which is why `Margin` carries `spanned` next to the number.
