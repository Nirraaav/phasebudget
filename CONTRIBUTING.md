# Contributing to PhaseBudget

PhaseBudget studies a trade: haplotype assembly is exactly solvable only at
bounded sequencing coverage, and bounding the coverage throws away reads. The
objective, the selection problem, and the accuracy metric are three different
things, and most of the results in this repository come from keeping them
apart. The repository map is in `README.md`.

## Correctness contract

1. Variant indices are global. `FragmentMatrix.subset` keeps `n_variants`, so a
   boundary index means the same thing before and after selection.
2. A read's span is `[start, end)` and is trimmed to its called alleles. An
   internal gap carries weight zero and still occupies coverage.
3. The diploid model is complementary: the second haplotype is the bitwise
   complement of the first. The general bipartition model exists in
   `exact.solve_bipartition_bruteforce` as a comparison, not as the contract.
4. Optimal solutions are not unique. Compare costs, never haplotypes, unless
   the instance is error-free.
5. A margin is a difference of two constrained optima, and its definition is
   `exact.bruteforce_margins`. Anything faster answers to that.
6. `dp`, `ilp`, and `exact` must agree on every instance all three can solve.
7. Coverage-constrained selection is exact. A heuristic may choose reads; it
   may not be described as solving the selection problem.

## Adding an optimization

1. State what it optimizes, in `docs/design/`, precisely enough to be wrong.
2. Add an independent oracle, or reuse the one in `exact`.
3. Add fixed examples and a randomized differential test against the oracle.
4. Implement it.
5. Run the full suite (`python3 -m pytest`), then measure. Accuracy and
   objective value are separate columns and must stay separate.
6. Report the numbers, including the ones that went the wrong way.

Use Python 3.10 or newer, a `src/` layout, type annotations, immutable result
objects, and small public APIs. New dependencies need a reason; the library
core has none.

## Measurement rules

A benchmark that reports accuracy must also report what the method optimized,
so the two can be compared. A benchmark that reports a ratio must report the
worst case as well as the mean. A run whose parameters are not in the printed
header did not happen.

Selection objectives are only comparable under the same weights. The adaptive
loop reweights, so its selections are scored under the base weights when it is
compared with anything else.

A test of a claim states its data, thresholds and decision rule in
`docs/design/` before the test data are read; a later change is a dated
amendment that says whether results had been seen.

## Claims about the literature

Every row of `docs/literature/01-prior-art-matrix.md` carries a status:
`recalled`, `abstract`, or `paper`. A `recalled` row may not support a
novelty claim, and no claim may depend on a negative claim about the
literature.
