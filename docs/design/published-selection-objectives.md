# Published read-selection objectives

The three selectors in `phasebudget.selection` do not solve the same
optimization problem. This note fixes the comparison before any benchmark is
interpreted.

## Fischer and Marschall 2016

Fischer and Marschall's [preprint](https://doi.org/10.1101/046771) first
removes reads that call fewer than two heterozygous variants. It then builds
disjoint slices. For a read `r`, let `V(r)` be its called variants, `H(r)` its
uncalled positions strictly inside its span, and `V(S)` the variants covered
by the reads already chosen for the current slice. The priority key, in
descending lexicographic order, is

```text
( |V(r) \ V(S)| - |H(r)|, |V(r)| - |H(r)|, min_{v in V(r)} quality(r, v) ).
```

The first phase admits a feasible read only when it covers a variant new to
the current slice. The optional second phase uses union-find components induced
by the first phase and adds feasible remaining reads that join two components.
The stated goals are broad coverage, connectivity, independent connections,
many variants per read, and high read quality. They are not an additive
per-read objective.

Algorithm 1 says slices are disjoint and its prose says reads move from the
undecided set when selected or discarded. Its bridging pseudocode adds bridge
reads but omits the corresponding removal from the undecided set. The
implementation removes those reads after bridging. Otherwise a later slice
can select a bridge read again, contradicting the stated disjointness and the
coverage cap.

`select_fischer_marschall` implements this procedure. Its `objective` field is
still the caller-supplied additive score because `SelectionResult` needs a
common reporting field. That field is a cross-objective measurement, not an
optimality certificate for this selector. `distinct_variants_covered` and
`component_count` report two outcomes that the procedure is designed to
affect.

## What WhatsHap 2.8 runs

The procedure above is the paper's. The code WhatsHap ships,
`whatshap/readselect.pyx` in release 2.8, departs from it, and on 1,500
simulated instance and cap pairs the paper's procedure selected the same
reads as WhatsHap in 178. `select_whatshap` is a line-by-line port of the
shipped code and selected the same reads as WhatsHap in all 1,500 pairs and
in 36 of 36 chromosome-scale pairs (2,000 variants, caps 5, 10 and 15);
`tests/test_whatshap_port.py` holds 160 of those pairs as a fixture and
checks them in CI. The departures:

1. **The first score falls by the wrong count.** When a read is admitted to
   a slice, every queued read that shares a newly covered variant has its
   first score lowered by `whatshap_score_drop`: one for each of its variants
   that the admitted read did not newly cover. The paper's `|V(r) \ V(S)|`
   would lower it by the variants that were newly covered. The docstring of
   `_compute_score_for_read` describes the paper's score, so the code is
   probably not doing what its authors meant. The effect is to push down
   queued reads that extend past the admitted read, which are the reads that
   would cover new variants.
2. **Scores change inside the heap.** WhatsHap's `PriorityQueue` is a binary
   max-heap with an index; a changed score is sifted in place and ties
   leave in heap order, which depends on insertion order and on the order of
   earlier changes. That order in turn depends on iterating a C++
   `unordered_set<int>` of genomic positions. The port reproduces both,
   using libstdc++'s bucket growth measured with g++ 13.3.
3. **Bridging checks the cap first.** In the bridging phase a read that would
   breach the cap is removed from the undecided set before the component
   test, so it is gone for later slices even if it bridges nothing now. The
   paper does not say.
4. **Gaps and coverage are counted over present positions.** Holes are
   counted over the positions some candidate read calls, not over all
   variants. This changes nothing when every variant is called by some read.

Departure 1 is the one that matters. On 2,000-variant ONT-like instances at
cap 3, the shipped rule left 22 to 27 phase blocks where the paper's rule,
substituted into the port and nothing else changed, left 1 or 2.
`benchmarks/diagnose_whatshap_score_rule.py` measures the cost in switch
error and block length; the development notebook has the numbers.

Every benchmark in this repository now reports both selectors:
`fischer_marschall` is the paper, `whatshap` is the shipped code.

## Mäkinen, Staneva, Tomescu, and Valenzuela 2015

Mäkinen et al.'s [paper](https://arxiv.org/abs/1508.07820) solves a different
interval problem. Given a cap `c`, choose a feasible read set `S` that
maximizes

```text
min_j coverage_S(j), subject to max_j coverage_S(j) <= c.
```

For a candidate floor `t`, their Theorem 1 reduces feasibility to a flow
network. Interior backbone arcs have capacity `c - t`, the two end arcs have
capacity `c`, and every read supplies a unit-capacity interval arc. A flow of
value `c` exists exactly when a selection has maximum coverage at most `c` and
minimum coverage at least `t`. Binary search over `t` gives the optimum floor.

`select_max_min_coverage` uses that reduction. It can use the additive score
only to break ties after it has fixed the best floor. The default excludes
reads that are not phase informative, matching the project's wMEC contract.
Passing `informative_only=False` implements the interval problem as stated in
the paper, including one-variant intervals.

## Comparison rule

`select_optimal` maximizes additive weight, `select_max_min_coverage` maximizes
the coverage floor, and `select_fischer_marschall` is a heuristic for coverage
and connectivity. A ratio to the additive optimum is meaningful only for the
first objective. The selection benchmark therefore prints additive weight,
minimum coverage, covered variants, component count, read count, and runtime
for every selector.
