# The coverage convention

A read that calls alleles at variants 3 and 7 and nothing in between occupies
five variants or two, depending on what "occupies" means. Every coverage
number in this repository, every feasibility check, and the tractability of
the selection problem depend on which of the two is meant. This page fixes it
and gives the reason.

## 1. The two conventions

Let read `r` call alleles at the set `C_r` of variants, and let its span be
`[l_r, u_r]` with `l_r = min C_r` and `u_r = max C_r`.

**Span convention.** Read `r` consumes coverage at every variant in
`[l_r, u_r]`, holes included. Per-variant coverage is

```text
cov_span(j) = |{ r selected : l_r <= j <= u_r }|.
```

**Called convention.** Read `r` consumes coverage only where it calls an
allele. Per-variant coverage is

```text
cov_called(j) = |{ r selected : j in C_r }|.
```

The span convention is the contract. `FragmentMatrix.coverage`,
`FragmentMatrix.active`, `select_optimal`, and every benchmark in this
repository use it. `phasebudget.gapped` implements the called convention for
the theory and for the measurement in section 6, and no part of the pipeline
calls it.

## 2. Why the span convention

### It bounds the quantity the solver pays for

`CoverageDP` keeps one state per labelled bipartition of the fragments active
at a variant, where active means the span contains the variant. A read with an
internal hole is still active across the hole: it still ties the variants on
either side into one phase block, and the dynamic program still carries its
side in the state. The table at variant `j` therefore has `2^(cov_span(j)+1)`
entries, not `2^(cov_called(j)+1)`.

A cap on `cov_called` bounds a quantity the solver does not pay for. Take one
read per unit of an instance with span `[0, m)` and calls only at `0` and
`m-1`. Under the called convention `k` such reads give every interior variant
coverage zero, so any cap admits all of them, and the dynamic program then
holds `2^(k+1)` states at every interior variant. The called convention is not
merely harder to solve, it constrains the wrong thing.

This is the decisive argument and it is checkable here rather than in a paper:
`CoverageDP` raises `CoverageLimitExceeded` off `matrix.max_coverage()`, which
is span coverage.

### It keeps the selection problem in P

Under the span convention the constraint matrix has consecutive ones in every
column, so it is an interval matrix, so it is totally unimodular, so the
linear relaxation is integral and a min-cost flow solves the problem.
`docs/design/selection-integrality.md` writes the flow out. Sections 4 and 5
below show that all of that fails under the called convention.

### It agrees with WhatsHap's accounting

Patterson et al. 2015 define `F(j)` as the fragments covering column `j` and
store `2^(|F(j)|-1)` entries for it. The deployed selector agrees, and this is
checked against source rather than read off a paper: in `whatshap/readselect.pyx`
at commit `1f92cf6` (22 Sep 2026) a read is admitted only if
`coverages.max_coverage_in_range(begin, end) < max_cov`, where `begin` and `end`
are its first and last variant, and `CovMonitor.add_read` in
`whatshap/coverage.py` increments every position of `range(begin, end)`. That
is span coverage, holes included. This is a consistency note, not the reason:
the reason is the previous two subsections.

## 3. The two conventions are ordered

`cov_called(j) <= cov_span(j)` at every variant and for every selection,
because `C_r` is a subset of `[l_r, u_r]`. So a selection feasible under the
span convention is feasible under the called convention, and

```text
OPT_span(c) <= OPT_called(c)
```

for every cap `c` and every non-negative additive score. The span convention
is the conservative one. It never admits a selection the dynamic program
cannot afford, and it does refuse selections a called-convention solver would
accept. Section 6 measures how much weight that refusal costs.

## 4. The called convention breaks integrality

Three reads over three variants, cap `c = 1`, unit weights:

```text
read 0:  1 1 -
read 1:  1 - 1
read 2:  - 1 1
```

Read 1 has a hole at variant 1. Under the called convention the constraints
are `x0 + x1 <= 1` at variant 0, `x0 + x2 <= 1` at variant 1, and
`x1 + x2 <= 1` at variant 2. Any two reads collide, so the integer optimum is
1. The point `x = (1/2, 1/2, 1/2)` satisfies all three constraints and scores
3/2.

The linear relaxation is therefore strictly above the integer optimum, the
constraint matrix is not totally unimodular, and no flow argument of the kind
in `docs/design/selection-integrality.md` is available.
`tests/test_gapped.py` solves both programs with PuLP and checks the two
values.

This instance is the triangle. It is the smallest case of the general
reduction below, and it is a hole in read 1 that creates it: under the span
convention read 1 covers variant 1, all three reads are then active there, and
the cap already forbids any two of them.

## 5. The called convention is NP-hard

**Claim.** Deciding whether a selection of total weight at least `W` exists
under the called convention, given a cap `c`, is NP-hard.

**Reduction.** Let `G = (V, E)` be a graph with no isolated vertex, an
instance of maximum independent set. Index the edges `0, ..., |E|-1`. Build
one variant per edge and one read per vertex: read `v` calls allele `0` at
every variant whose edge is incident to `v`, with weight 1, and has holes
everywhere else inside its span. Set `c = 1` and `W = k`.

Every variant `j` is an edge `{u, v}`, and the reads calling `j` are exactly
`u` and `v`. The constraint at `j` is `x_u + x_v <= 1`. A selection is
feasible exactly when the chosen vertices are pairwise non-adjacent, and its
weight is the number of them. So `G` has an independent set of size `k`
exactly when the selection instance has a feasible selection of weight `k`.
The construction is linear in `|V| + |E|`, and maximum independent set is
NP-hard (Karp 1972).

`phasebudget.gapped.independent_set_instance` builds the reduction and
`tests/test_gapped.py` checks the equivalence against a brute-force maximum
independent set on every graph up to five vertices, both directions: optimal
weight equals independence number, and every feasible selection is an
independent set.

**How many holes this needs: one.** The reduction above needs up to
`deg(v) - 1` holes per read. The literature settles the rest.
Bar-Yehuda, Halldórsson, Naor, Shachnai and Shapira, "Scheduling split
intervals", SIAM J. Comput. 36(1):1-15, 2006, Corollary 2.3: maximum
independent set is APX-hard on unweighted (2,2)-union graphs, where every
vertex is one segment of length 2 on each of two lines. Lay the second line to
the right of the first. Each vertex becomes a read that calls two adjacent
variants, skips a gap, and calls two more, which is exactly one hole run, and
two reads share a called variant exactly when the vertices conflict.
`phasebudget.gapped.two_union_instance` builds the embedding and
`tests/test_gapped.py` checks conflict preservation on 300 random instances.

So the called convention is APX-hard, not merely NP-hard, already when every
read has one hole and calls four alleles. Earlier NP-hardness for 2-interval
independent set is Bafna, Narayanan and Ravi 1996 and Vialette 2004, cited
through Bar-Yehuda et al. The hardness at cap 1 carries to every fixed cap
`c`: add `c - 1` hole-free reads spanning everything with weight larger than
the total, which any optimum keeps and which use up `c - 1` of the cap.

This means the hardness result in this section is a short corollary of known
results, not a theorem this project owns. The earlier version of this page
called the one-hole case open; it was not, and the review that found this is
recorded in the development notebook under 24 Sep.

## 6. What is left

**Algorithms, not complexity.** Bar-Yehuda et al. give a `2t`-approximation
for independent set in `t`-interval graphs, which is the cap-1 case with
`t - 1` holes per read, and "Scheduling split intervals with non-uniform
demands" (Discrete Optimization, 2020, status `recalled`) extends that line
towards capacities. Whether any of it beats the span convention's exact
answer in weight on phasing-shaped instances is a measurement, and it is not
taken.

**Unimplemented, not open: recognizing the easy instances.** Total
unimodularity of a 0-1 matrix is decidable in polynomial time by Seymour's
decomposition theorem (status `recalled`), and consecutive ones under some
permutation by PQ-trees (Booth and Lueker 1976, status `recalled`). Neither is
implemented and no deliverable needs either.

**Measured on 26 September 2026.** By section 3 the span convention leaves
`OPT_called(c) - OPT_span(c)` of weight on the table.
`benchmarks/benchmark_convention_cost.py`, 60 variants, depth 20, span 6,
error 0.05, 20 replicates per cell, called optimum by 0-1 program under CBC:

| gap rate | hole share | cap 3 | cap 5 | cap 8 | span coverage of the called optimum at cap 8, mean and max |
|---:|---:|---:|---:|---:|---|
| 0.10 | 0.066 | 0.983 | 0.974 | 0.963 | 10.1, 11 |
| 0.20 | 0.125 | 0.960 | 0.948 | 0.934 | 11.0, 12 |
| 0.30 | 0.194 | 0.924 | 0.910 | 0.893 | 12.2, 14 |

The columns are `OPT_span / OPT_called`; worst replicates are 0.944, 0.916 and
0.866 at cap 8. The span convention gives up 2 to 11 per cent of the weight as
the share of holes rises to a fifth. The called optimum does not come free: it
puts up to 14 spanning reads over one variant at cap 8, and the dynamic program
pays for span coverage, so phasing it exactly costs the state count of cap 14.
The weight the span convention gives up is weight the solver could not afford
at the same cap.

## 7. The rule this leaves

1. Coverage means span coverage. `FragmentMatrix.coverage` is the definition.
2. A hole carries weight zero, constrains no allele, and still consumes
   coverage. The three are independent statements and all three hold.
3. The called convention lives in `phasebudget.gapped`. It exists to support
   sections 4 to 6 and is not part of any pipeline. A benchmark using it says
   so in its header.
