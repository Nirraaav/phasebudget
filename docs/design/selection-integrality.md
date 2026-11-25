# Why coverage-constrained read selection is easy

## The program

Reads `r = 0, ..., n-1` occupy intervals `[l_r, u_r]` over variants
`0, ..., m-1` and carry scores `w_r >= 0`. For a cap `c`,

```text
maximize    sum_r w_r x_r
subject to  sum_(r : l_r <= j <= u_r) x_r <= c    for j = 0, ..., m-1
            x_r in {0, 1}.
```

## Integrality

Write the constraint matrix `A` with `A_(j,r) = 1` when read `r` covers variant
`j`. Because a read covers a contiguous run of variants, the ones in every
column of `A` are consecutive. Arkin and Silverberg 1987 prove, in the Lemma
of their section 2, that a matrix of this form is totally unimodular; theirs is
the clique-by-job matrix of an interval graph, which is `A` with maximal
cliques in place of variants. The row form of the same fact is Corollary 8.4c
of Schrijver's lecture notes *A Course in Combinatorial Optimization*, and
total unimodularity is closed under transposition. Every vertex of `{x : Ax <= c1, 0 <= x <= 1}` is
therefore integral, and solving the linear relaxation solves the integer
program.

That is a statement about the polytope, not an algorithm. The algorithm comes
from the flow formulation below, which is what the code uses.

## The flow network

Nodes are the `m+1` variant boundaries `0, ..., m`. Arcs are:

* a **slack arc** `j -> j+1` of capacity `c` and cost `0`, for each `j < m`;
* a **read arc** `l_r -> u_r + 1` of capacity `1` and cost `-w_r`, for each `r`.

Push `c` units from node `0` to node `m` at minimum cost. Reads whose arc
carries a unit are selected. This is the network of Lemma 4 of Carlisle and
Lloyd 1995, *On the k-coloring of intervals*, with variants in place of interval
endpoints; their Theorem 3 gives the `O(c S(n))` bound.

**The cut argument.** The arcs crossing the cut between nodes `j` and `j+1` are
the slack arc and exactly the read arcs with `l_r <= j <= u_r`, which is to say
the reads covering variant `j`. Every unit of the `c` units crosses that cut
once, so the number of selected reads covering variant `j` is `c` minus the
flow on the slack arc, which is at most `c`. Feasibility of the flow and
feasibility of the selection are the same condition.

**The converse.** Given a feasible selection `S`, put one unit on each selected
read arc and `c - |{r in S : l_r <= j <= u_r}|` units on slack arc `j`. Flow
conservation at node `j` holds because the coverage at `j` differs from the
coverage at `j-1` by exactly the reads starting at `j` minus the reads ending
before `j`. The cost of that flow is `-sum_(r in S) w_r`.

The two directions together say the minimum-cost flow of value `c` is the
negative of the maximum-weight feasible selection. Since the network has
integer capacities, an integral optimum exists and successive shortest paths
finds one.

## Implementation

`phasebudget.selection.select_optimal` builds the network above and runs
successive shortest paths with potentials. Every arc runs from a lower node
index to a higher one, so the initial potentials come from one relaxation pass
in node order rather than from Bellman-Ford, and the negative read costs are
handled without a negative cycle ever appearing. Augmentation is by the
bottleneck capacity, so at most `c` iterations run.

Reads with `w_r <= 0` are left out of the network. Including them would not
change the objective, since the slack arcs can carry the flow instead, but a
zero-weight read admitted on a tie would still consume coverage that a later
stage might have used.

Complexity is `O(c * (n + m) log(n + m))` with a binary heap, against
`O(n log n + n m)` for the greedy baseline. The measured cost is in
`benchmarks/benchmark_selection.py`; at 200 variants and 40-fold coverage the
exact selector takes single-digit milliseconds.

## What this does not cover

The argument needs the objective to be additive over reads. Three plausible
objectives are not:

1. **Redundancy-sensitive scores.** Paying once for covering a boundary rather
   than once per read that covers it turns the objective into a coverage
   function, and maximizing a submodular function under these constraints is no
   longer obviously in P.
2. **Connectivity requirements.** Requiring the selected reads to leave the
   variant graph connected adds a constraint that is not an interval
   constraint. The published WhatsHap selection procedure attends to
   connectivity explicitly, which is one reason its objective is not the one
   solved here.
3. **Position-dependent caps.** A cap `c_j` that varies by variant keeps the
   flow argument intact, since the slack arc capacity may vary. This one is
   still easy and is worth implementing.

Item 1 is the optional theoretical extension in the proposal. Nothing in the
deliverable depends on its resolution.
