# Problem statement: coverage-constrained haplotype assembly

## 1. Objects

A diploid individual carries two copies of each chromosome. At a heterozygous
site the two copies differ, and after genotyping the two alleles are labelled
`0` and `1` without knowing which copy carries which. Phasing is the assignment
of alleles to copies across many sites at once.

Let the heterozygous variants be indexed `0, ..., m-1` in reference order. A
sequencing read that covers some of them contributes a row of the **fragment
matrix**

```text
R in {0, 1, -}^(n x m),
```

with `-` for a variant the read does not call. Row `r` carries a weight
`w_(r,j) >= 0` for every called allele, the confidence that the call is
correct. The standard choice is the log-odds weight
`w = log((1 - e) / e)` for a base with error probability `e`, which is zero for
a coin flip and grows without bound as the call becomes certain.

Every read is trimmed to its called alleles, so its span is the interval

```text
[l_r, u_r],  l_r = min {j : R_(r,j) != -},  u_r = max {j : R_(r,j) != -}.
```

A `-` strictly inside the span is an internal gap. The span, gaps included, is
what consumes coverage: a read that skips a variant in the middle still ties
the variants on either side of it together, and the dynamic program still
carries its state there. `docs/design/coverage-convention.md` states the
alternative convention, shows that it makes the linear relaxation fractional
and the selection problem NP-hard, and fixes this one as the contract.

## 2. The objective

Under the diploid model at heterozygous sites the two haplotypes are
complementary, so a solution is one binary string `h in {0,1}^m` together with
a side `s_r in {0,1}` for every read. Weighted minimum error correction is

```text
minimize   sum_r sum_(j: R_(r,j) != -) w_(r,j) * [ R_(r,j) != h_j XOR s_r ].
```

For a fixed `h` the inner minimization over `s_r` is independent per read, so

```text
wMEC(R) = min_h sum_r min( d(r, h), d(r, complement of h) ).
```

`h` and its complement have equal cost. `phasebudget.exact` breaks that
symmetry by fixing `h_0 = 0`.

### 2.1 What the complementary assumption costs

The classical MEC formulation splits the reads into two sets and corrects each
to its own column consensus, with no constraint linking the two consensus
values at a variant. That model is strictly cheaper: a column where every read
agrees costs nothing under it and costs the minority weight under the
complementary model. `phasebudget.exact.solve_bipartition_bruteforce`
implements it and `tests/test_exact.py` carries an instance where the two
disagree. The complementary model is the one used here, because a site that is
called heterozygous is heterozygous by assumption, and a genotyping error at
that site is a different problem from a sequencing error.

## 3. The coverage constraint

The dynamic program in `phasebudget.dp` keeps one state per labelled
bipartition of the reads active at a variant, so its cost is exponential in the
active coverage `c` and linear in `m`. Sequencing depth is not bounded by
anything the assembler controls, so reads must be dropped before solving.
Dropping them discards information. That trade is the object of study.

For a chosen cap `c`, a selection `x in {0,1}^n` is feasible when

```text
sum_(r : j in [l_r, u_r]) x_r <= c    for every variant j.
```

For an additive score `w_r` the selection problem is

```text
maximize  sum_r w_r x_r   subject to feasibility.
```

The constraint matrix has a one in row `j` and column `r` exactly when
`l_r <= j <= u_r`, so its ones are consecutive in every column. A matrix with
the consecutive-ones property is an interval matrix and therefore totally
unimodular, so the linear relaxation has an integral optimum and the problem is
in P. `docs/design/selection-integrality.md` writes out the equivalent
min-cost-flow network, which is what `phasebudget.selection.select_optimal`
solves.

## 4. Margins

For adjacent variants `j` and `j+1` the **phase relation** is
`d_j = h_j XOR h_(j+1)`. It is the only part of the phasing that survives the
complement symmetry, and it is what a switch error gets wrong. The **MEC
margin** of that decision is

```text
margin_j = ( min over h with d_j reversed of cost(h) ) - wMEC(R).
```

A margin of zero means the objective is indifferent to the decision. That
happens for a structural reason when no read spans the boundary, and it also
happens when the reads that do span it carry no net evidence, which is a
different situation with the same number attached to it.

Margins are read off the forward and backward tables of the same dynamic
program, at both settings of `d_j`, so all `m-1` of them cost one extra
backward pass rather than `m-1` constrained re-solves.
`phasebudget.exact.bruteforce_margins` performs the re-solves and is what the
block-wise extraction is tested against.

## 5. Competitive-programming formulation

### Input

```text
n m c
l_0 u_0 R_0[l_0] ... R_0[u_0] w_0[l_0] ... w_0[u_0]
...
```

`-1` represents an uncalled allele inside a span. Weights are non-negative
reals.

### Exact constraints for correctness testing

```text
1 <= n <= 20
1 <= m <= 20
0 <= c <= n
0 <= w <= 1000
```

Every solver in this repository answers these instances, and
`phasebudget.exact` answers them by enumeration, so a checker can compare
optimal cost exactly. Optimal solutions are not unique, so the checker compares
cost and verifies the reported assignment reproduces the reported cost. It does
not compare haplotypes.

### Performance division

The production track raises `m` to the number of heterozygous variants on a
human chromosome and `n` to sequencing depth times chromosome length divided by
read span. Submissions are measured on:

* switch error rate and flip error count against a truth set;
* phased fraction, block N50, AN50, and QAN50;
* wMEC cost of the reported haplotype under the full read set;
* selection objective as a fraction of the exact optimum;
* wall time and peak resident memory, reported separately for selection and
  for solving.

No composite score is defined. Accuracy and completeness trade against each
other, and a single number would hide which of the two a method bought.

## 6. Research questions

1. Over what range of coverage, read span, error rate, and heterozygosity does
   the wMEC optimum coincide with the true phasing, and where does it stop?
2. How large is the optimality gap of deployed read selection heuristics on the
   additive objective, and does closing that gap improve phasing accuracy?
3. Does the additive objective remain in P for redundancy-sensitive scores,
   such as one that pays for covering a boundary rather than for each read
   independently?
4. Do MEC margins identify the boundaries that carry switch errors, and does
   that signal decay as the objective stops tracking the truth?
5. Can reweighting by margin beat exact selection under the unweighted
   objective, and is the improvement stable or an artefact of the iteration cap?

Question 4 is where the first and third components of the project meet, and
question 5 is the deliverable that depends on both.
