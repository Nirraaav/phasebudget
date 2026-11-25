# What the complementary-haplotype restriction costs

## The restriction, and why it is not free

`phasebudget.dp` treats every column as a heterozygous site, so the second
haplotype is the bitwise complement of the first. `docs/design/problem-statement.md`
records this as a decision. Reading Patterson et al. 2015 changed its standing:
WhatsHap states the all-heterozygous assumption and then declines it, "Some
methods need to make this assumption for computational reasons, and this is why
we mention it, although our tool does not need it." So the restriction is a
divergence from the reference algorithm rather than a shared premise, and the
question of what it costs is OQ-7.

The case where it bites is a genotyping error: a site called heterozygous whose
two copies actually carry the same allele. Reads from both copies then agree
there. The complementary model cannot represent that column. It must claim the
copies differ and pay for the difference by mis-correcting reads, and the
mis-correction can propagate into the phase of neighbouring sites, which is the
part that matters.

## The two models

The general model is the classical MEC formulation, implemented as the oracle
`exact.solve_bipartition_bruteforce`: split the reads into two parts and
correct each part to its own weighted column consensus, with nothing tying the
two consensus values at a variant. A homozygous-looking column costs it
nothing, because both parts can agree.

The general model is a relaxation, so its optimum is a lower bound on the
complementary optimum. That inequality is already tested. What is not measured
is whether the extra freedom buys accuracy, which is a different question:
a relaxation can lower the objective and leave the answer no better.

## The common metric

The two models return different objects. The complementary model returns a
haplotype plus a read assignment; the general model returns only a read
bipartition. Comparing them needs a quantity both produce, and the read
bipartition is the one:

```text
bipartition error = min( hamming(a, origin), hamming(a, 1 - origin) ) / n_reads,
```

where `origin` is the copy each read was drawn from and the minimum removes the
global label flip, which carries no information in either model.

Phase accuracy is then derived the same way for both. From a bipartition, take
the weighted consensus of part 0 at every variant to get a predicted haplotype.
Phase is only defined between sites that really are heterozygous, so the
relations compared are those between consecutive members of
`heterozygous_sites()`, and a switch is a relation that disagrees with the
truth. Scoring a phase relation at a miscalled site would score a quantity that
does not exist.

## What each outcome would mean

- The general model has the lower bipartition error and the lower switch error.
  Then the restriction costs accuracy, the size of the cost is the answer to
  OQ-7, and the report has to state it as a limit of every number in the
  repository.
- The two tie on accuracy while the general model has the lower cost. Then the
  restriction is free in the regime tested, and the relaxation only buys
  objective value, which is the outcome the surrogate results elsewhere in this
  project would predict.
- The complementary model does better. Then the restriction is acting as a
  correct prior, which it partly is: the sites were called heterozygous for a
  reason and most of them are.

## The general-model dynamic program

`solve_bipartition_bruteforce` enumerates all `2^n` read bipartitions, so the
first comparison ran on ten to sixteen reads. `phasebudget.general` removes
that limit with the WhatsHap recurrence in its general form.

The state at variant `j` is a bipartition of the reads active at `j`, a mask
`S` whose set bits are the reads on side 1; there is no allele bit, because the
two sides choose their alleles independently. For a read `i` active at `j`,
`c0[i]` is the weight of changing its call to 0 and `c1[i]` the weight of
changing it to 1, both zero at an internal gap. The column cost is

```text
cost_j(S) = min(C0 - A0(S), C1 - A1(S)) + min(A0(S), A1(S)),
A0(S) = sum over i in S of c0[i],   A1(S) = sum over i in S of c1[i],
```

with `C0` and `C1` the sums over all active reads: side 0 is corrected to its
cheaper consensus and so is side 1. The recurrence is the complementary one
with this cost:

```text
F_j(S) = cost_j(S) + min over S' at j-1 agreeing with S on shared reads of F_(j-1)(S').
```

A bipartition and its complement cost the same, so the optimum is found twice
and the traceback returns one of them. The oracle is
`solve_bipartition_bruteforce`; `tests/test_general.py` requires equal costs
on random instances, including instances with genotype errors, and checks that
the returned assignment achieves the cost. It is `O(n 2^k)` for `n` variants
and `k` reads active at once, like the complementary pass, and is written in
Python only, which is enough for the coverages OQ-7 needs.
