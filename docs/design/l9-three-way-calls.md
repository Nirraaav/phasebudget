# L9: three-way end-variant calls at r = 1e-3 (development, rules stated before computing)

L2 certified within-gene cis or trans calls at r = 1e-2 when the assembly gives
the phaser's relation at the two end variants. At r = 1e-3 the 10 to 100 kb
bins certified nothing: the pool had 9 to 14 wrong picks of about 11,500 in
those bins, and certifying 1e-3 at delta 0.05 with about 11,500 picks needs
roughly 5 or fewer. L9 asks whether requiring a third source closes that gap.

L9 is a development analysis. NA12877 and HG005 were scored in L2, so their
numbers here are not a blind test. NA12878 is not read (it is L7's blind test).

## Claim

A call is wrong only if the reported relation differs from the truth at the
two end variants. Under the three-way rule a pair is called only when
whatshap, LongPhase and the assembly all relate it, each within one block or
contig, and all three give the same relation; the call is that relation. For
the call to be wrong, all three sources must be wrong in the same way at the
two ends. If the phasers' long-range errors are partly independent of each
other, the three-way rule has fewer wrong picks than either two-way rule.

## Computation

Everything is L2's (`benchmarks/l2_end_calls.py`) except the call rule: same
cached per-variant tables (`data/l2/`), pool (HG002 and the five Pedigree
genomes, 30x), tests (NA12877, HG005), clinical gene list, pairs at most 100 kb
apart, groups (distance bin 1, 10, 100 kb by easy or hard), one CRC32 pick per
(pool genome, gene, group), a group certified at r when
P(X <= k | n, r) <= 0.05, targets 1e-3 and 1e-2, test validity level
0.05/8 as in L2. The two L2 assembly arms (whatshap and LongPhase each with
the assembly) are recomputed beside the three-way arm for comparison.

Reported per arm and test genome, beside per-pair shares: the share of genes
with at least one within-100 kb pair that have at least one certified pair, and
the median over those genes of the per-gene certified share.

## Decision rules

- **L9-1, pool.** Under the three-way rule both 10 to 100 kb groups (easy and
  hard) certify at r = 1e-3.
- **L9-2, NA12877.** Under the three-way rule at r = 1e-3, at least 0.5 of
  NA12877's within-gene pairs 10 to 100 kb apart are certified.
- **L9-3, validity.** No three-way test cell (NA12877, HG005 by r) has
  P(X >= k | n, r) below 0.05/8 on its own picks.

L9 succeeds when all three hold. A success is a development result: it
earns a registered test on a genome not yet read, not a claim.

## Result (22 November 2026, after scoring)

`benchmarks/results/runs/2026-11-22_l9_three_way.txt`, 364 s on the server.
L9 does not succeed: L9-1 and L9-2 fail, L9-3 holds.

- Pool, 10 to 100 kb: the three-way rule has 6 wrong picks of 11,472 (`easy`)
  and 6 of 11,632 (`hard`), against 9 and 14 for whatshap with the assembly
  and 13 and 12 for LongPhase with the assembly. Certifying 1e-3 with about
  11,500 picks allows at most 5. Both groups miss by one pick.
- The three-way rule cuts wrong called pairs in that bin by about a quarter
  to a half (`easy` 20,322 whatshap and 20,721 LongPhase to 14,817; `hard`
  20,709 and 16,184 to 10,877) and costs under 1% of called pairs.
- NA12877 and HG005 at 1e-3: 0 certified pairs 10 to 100 kb apart under every
  arm. At 1e-2, three-way certifies 0.883 of NA12877's pairs with 88 wrong
  certified pairs, against 498 (whatshap) and 2,081 (LongPhase) under L2's
  arms; HG005 0.946 with 7 wrong, against 31 and 11.
- Per gene at 1e-2 (three-way): NA12877 0.821 of 2,672 genes have a certified
  pair, median per-gene share 0.929; HG005 0.968 of 2,498, median 1.000.

Reading: the third source helps, but at r = 1e-3 the limit is the number of
independent picks (one per pool genome, gene and group). Six pool genomes give
about 11,500 picks, too few to separate 6 wrong picks from 1e-3. More pool
genomes or a pick unit finer than the gene would be needed; neither was tried.
