# Two-level exact phasing

## What it optimizes

The coverage-parameterized DP phases a selected subset of the reads exactly.
Where the subset leaves a boundary unspanned, the subset's optimum is
defined separately on each side: its phase blocks `B_1, ..., B_m` each carry
a haplotype that is determined only up to complement. The relation between
blocks is then arbitrary, and the reads the selection discarded are never
consulted about it. OQ-13 found that the adaptive loop's gain comes from
consulting all reads, and OQ-15 found that the stitcher consulted them in
the wrong places.

Two-level phasing keeps the subset's haplotype inside every block and
chooses one orientation `o_i` in {0, 1} per block, by minimizing the
full-matrix wMEC:

```text
minimize over o_2, ..., o_m in {0,1}   sum over reads r in the full matrix of
    min over s in {0,1} of  sum over blocks i the read calls alleles in of
        cost_r,i(o_i XOR s)
```

where `cost_r,i(0)` is the weight of read `r`'s alleles in block `i` that
disagree with the block's haplotype, `cost_r,i(1)` the weight of those that
agree, and `o_1 = 0` fixes the global complement. Only alleles at variants
inside a block are counted. A variant no selected read covers keeps an
arbitrary allele in the subset's solution and is reported unphased, so
scoring it would let that arbitrary allele pull on the orientation. The
first version counted it, and the exhaustive test caught the orientation
raising the cost. The returned haplotype is
the subset optimum with every block complemented where `o_i = 1`. Variants
no selected read covers are left unphased, as before.

This is the full-matrix wMEC, over the phased variants, restricted to
haplotypes that agree with the subset optimum inside each block. Because
the identity orientation is a candidate, that cost never rises. It loses nothing when every
block haplotype is the full-matrix optimum restricted to that block; it is
not in general the full-matrix optimum.

## How it is solved exactly

A read that calls alleles in one block contributes a constant. A read that
calls alleles in blocks `a < ... < b` is a factor over the orientations of
the blocks from `a` to `b`. The blocks form a chain, and every factor lies
in a window of consecutive blocks, so a dynamic program over blocks, whose
state is the orientation of the last `W - 1` blocks, is exact, where `W` is
the widest window of any factor. It costs `O(m 2^W)` plus the factor
evaluations. When `W` exceeds a limit, the widest reads are dropped and the
number dropped is reported, so the result is exact for the reads kept.
That case has to be visible, not silent.

## Blocks of the result

After orientation, two blocks belong to the same output block when some read
of the full matrix calls alleles in both, since that read decided their
relation. Output blocks are therefore unions of subset blocks, joined along
the chain of factors.

## Oracle and tests

`tests/test_two_level.py`:
- enumerates all `2^(m-1)` orientations on small instances and checks that
  the DP attains the minimum full-matrix cost;
- checks that the cost never exceeds the one-level cost;
- reproduces the HG002 failure of stitched acceptance: a block break inside
  a window, with crossing reads that favour one orientation.

## Relation to stitched acceptance

The stitched rule orients at fixed window junctions, greedily, using the
reads that cross each junction, and can take a round's haplotype per window.
Two-level phasing orients at exactly the places where orientation is free,
the block breaks, jointly and exactly, from one round. It has no window
parameter and no reweighting.
