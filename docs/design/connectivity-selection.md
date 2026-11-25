# Coverage-constrained selection that keeps blocks connected

## What it optimizes

The additive selector maximizes the summed weight of the kept reads subject
to the coverage cap `c` at every variant. It pays nothing for connectivity,
and at small caps it leaves boundaries unspanned (M2: block N50 halves at
cap 3). The connected selector adds a requirement on boundaries. Boundary
`b` lies between variants `b` and `b + 1`, and a read spans it when its span
contains both. For every boundary that at least one read of the matrix
spans:

```text
maximize    sum_r w_r x_r  -  lambda * sum_b s_b
subject to  sum_{r covers j} x_r        <= c        for every variant j
            sum_{r spans b} x_r + s_b   >= k        for every spannable boundary b
            x_r in {0, 1},  s_b >= 0
```

`k` asks for `k` spanning reads, and `s_b`, priced at `lambda` per missing
read, makes the requirement soft, so the problem is always feasible. With
`lambda` above the total read weight the requirement is hard wherever it can
be met. `k` is at most `c`.

## Why it is exact

The Carlisle and Lloyd network of `selection-integrality.md` has one node
between consecutive variants and one arc of capacity `c` per variant; `c`
units flow from the first node to the last, and each read is an arc of
capacity 1 and cost `-w_r` that jumps over the variants it covers. A flow
decomposes into `c` chains, each a set of reads with disjoint spans, so the
coverage of variant `j` is `c` minus the flow on its path arc.

Here every node between variants `b` and `b + 1` is split in two, joined by
a boundary arc. A read covering variants `start..end-1` jumps from before
`start` to after `end - 1`, so it passes over the boundary arcs of the
boundaries it spans and over no others. Two consecutive reads of one chain,
one ending at variant `b` and the next starting at `b + 1`, force their unit
through boundary `b`'s arc, which is right: neither spans it. The number of
selected reads spanning boundary `b` is therefore `c` minus the flow on its
boundary arcs. Requiring at least `k` of them is a capacity of `c - k` on
the boundary arc. The soft version adds a parallel arc of capacity `k` and
cost `lambda`; the flow on it is exactly `s_b`.

The network has integer capacities, so min-cost flow returns an integral
optimum, which is the optimum of the program above. This is the same
argument as for the additive selector, applied to the interleaved rows
`v_1, b_1, v_2, ..., v_n`, in whose order every read's rows are consecutive.

## Novelty

Mäkinen, Staneva, Tomescu and Valenzuela 2015 solve the max-min coverage
version by flow; the combination of packing rows (variants) and covering
rows (boundaries) on an interval matrix may be known in the scheduling
literature. Until a search has been recorded in the prior-art matrix, this
is claimed as an application, not as a new algorithm.

## Oracle and tests

`tests/test_connected_selection.py` enumerates every subset of up to 12
reads, keeps those within the cap, scores them by the objective above, and
checks that the flow reaches the maximum. It also checks that `lambda = 0`
returns the additive optimum's objective, that with a hard `lambda` every
spannable boundary the cap allows is spanned, and that the cap holds.
