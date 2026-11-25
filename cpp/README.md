# C++

`dp.cpp` is the coverage-parameterized dynamic program of `phasebudget.dp` in
C++. It moved because the Python pass could not reach the caps the deployed
tools use: on 60 variants, a solve with all margins took 206 ms at cap 8 and
3.3 s at cap 12 in Python, against 7 ms and 41 ms here, and 188 ms at cap 15.

Build it with

```
make -C cpp
```

which writes `src/phasebudget/_pbdp.so`. `CoverageDP(backend="auto")`, the
default, uses it when it exists; `PHASEBUDGET_NATIVE=0` turns it off, and the
library has no runtime dependency either way. It is loaded through `ctypes`.

## What it answers to

The native pass repeats the Python arithmetic in the same order and breaks
ties the same way, so the two return equal floats, not close ones.
`tests/test_native.py` checks cost, haplotype, read assignment and every
boundary cost for equality on 120 random instances, at every checkpoint
spacing, and checks costs against the exhaustive oracle in `phasebudget.exact`.
Two details make the equality hold: the build uses `-ffp-contract=off`, so no
multiply and add are fused, and the column totals use
`phasebudget.dp.column_sum`, which is CPython 3.12's compensated `sum`.

Memory is `O(sqrt(n) 2^k)` for `n` variants and cap `k`: the forward pass
keeps one table in every `sqrt(n)`, and a single descending sweep recomputes
each segment to run the traceback, the backward pass and the boundary costs.

The min-cost flow stays in Python. It takes about a millisecond per instance
and is not the bottleneck.
