"""Run `whatshap phase` with its read selection replaced by the exact optimum.

Everything else is WhatsHap 2.8 as shipped: read extraction, realignment,
genotype handling, the DP, block construction and VCF output. Only
`whatshap.cli.phase.readselection` is swapped for a selector that maximizes
the total allele quality of the selected reads under the same coverage cap,
by `phasebudget.selection.select_optimal`. Reads calling fewer than two
variants are never selected, as in WhatsHap.

Usage: `PYTHONPATH=src .venv/bin/python benchmarks/whatshap_exact.py phase [whatshap phase args]`
"""

from __future__ import annotations

import sys

from phasebudget._types import Fragment
from phasebudget.fragments import FragmentMatrix
from phasebudget.selection import read_weights, select_optimal


def exact_readselection(readset, max_coverage, preferred_source_ids=None):
    """Indices into `readset` of the maximum-quality selection under `max_coverage`."""
    positions = sorted({variant.position for read in readset for variant in read})
    index = {p: j for j, p in enumerate(positions)}
    fragments = []
    for k, read in enumerate(readset):
        calls = sorted((index[v.position], v.allele, v.quality) for v in read)
        if len(calls) < 2:
            continue
        start, last = calls[0][0], calls[-1][0]
        alleles = [None] * (last - start + 1)
        weights = [0.0] * (last - start + 1)
        for j, allele, quality in calls:
            alleles[j - start] = allele
            weights[j - start] = float(quality)
        fragments.append(Fragment(k, start, tuple(alleles), tuple(weights)))
    if not fragments:
        return set()
    matrix = FragmentMatrix(fragments, len(positions))
    return set(select_optimal(matrix, max_coverage, read_weights(matrix)).selected)


def main() -> None:
    import whatshap.cli.phase as phase_module
    from whatshap.__main__ import main as whatshap_main

    phase_module.readselection = exact_readselection
    sys.argv = ["whatshap", *sys.argv[1:]]
    whatshap_main()


if __name__ == "__main__":
    main()
