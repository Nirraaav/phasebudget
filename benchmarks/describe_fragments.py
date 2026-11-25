"""Describe an extracted fragment matrix: coverage, read content, qualities.

This is the OQ-9 check of `docs/design/real-data-plan.md`: whether the
registered caps fit the state limit and how much the matrix needs selection.
Coverage is span coverage, the convention every selector enforces, reported
over all reads and over reads calling two or more variants (the ones any
selector can choose).
"""

from __future__ import annotations

import argparse
import importlib.util
import statistics
from pathlib import Path

from phasebudget.dp import MAX_ACTIVE_READS


def _extractor():
    path = Path(__file__).with_name("extract_fragments.py")
    spec = importlib.util.spec_from_file_location("extract_fragments", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def quantiles(values, qs=(0.05, 0.25, 0.5, 0.75, 0.95, 1.0)):
    ordered = sorted(values)
    return [ordered[min(len(ordered) - 1, int(q * (len(ordered) - 1) + 0.5))] for q in qs]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fragments")
    parser.add_argument("--caps", type=int, nargs="+", default=[3, 5, 10, 15])
    args = parser.parse_args()
    matrix, positions, names, provenance = _extractor().load_fragments(args.fragments)
    informative = matrix.subset([f.read_id for f in matrix if f.n_observed >= 2])
    print("provenance: " + ", ".join(f"{k} {v}" for k, v in provenance.items()))
    print(f"variants {matrix.n_variants}; reads {len(matrix)}; calling two or more "
          f"{len(informative)}; span {positions[0] + 1} to {positions[-1] + 1}")
    labels = "5% 25% 50% 75% 95% max"
    for name, m in (("all reads", matrix), ("two or more", informative)):
        coverage = m.coverage()
        print(f"span coverage per variant, {name} ({labels}): "
              + " ".join(map(str, quantiles(coverage)))
              + f"; mean {statistics.mean(coverage):.1f}; variants with none "
              f"{sum(1 for c in coverage if c == 0)}")
    called = [f.n_observed for f in informative]
    spans = [len(f) for f in informative]
    holes = [len(f) - f.n_observed for f in informative]
    print(f"variants called per read ({labels}): " + " ".join(map(str, quantiles(called)))
          + f"; mean {statistics.mean(called):.1f}")
    print(f"variants spanned per read ({labels}): " + " ".join(map(str, quantiles(spans))))
    print(f"reads with a hole: {sum(1 for h in holes if h) / len(holes):.3f}; "
          f"holes per read, mean {statistics.mean(holes):.2f}")
    weights = [w for f in informative for _, _, w in f.observed()]
    print(f"allele quality ({labels}): " + " ".join(f"{w:g}" for w in quantiles(weights))
          + f"; zero {sum(1 for w in weights if w == 0)} of {len(weights)}")
    coverage = informative.coverage()
    for cap in args.caps:
        above = sum(1 for c in coverage if c > cap) / len(coverage)
        print(f"cap {cap}: fits the state limit of {MAX_ACTIVE_READS} "
              f"{'yes' if cap <= MAX_ACTIVE_READS else 'no'}; share of variants whose "
              f"coverage the cap binds: {above:.3f}")
    blocks = informative.blocks()
    print(f"blocks with every read: {len({b for b in blocks if b is not None})}; "
          f"variants in no block: {sum(1 for b in blocks if b is None)}")


if __name__ == "__main__":
    main()
