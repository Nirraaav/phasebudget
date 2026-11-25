"""What the span convention gives up: OPT_span(c) against OPT_called(c).

`docs/design/coverage-convention.md` fixes coverage as span coverage, which
counts a read over its internal gaps, and shows the called convention, which
counts only called sites, is NP-hard. Section 3 proves OPT_span <= OPT_called;
section 6 names the difference as unmeasured. This measures it: the certified
span optimum by min-cost flow, the called optimum by 0-1 program (CBC through
PuLP), their ratio, and the largest span coverage the called optimum would put
on one variant, which is the state count the dynamic program would face to
phase it.
"""

from __future__ import annotations

import argparse
import statistics

from phasebudget.gapped import select_called_ilp
from phasebudget.selection import read_weights, select_optimal
from phasebudget.simulate import simulate


def run(caps, gap_rates, replicates, n_variants, coverage, span):
    rows = []
    for gap_rate in gap_rates:
        for cap in caps:
            ratios, spans, holes = [], [], []
            for replicate in range(replicates):
                matrix = simulate(
                    n_variants, coverage=coverage, read_span=span, error_rate=0.05,
                    gap_rate=gap_rate, seed=6100 + 19 * replicate,
                ).matrix
                weights = read_weights(matrix)
                span_opt = select_optimal(matrix, cap, weights)
                called = select_called_ilp(matrix, cap, weights)
                ratios.append(span_opt.objective / called.objective if called.objective else 1.0)
                spans.append(matrix.subset(called.selected).max_coverage())
                holes.append(sum(a is None for f in matrix for a in f.alleles) /
                             max(1, sum(len(f) for f in matrix)))
            rows.append({
                "gap_rate": gap_rate,
                "cap": cap,
                "ratio": statistics.mean(ratios),
                "worst": min(ratios),
                "span_coverage": statistics.mean(spans),
                "span_coverage_max": max(spans),
                "holes": statistics.mean(holes),
            })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=int, default=60)
    parser.add_argument("--coverage", type=float, default=20.0)
    parser.add_argument("--span", type=int, default=6)
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--caps", type=int, nargs="+", default=[3, 5, 8])
    parser.add_argument("--gap-rates", type=float, nargs="+", default=[0.0, 0.1, 0.2, 0.3])
    args = parser.parse_args()
    print(f"{args.variants} variants, depth {args.coverage:g}, span {args.span}, error 0.05, "
          f"{args.replicates} replicates, seed 6100 + 19 * replicate; simulated data")
    print("ratio is OPT_span / OPT_called under read_weights; span cov is the largest span "
          "coverage of the called optimum, the state count the DP would need")
    header = (f"{'gaps':>5} {'cap':>4} {'hole share':>11} {'ratio':>7} {'worst':>7} "
              f"{'span cov':>9} {'max':>4}")
    print(header)
    print("-" * len(header))
    for row in run(args.caps, args.gap_rates, args.replicates, args.variants,
                   args.coverage, args.span):
        print(f"{row['gap_rate']:>5.2f} {row['cap']:>4} {row['holes']:>11.3f} "
              f"{row['ratio']:>7.4f} {row['worst']:>7.4f} {row['span_coverage']:>9.1f} "
              f"{row['span_coverage_max']:>4}")


if __name__ == "__main__":
    main()
