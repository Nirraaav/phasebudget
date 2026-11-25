"""Milestone 1: where does the wMEC optimum coincide with the truth?

One grid over depth, error rate, variant density and base-quality spread, on
the genomic simulator (`docs/design/genomic-simulator.md`), with every read
given to the solver. Per cell:

  * truth_optimal: the fraction of instances whose true haplotype is an
    optimum of wMEC, the event of Majidian et al.'s Theorem 2;
  * gap: the wMEC cost of the truth minus the optimum, the objective-truth gap;
  * switch: the switch error rate of the optimal phasing;
  * margin_info and posterior_info: the rank probability that a boundary
    carrying a switch has a smaller MEC margin, or posterior margin, than one
    that does not; 0.5 is no information (`docs/design/posterior-margins.md`).

Rows go to a CSV whose comment lines record every parameter, so the table and
the figures made from it carry their own provenance.

An instance whose active coverage somewhere exceeds the exact solver's state
limit cannot be given every read, so it is skipped and another seed drawn;
`skipped` counts them per cell. Skipping favours shallower instances, so a
cell with many skips is biased toward the truth being recovered.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import itertools
import math
import os
import platform
import statistics
import subprocess
from multiprocessing import Pool
from pathlib import Path

from phasebudget._validation import CoverageLimitExceeded
from phasebudget.dp import CoverageDP
from phasebudget.metrics import phasing_from_solution, switch_positions
from phasebudget.posterior import posterior_margins
from phasebudget.simulate import simulate_genomic_connected

READ_LENGTH = 10_000


def _rank_probability():
    path = Path(__file__).with_name("benchmark_objective_truth.py")
    spec = importlib.util.spec_from_file_location("objective_truth", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.rank_probability


def cell(arguments):
    coverage, error_rate, heterozygosity, spread, replicates, n_variants = arguments
    rank_probability = _rank_probability()
    optimal = skipped = 0
    gaps, switch_rates, hard, soft = [], [], [], []
    draw = 0
    accepted = 0
    while accepted < replicates and draw < 20 * replicates:
        instance = simulate_genomic_connected(
            n_variants,
            coverage=coverage,
            read_length=READ_LENGTH,
            heterozygosity=heterozygosity,
            error_rate=error_rate,
            quality_spread=spread,
            seed=101 * draw + 7,
        )
        draw += 1
        matrix = instance.matrix
        try:
            dp = CoverageDP(matrix)
        except CoverageLimitExceeded:
            skipped += 1
            continue
        accepted += 1
        solution = dp.solve()
        gap = matrix.cost(instance.haplotype) - solution.cost
        gaps.append(gap)
        optimal += gap <= 1e-9 * max(1.0, solution.cost)
        phasing = phasing_from_solution(matrix, solution)
        switches = set(switch_positions(instance.haplotype, phasing))
        margins = dp.margins()
        switch_rates.append(len(switches) / len(margins) if margins else 0.0)
        for signal, scores in ((margins, hard), (posterior_margins(matrix), soft)):
            bad = [m.margin for m in signal if m.boundary + 1 in switches]
            good = [m.margin for m in signal if m.boundary + 1 not in switches]
            probability = rank_probability(bad, good)
            if not math.isnan(probability):
                scores.append(probability)
    return {
        "coverage": coverage,
        "error_rate": error_rate,
        "heterozygosity": heterozygosity,
        "variants_per_read": heterozygosity * READ_LENGTH,
        "quality_spread": spread,
        "truth_optimal": optimal / accepted if accepted else float("nan"),
        "gap": statistics.mean(gaps),
        "gap_worst": max(gaps),
        "switch": statistics.mean(switch_rates),
        "margin_info": statistics.mean(hard) if hard else float("nan"),
        "posterior_info": statistics.mean(soft) if soft else float("nan"),
        "instances_with_switch": len(hard),
        "replicates": accepted,
        "skipped": skipped,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=int, default=40)
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--coverages", type=float, nargs="+", default=[4, 6, 8, 10, 12])
    parser.add_argument(
        "--error-rates", type=float, nargs="+",
        default=[0.01, 0.03, 0.05, 0.10, 0.15, 0.20, 0.30],
    )
    parser.add_argument(
        "--heterozygosities", type=float, nargs="+", default=[0.0005, 0.001, 0.002]
    )
    parser.add_argument("--quality-spreads", type=float, nargs="+", default=[0.0, 2.0])
    parser.add_argument("--processes", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    parser.add_argument("--out", type=Path, default=Path("benchmarks/results/phase_diagram.csv"))
    args = parser.parse_args()

    grid = list(itertools.product(
        args.coverages, args.error_rates, args.heterozygosities, args.quality_spreads
    ))
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        commit += " plus uncommitted changes" if dirty else ""
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    header = [
        f"benchmark_phase_diagram.py; commit {commit}; {platform.system()} "
        f"{platform.release()} {platform.machine()}; Python {platform.python_version()}",
        f"{args.variants} variants; read length {READ_LENGTH}; {args.replicates} replicates "
        "per cell; seeds 101 * draw + 7, skipping draws above the state limit; every read "
        "to the solver; simulated data",
        f"coverages {args.coverages}; error rates {args.error_rates}; heterozygosities "
        f"{args.heterozygosities}; quality spreads {args.quality_spreads}",
    ]
    with Pool(args.processes) as pool:
        rows = pool.map(
            cell,
            [(c, e, h, q, args.replicates, args.variants) for c, e, h, q in grid],
            chunksize=1,
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as handle:
        for line in header:
            handle.write(f"# {line}\n")
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for line in header:
        print(line)
    print(f"{len(rows)} cells written to {args.out}")


if __name__ == "__main__":
    main()
