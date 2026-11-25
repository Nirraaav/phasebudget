"""The D1 pilot: error-optimal selection on the genome-wide HiFi fragments.

`docs/design/error-optimal-selection.md` fixes the protocol: three scales of
the switch-model cost at caps 3 and 5, tuned on the odd confirmatory
chromosomes and tested on the even ones against `optimal` by the paired
sign test. The run files of `optimal`, `connected` and `whatshap` come from
the registered genome-wide run; the new ones are written beside them by
`hg002_pipelines.py`, one process per job, and skipped when present.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/benchmark_error_optimal.py`
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402
import score_genome as sg  # noqa: E402

from phasebudget.metrics import block_lengths, comparable_pairs, n50, switch_positions  # noqa: E402
from phasebudget.paired import common_pairs, sign_test  # noqa: E402

SCALES = ("1e2", "1e4", "1e6")
CAPS = (3, 5)
BASELINES = ("optimal", "connected", "whatshap")


def run_file(chrom, depth, pipeline, cap):
    return common.chrom_dir(chrom) / "runs" / f"lib_{depth}_{pipeline}_c{cap}.json.gz"


def run_all(jobs: int) -> None:
    todo = [(c, d, f"error_optimal_{s}", cap) for c in common.CONFIRMATORY for d in common.DEPTHS
            for s in SCALES for cap in CAPS if not run_file(c, d, f"error_optimal_{s}", cap).exists()]
    print(f"{len(todo)} pilot jobs to run", flush=True)

    def one(job):
        chrom, depth, pipeline, cap = job
        subprocess.run(["nice", "-n", "19", sys.executable, str(HERE / "hg002_pipelines.py"),
                        "--fragments", str(common.fragments_path(chrom, depth, "realign")),
                        "--pipeline", pipeline, "--cap", str(cap), "--weights", "logodds",
                        "--out", str(run_file(chrom, depth, pipeline, cap))],
                       check=True, capture_output=True)

    with ThreadPoolExecutor(jobs) as pool:
        list(pool.map(one, todo))


def score(chroms):
    """Pooled switches, pairs, N50 per pipeline, and discordant pairs against `optimal`."""
    names = BASELINES + tuple(f"error_optimal_{s}" for s in SCALES)
    pooled, paired = {}, {}
    for chrom in chroms:
        truth, positions = sg.TRUTH[chrom]["filtered"], sg.POSITIONS[chrom]
        for depth in common.DEPTHS:
            for cap in CAPS:
                phasings = {n: sg.load_phasing("lib", chrom, depth, n, cap)[0] for n in names}
                for n, phasing in phasings.items():
                    p = pooled.setdefault((depth, cap, n), {"switches": 0, "pairs": 0,
                                                            "blocks": []})
                    p["switches"] += len(switch_positions(truth, phasing))
                    p["pairs"] += comparable_pairs(phasing, truth)
                    p["blocks"].extend(block_lengths(phasing, positions))
                    if n == "optimal":
                        continue
                    rows = common_pairs(truth, phasings[n], phasings["optimal"])
                    q = paired.setdefault((depth, cap, n), [0, 0, 0])
                    q[0] += len(rows)
                    q[1] += sum(1 for _, a, b in rows if a and not b)
                    q[2] += sum(1 for _, a, b in rows if b and not a)
    return pooled, paired


def report(label, pooled, paired):
    print(f"\n{label}: pooled against the filtered truth; 'vs optimal' counts common pairs "
          "where only this pipeline / only optimal is wrong")
    print(f"{'depth':>5} {'cap':>3} {'pipeline':>18} {'switch':>9} {'sw/pairs':>15} "
          f"{'N50 kb':>7} {'vs optimal':>11} {'p':>9}")
    for depth in common.DEPTHS:
        for cap in CAPS:
            for (d, c, n), p in pooled.items():
                if (d, c) != (depth, cap):
                    continue
                rate = p["switches"] / p["pairs"] if p["pairs"] else 0.0
                q = paired.get((d, c, n))
                versus = f"{q[1]:>5}/{q[2]:<5}" if q else ""
                pv = f"{sign_test(q[1], q[1] + q[2]):9.2e}" if q else ""
                print(f"{depth:>5} {cap:>3} {n:>18} {rate:9.6f} {p['switches']:>6}/{p['pairs']:<8} "
                      f"{n50(p['blocks']) / 1000:>7.0f} {versus:>11} {pv:>9}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=6)
    args = parser.parse_args()
    run_all(args.jobs)
    sg.load_all_truths()
    print("D1 pilot (docs/design/error-optimal-selection.md): switch-model costs, allele error "
          f"0.00117, scales {', '.join(SCALES)}, caps {CAPS}, depths {', '.join(common.DEPTHS)}; "
          "realign fragments, log-odds weights; CoverageDP phasing")

    pooled, paired = score(common.ODD)
    report("tuning, odd chromosomes " + ", ".join(common.ODD), pooled, paired)
    rates = {}
    for s in SCALES:
        keys = [k for k in pooled if k[2] == f"error_optimal_{s}"]
        rates[s] = sum(pooled[k]["switches"] for k in keys) / sum(pooled[k]["pairs"] for k in keys)
    best = min(SCALES, key=lambda s: (rates[s], float(s)))
    print("tuning: pooled switch error by scale, caps 3 and 5, all depths: "
          + ", ".join(f"{s} {rates[s]:.7f}" for s in SCALES) + f" -> {best}")

    pooled, paired = score(common.EVEN)
    report("test, even chromosomes " + ", ".join(common.EVEN), pooled, paired)
    name = f"error_optimal_{best}"
    counts = [paired[(d, 3, name)] for d in common.DEPTHS]
    ours, theirs = sum(c[1] for c in counts), sum(c[2] for c in counts)
    every = all(c[1] < c[2] for c in counts)
    total = sign_test(ours, ours + theirs)
    print(f"\ngate, cap 3, {name} vs optimal on the even chromosomes: discordant "
          + ", ".join(f"{d}x {c[1]}/{c[2]}" for d, c in zip(common.DEPTHS, counts))
          + f"; total {ours}/{theirs}, p {total:.2e}; fewer at every depth: {every} -> "
          + ("D1 goes on" if every and total < 0.05 and ours < theirs else "gate not met"))


if __name__ == "__main__":
    main()
