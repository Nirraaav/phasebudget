"""Does Fischer and Marschall's published selection phase better than WhatsHap's?

WhatsHap 2.8 ships a selector whose score update departs from the published
procedure (OQ-14). Its extra errors at cap 3 sit at boundaries it leaves with
a single spanning read, and its selection is already maximal, so no read can
be added afterwards. This runs the published procedure (`fm_paper`) through
the same phasing as every library pipeline and compares it with `whatshap`
and `optimal` on common pairs against GIAB v5.0q.

On HiFi this is exploratory, since the genome-wide results were seen. The
platform follows `PHASEBUDGET_PLATFORM`. Usage:
`PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/benchmark_selector_fix.py`
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

NAMES = ("whatshap", "fm_paper", "optimal")


def run_file(chrom, depth, pipeline, cap):
    return common.chrom_dir(chrom) / "runs" / f"lib_{depth}_{pipeline}_c{cap}.json.gz"


def run_all(chroms, caps, jobs):
    todo = [(c, d, cap) for c in chroms for d in common.DEPTHS for cap in caps
            if not run_file(c, d, "fm_paper", cap).exists()]
    print(f"{len(todo)} fm_paper runs to do", flush=True)

    def one(job):
        chrom, depth, cap = job
        subprocess.run(["nice", "-n", "19", sys.executable, str(HERE / "hg002_pipelines.py"),
                        "--fragments", str(common.fragments_path(chrom, depth, "realign")),
                        "--pipeline", "fm_paper", "--cap", str(cap), "--weights", "logodds",
                        "--out", str(run_file(chrom, depth, "fm_paper", cap))],
                       check=True, capture_output=True)

    with ThreadPoolExecutor(jobs) as pool:
        list(pool.map(one, todo))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chromosomes", nargs="+", default=list(common.CONFIRMATORY))
    parser.add_argument("--caps", nargs="+", type=int, default=[3, 5])
    parser.add_argument("--jobs", type=int, default=3)
    args = parser.parse_args()
    run_all(args.chromosomes, args.caps, args.jobs)
    sg.CHROMS = tuple(args.chromosomes)
    sg.load_all_truths()
    print(f"selector fix, platform {common.PLATFORM}, truth GIAB v5.0q, "
          f"{len(args.chromosomes)} chromosomes; paired columns count common pairs where only "
          "the first / only the second is wrong")
    print(f"{'depth':>5} {'cap':>3} {'pipeline':>9} {'switches':>9} {'pairs':>9} {'N50 kb':>7} | "
          f"{'vs whatshap':>12} {'p':>9} | {'vs optimal':>12} {'p':>9}")
    for depth in common.DEPTHS:
        for cap in args.caps:
            pooled = {n: [0, 0, []] for n in NAMES}
            paired = {("fm_paper", "whatshap"): [0, 0], ("fm_paper", "optimal"): [0, 0]}
            for chrom in args.chromosomes:
                truth, positions = sg.TRUTH[chrom]["v5"], sg.POSITIONS[chrom]
                ph = {n: sg.load_phasing("lib", chrom, depth, n, cap)[0] for n in NAMES}
                for n in NAMES:
                    pooled[n][0] += len(switch_positions(truth, ph[n]))
                    pooled[n][1] += comparable_pairs(ph[n], truth)
                    pooled[n][2].extend(block_lengths(ph[n], positions))
                for (a, b), acc in paired.items():
                    rows = common_pairs(truth, ph[a], ph[b])
                    acc[0] += sum(1 for _, x, y in rows if x and not y)
                    acc[1] += sum(1 for _, x, y in rows if y and not x)
            for n in NAMES:
                s, pairs, blocks = pooled[n]
                extra = ""
                if n == "fm_paper":
                    a = paired[("fm_paper", "whatshap")]
                    b = paired[("fm_paper", "optimal")]
                    extra = (f"| {a[0]:>5}/{a[1]:<6} {sign_test(a[0], a[0] + a[1]):>9.2e} | "
                             f"{b[0]:>5}/{b[1]:<6} {sign_test(b[0], b[0] + b[1]):>9.2e}")
                print(f"{depth:>5} {cap:>3} {n:>9} {s:>9} {pairs:>9} {n50(blocks) / 1e3:>7.0f} {extra}")
            print()


if __name__ == "__main__":
    main()
