"""Greedy (`whatshap`) against `optimal` selection on the extra subsample draws.

Run with the same `PHASEBUDGET_PLATFORM` and `PHASEBUDGET_SEEDS` as
`genome_run.py`. Reuses `score_genome.score_group` and its paired sign test and
window-bootstrap interval, restricted to the one comparison the extra draws were
made for. The registered seed-0 rows are in the original score files.
"""

from __future__ import annotations

import argparse
import multiprocessing
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402
import score_genome as sg  # noqa: E402

from phasebudget.paired import sign_test  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=12)
    parser.add_argument("--truth", choices=("filtered", "v5", "transmission"),
                        default="transmission" if common.PLATFORM == "hg005" else "v5")
    args = parser.parse_args()
    if not common.SEEDS:
        raise SystemExit("set PHASEBUDGET_SEEDS")
    sg.PRIMARY = args.truth
    sg.COMPARISONS = (("whatshap", "optimal", sg.CAPS),)
    sg.load_all_truths()
    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                            text=True, cwd=common.ROOT).stdout.strip()
    print(f"extra subsample draws, greedy (whatshap) against optimal; commit {commit}")
    print(f"platform {common.PLATFORM}; seeds {', '.join(map(str, common.SEEDS))}; truth "
          f"{args.truth}; {len(sg.CHROMS)} confirmatory chromosomes; diff = greedy minus "
          f"optimal per 10,000 common pairs, {sg.RESAMPLES}-resample 1 Mb window bootstrap")
    groups = [("lib", d, c) for d in common.DEPTHS for c in sg.CAPS]
    with multiprocessing.get_context("fork").Pool(args.jobs) as pool:
        results = pool.map(sg.score_group, groups, chunksize=1)
    print(f"{'depth':>6} {'seed':>4} {'cap':>3} {'pairs':>8} {'greedy':>6} {'optimal':>7} "
          f"{'p':>9} {'diff per 10k [95%]':>26}")
    for (_, depth, cap), _, paired, _, _ in results:
        q = paired[("whatshap", "optimal", args.truth)]
        p = sign_test(q["first"], q["first"] + q["second"])
        target, seed = common.depth_seed(depth)
        print(f"{target:>6g} {seed:>4} {cap:>3} {q['pairs']:>8} {q['first']:>6} "
              f"{q['second']:>7} {p:>9.2e} {sg.fmt(sg.difference_interval(q['windows'])):>26}")


if __name__ == "__main__":
    main()
