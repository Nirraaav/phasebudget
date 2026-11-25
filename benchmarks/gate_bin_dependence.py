"""Follow-up: do wrong 10 kb bins cluster? A check of the independence G-D assumes.

For HG002 HiFi `whatshap phase` at cap 15 (21 chromosomes), it lists the
bins holding a wrong pair and counts pairs of wrong bins within 1, 2, 5, 10
and 100 bins of each other. It compares each count with its expectation if
the same number of wrong bins were spread uniformly over the bins that hold
compared pairs on each chromosome. Development data. Usage:
`PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/gate_bin_dependence.py`
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gate_d_feasibility import BIN, pairs_with_positions  # noqa: E402

DISTANCES = (1, 2, 5, 10, 100)


def close_pairs(wrong_bins, d):
    count = 0
    for chrom, bins in wrong_bins.items():
        s = sorted(bins)
        for i, a in enumerate(s):
            for b in s[i + 1:]:
                if b - a > d:
                    break
                count += 1
    return count


def main() -> None:
    rows = pairs_with_positions()
    occupied, wrong = {}, {}
    for chrom, pos, _, w in rows:
        occupied.setdefault(chrom, set()).add(pos // BIN)
        if w:
            wrong.setdefault(chrom, set()).add(pos // BIN)
    n_wrong = sum(len(v) for v in wrong.values())
    observed = {d: close_pairs(wrong, d) for d in DISTANCES}
    rng = random.Random(0)
    sims = {d: [] for d in DISTANCES}
    for _ in range(2000):
        fake = {c: set(rng.sample(sorted(occupied[c]), len(wrong[c]))) for c in wrong}
        for d in DISTANCES:
            sims[d].append(close_pairs(fake, d))
    from gate_d_feasibility import DEPTH
    import genome_common as common

    print(f"wrong-bin clustering, HG002 {common.PLATFORM} whatshap phase, depth {DEPTH}, {BIN // 1000} kb bins: "
          f"{n_wrong} wrong bins of {sum(len(v) for v in occupied.values())}; 2,000 uniform "
          "placements per chromosome, seed 0")
    for d in DISTANCES:
        mean = sum(sims[d]) / len(sims[d])
        p = (1 + sum(x >= observed[d] for x in sims[d])) / (1 + len(sims[d]))
        print(f"  wrong-bin pairs within {d:>3} bins: observed {observed[d]}, expected "
              f"{mean:.2f}, p {p:.4f}")


if __name__ == "__main__":
    main()
