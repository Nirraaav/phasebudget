"""Follow-up gate G-D (`docs/followup/error-controlled-phasing.md`): can a guarantee certify small rates?

Units are 10 kb bins. A bin's loss is 1 if any trusted pair in it is wrong.
For a threshold on the Wilson score (`gate_evidence_score.SCORES`), Learn
then Test scans thresholds from strict (half the pairs trusted) to loose. At each one it tests
H0: P(bin loss = 1) > beta with an exact binomial p-value over calibration
bins, and stops at the first p-value above delta. The last accepted
threshold is used on the test bins, where the realized bin-error rate, the
pair-level error rate and the share of pairs trusted are reported.

Calibration uses the odd confirmatory chromosomes and the test uses the even
ones; HG002 HiFi `whatshap phase` at cap 15, full depth. These data were
seen, so this checks feasibility only. Bins are assumed independent, and
that assumption is the thing to examine next. Usage:
`PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/gate_d_feasibility.py`
"""

from __future__ import annotations

import bisect
import math
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402
from gate_evidence_score import SCORES  # noqa: E402
from gate_truth_check import TOOLS, existing  # noqa: E402

extractor = common._module("extract_fragments.py")
score = common._module("score_hg002.py")
BIN = int(os.environ.get("PHASEBUDGET_BIN", "10000"))
DELTA = 0.05
# `PHASEBUDGET_DEPTH` selects the ONT depth; HiFi uses the full-depth tool-level run.
import os  # noqa: E402

DEPTH = os.environ.get("PHASEBUDGET_DEPTH", "full")


def binom_tail_le(k: int, n: int, p: float) -> float:
    """P(X <= k) for X ~ Binomial(n, p), summed in log space."""
    total = 0.0
    for i in range(k + 1):
        log_term = (math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
                    + i * math.log(p) + (n - i) * math.log1p(-p))
        total += math.exp(log_term)
    return min(1.0, total)


def pairs_with_positions():
    """`[(chrom, position, score key, wrong)]` for every compared pair."""
    from gate_evidence_score import pair_scores

    key = SCORES["wilson"]
    out = []
    for chrom in common.CONFIRMATORY:
        variants = score.variant_table(common.hetsnv_path(chrom))
        v5 = common.v5_truth(chrom, variants)
        if common.PLATFORM == "hifi":
            vcf = TOOLS / chrom / "stock_c15.vcf"
        else:
            # The genome run's `whatshap phase` at default settings, at this depth.
            vcf = common.chrom_dir(chrom) / "external" / f"whatshap_{DEPTH}.vcf"
        wh = score.whatshap_phasing(existing(vcf), variants)
        matrix, _, _, _ = extractor.load_fragments(
            str(common.fragments_path(chrom, DEPTH, "realign")))
        rows = pair_scores(matrix.fragments, wh.alleles, wh.blocks, v5)
        phased = [j for j, a in enumerate(wh.alleles) if a is not None and v5[j] is not None]
        starts = [a for a, b in zip(phased, phased[1:]) if wh.blocks[a] == wh.blocks[b]]
        for j, r in zip(starts, rows):
            out.append((chrom, variants[j][0], key(r), r[2]))
    return out


class Bins:
    """Per 10 kb bin, the best score of any pair and of any wrong pair, sorted for counting."""

    def __init__(self, rows):
        best, worst = {}, {}
        for chrom, pos, k, wrong in rows:
            b = (chrom, pos // BIN)
            best[b] = max(best.get(b, k), k)
            if wrong:
                worst[b] = max(worst.get(b, k), k)
        self.all = sorted(best.values())
        self.wrong = sorted(worst.values())

    def at(self, threshold):
        """`(bins with a trusted pair, bins with a trusted wrong pair)` at `threshold`."""
        return (len(self.all) - bisect.bisect_left(self.all, threshold),
                len(self.wrong) - bisect.bisect_left(self.wrong, threshold))


def main() -> None:
    rows = pairs_with_positions()
    odd = {f"chr{i}" for i in range(1, 23, 2)}
    cal = [r for r in rows if r[0] in odd]
    test = [r for r in rows if r[0] not in odd]
    keys = sorted((r[2] for r in cal), reverse=True)
    cal_bins, test_bins = Bins(cal), Bins(test)
    # Thresholds from strict to loose: the calibration score at which 50, 50.25, ..., 100 per cent
    # of pairs are trusted. A strict start with few bins could never reject, so it starts at half.
    grid = [keys[min(len(keys) - 1, int(c / 400 * len(keys)) - 1)] for c in range(200, 401)]
    print(f"gate G-D, HG002 {common.PLATFORM} whatshap phase, depth {DEPTH}; {BIN // 1000} kb bins; calibration "
          f"{len(cal)} pairs (odd chromosomes), test {len(test)} pairs (even); delta {DELTA}")
    print(f"{'beta':>8} | {'threshold':>22} | {'test bins wrong':>15} | {'test pair errors':>16} "
          f"| {'pairs trusted':>13}")
    for beta in (1e-4, 3e-4, 1e-3, 3e-3, 1e-2):
        chosen = None
        for t in grid:
            n, k = cal_bins.at(t)
            if binom_tail_le(k, n, beta) > DELTA:
                break
            chosen = t
        if chosen is None:
            print(f"{beta:>8.0e} | {'none certified':>22} |")
            continue
        tn, tk = test_bins.at(chosen)
        trusted = [r for r in test if r[2] >= chosen]
        wrong = sum(r[3] for r in trusted)
        print(f"{beta:>8.0e} | {str(tuple(round(x, 4) for x in chosen)):>22} | "
              f"{tk:>6}/{tn:<8} | {wrong:>5}/{len(trusted):<10} | "
              f"{len(trusted) / len(test):>13.4f}")


if __name__ == "__main__":
    main()
