"""Fit and test the posterior correction of `docs/design/calibration.md`.

`--fit` chromosomes give `T` and `eps` by maximum likelihood; `--test`
chromosomes are then scored with and without the correction. With `--T`
and `--eps` the fit is skipped. The platform follows `PHASEBUDGET_PLATFORM`.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/calibrate_posterior.py
--fit chr1 chr3 --test chr2 chr4`
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402

from phasebudget.bootstrap import statistic_interval  # noqa: E402

score = common._module("score_hg002.py")
WINDOW = 1_000_000


def logistic(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x)) if x >= 0 else math.exp(x) / (1.0 + math.exp(x))


def logit(p: float) -> float:
    p = min(max(p, 1e-300), 1 - 1e-16)
    return math.log(p) - math.log1p(-p)


def load(chrom: str, depth: str):
    """`[(window, p, wrong)]` for boundaries v5.0q phases on both sides."""
    path = common.chrom_dir(chrom) / "runs" / f"lib_{depth}_optimal_c5.json.gz"
    with gzip.open(path, "rt") as handle:
        signals = json.load(handle)["signals"]
    variants = score.variant_table(common.hetsnv_path(chrom))
    truth = common.v5_truth(chrom, variants)
    first = variants[0][0]
    rows = []
    for b, called, predicted, _, _ in signals:
        if truth[b] is None or truth[b + 1] is None:
            continue
        rows.append(((variants[b + 1][0] - first) // WINDOW, predicted,
                     called != (truth[b] ^ truth[b + 1])))
    return rows


def corrected(p: float, T: float, eps: float) -> float:
    return eps + (1 - eps) * logistic(T * logit(p))


def fit(rows):
    from scipy.optimize import minimize

    def loss(theta):
        T, eps = math.exp(theta[0]), logistic(theta[1])
        total = 0.0
        for _, p, wrong in rows:
            q = min(max(corrected(p, T, eps), 1e-300), 1 - 1e-16)
            total -= math.log(q) if wrong else math.log1p(-q)
        return total

    result = minimize(loss, [0.0, -12.0], method="Nelder-Mead",
                      options={"xatol": 1e-4, "fatol": 1e-4, "maxiter": 400})
    return math.exp(result.x[0]), logistic(result.x[1])


def evaluate(label, groups, T, eps):
    """`groups` is one row list per chromosome; windows are offset so they stay apart."""
    rows, offset = [], 0
    for group in groups:
        rows.extend((offset + w, p, wrong) for w, p, wrong in group)
        offset += max((w for w, _, _ in group), default=-1) + 1
    for name, (t, e) in (("uncorrected", (1.0, 0.0)), ("corrected", (T, eps))):
        probs = [(w, corrected(p, t, e), wrong) for w, p, wrong in rows]
        per_window = [[0.0, 0.0] for _ in range(offset)]
        for w, q, wrong in probs:
            per_window[w][0] += wrong
            per_window[w][1] += q

        def ratio(picks):
            predicted = sum(per_window[k][1] for k in picks)
            return None if predicted == 0 else sum(per_window[k][0] for k in picks) / predicted

        interval = statistic_interval(ratio, offset, resamples=2000, seed=0)
        bins = [[0, 0.0, 0.0] for _ in range(10)]
        for _, q, wrong in probs:
            k = min(int(q * 10), 9)
            bins[k][0] += 1
            bins[k][1] += wrong
            bins[k][2] += q
        ece = sum(abs(b[1] - b[2]) for b in bins) / len(probs)
        logloss = -sum(math.log(max(q, 1e-300)) if wrong else math.log1p(-min(q, 1 - 1e-16))
                       for _, q, wrong in probs) / len(probs)
        calibrated = interval.low <= 1 <= interval.high and ece < 0.01
        print(f"{label} {name:>11}: {len(probs)} boundaries, observed "
              f"{sum(w for _, _, w in probs)}, predicted {sum(q for _, q, _ in probs):.1f}, "
              f"ratio {interval.estimate:.2f} [{interval.low:.2f}, {interval.high:.2f}], ECE "
              f"{ece:.5f}, log-loss {logloss:.6f} -> "
              f"{'calibrated' if calibrated else 'not calibrated'}")
        decades = {}
        for _, q, wrong in probs:
            d = max(-8, math.floor(math.log10(max(q, 1e-300))))
            acc = decades.setdefault(d, [0, 0, 0.0])
            acc[0] += 1
            acc[1] += wrong
            acc[2] += q
        for d in sorted(decades):
            n, w, s = decades[d]
            print(f"    p' in [1e{d}, 1e{d + 1}){' or below' if d == -8 else ''}: {n:>8} "
                  f"boundaries, observed {w:>4}, predicted {s:9.2f}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fit", nargs="*", default=[])
    parser.add_argument("--test", nargs="*", default=[])
    parser.add_argument("--depths", nargs="+", default=["full"])
    parser.add_argument("--T", type=float)
    parser.add_argument("--eps", type=float)
    args = parser.parse_args()
    print(f"posterior calibration, platform {common.PLATFORM}, truth GIAB v5.0q, "
          "cap 5 optimal, realign fragments with log-odds weights")
    if args.T is None:
        rows = [r for c in args.fit for r in load(c, "full")]
        T, eps = fit(rows)
        print(f"fitted on {', '.join(args.fit)} at full depth ({len(rows)} boundaries, "
              f"{sum(w for _, _, w in rows)} wrong): T = {T:.4f}, eps = {eps:.3e}")
    else:
        T, eps = args.T, args.eps
        print(f"given: T = {T}, eps = {eps}")
    for depth in args.depths:
        if args.test:
            evaluate(f"test {depth}", [load(c, depth) for c in args.test], T, eps)


if __name__ == "__main__":
    main()
