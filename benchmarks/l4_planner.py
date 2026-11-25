"""L4 (`docs/design/l4-calibration-planner.md`): calibration-size planner and its check.

For the K1b pool (agreeing pairs) and the W1b pool (all pairs), 30x, both phasers, every
group: the unit rate `r` and units per genome at the full pool's first grid threshold, the
units and genomes needed for certification probability 0.8 and 0.95, and the bootstrap check
of predicted against observed certification frequency for pool sizes 1 to 6.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/l4_planner.py [--draws 400]`
"""

from __future__ import annotations

import argparse
import math
import random
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import confidence_score as cs  # noqa: E402

from phasebudget.confidence import (binomial_cdf, certification_probability,  # noqa: E402
                                    learn_then_test_by_group, units_needed)

DEPTH = "30"
GROUPS = ("easy", "difficult", "segdup", "lowmap", "tandem", "indel")
SIZES = range(1, 7)
POWERS = (0.8, 0.95)
SEED = 20261122


def k1b_rows(phaser: str) -> dict:
    """Per genome, K1b's agreeing `(group, unit, score, wrong)` rows."""
    import k1b_score as k
    from stratified_gate import strata

    stratum = strata(list(k.CHROMS))
    out = {}
    for sample, platform in k.POOL.items():
        with ThreadPoolExecutor(8) as ex:
            tables = dict(zip(k.CHROMS, ex.map(k.load_table, [(platform, c) for c in k.CHROMS])))
        rows = cs.load_rows(platform, DEPTH, phaser, list(k.CHROMS))
        out[sample] = [r[:4] for r in k.pair_rows(sample, rows, tables, phaser, stratum) if r[4]]
    return out


def w1b_rows(phaser: str) -> dict:
    """Per genome, W1b's pool rows (HG002 and the 20 HPRC pool genomes)."""
    import genome_common as common
    from stratified_gate import strata
    from stratified_score import grouped
    from w1_score import CAL
    from w1b_score import POOL

    chroms = list(common.CONFIRMATORY)
    stratum = strata(chroms)
    platforms = {"HG002": CAL, **{s: f"hprc_{s}" for s in POOL}}
    return {s: grouped(p, DEPTH, phaser, chroms, stratum) for s, p in platforms.items()}


def arrays(rows_by_genome: dict, group: str) -> list:
    """Per genome `(scores, best per unit, best wrong per unit)`, scores as integer ranks."""
    keys = sorted({k for rows in rows_by_genome.values() for g, _, k, _ in rows if g == group})
    rank = {k: i for i, k in enumerate(keys)}
    out = []
    for rows in rows_by_genome.values():
        best, worst, scores = {}, {}, []
        for g, u, k, w in rows:
            if g != group:
                continue
            x = rank[k]
            scores.append(x)
            best[u] = max(best.get(u, -1), x)
            if w:
                worst[u] = max(worst.get(u, -1), x)
        out.append((np.array(scores, dtype=np.int64), np.sort(np.fromiter(best.values(), np.int64)),
                    np.sort(np.fromiter(worst.values(), np.int64))))
    return out


def first_threshold(genomes: list) -> tuple:
    """`(t0, n per genome, k per genome)` at `coverage_grid`'s first threshold (half trusted)."""
    scores = np.concatenate([g[0] for g in genomes])
    total = len(scores)
    if total == 0:
        return None, [0] * len(genomes), [0] * len(genomes)
    idx = max(0, min(total, int(round(0.5 * total))) - 1)  # rank in descending order
    t0 = np.partition(scores, total - 1 - idx)[total - 1 - idx]
    n = [len(b) - int(np.searchsorted(b, t0, "left")) for _, b, _ in genomes]
    k = [len(w) - int(np.searchsorted(w, t0, "left")) for _, _, w in genomes]
    return t0, n, k


def certifies(n: int, k: int, beta: float) -> bool:
    return n > 0 and binomial_cdf(k, n, beta) <= cs.DELTA


def wilson(x: int, b: int, z: float = 1.959964) -> tuple:
    p = x / b
    c = (p + z * z / (2 * b)) / (1 + z * z / b)
    h = z * math.sqrt(p * (1 - p) / b + z * z / (4 * b * b)) / (1 + z * z / b)
    return c - h, c + h


def check_fast(rows_by_genome: dict, rng: random.Random) -> int:
    """Fast first-threshold verdicts against `learn_then_test_by_group` on a few draws."""
    names = list(rows_by_genome)
    checked = 0
    for m in (1, 2):
        draw = [rng.randrange(len(names)) for _ in range(m)]
        rows = [(g, (j, u), k, w) for j, i in enumerate(draw) for g, u, k, w in rows_by_genome[names[i]]]
        for beta in cs.BETAS:
            slow = learn_then_test_by_group(rows, beta, cs.DELTA)
            for group in GROUPS:
                genomes = arrays({j: rows_by_genome[names[i]] for j, i in enumerate(draw)}, group)
                _, n, k = first_threshold(genomes)
                fast = certifies(sum(n), sum(k), beta)
                assert fast == (slow.get(group) is not None), (draw, beta, group)
                checked += 1
    return checked


def run(item) -> list:
    pool, phaser, draws = item
    rows = (k1b_rows if pool == "K1b" else w1b_rows)(phaser)
    rng = random.Random(f"{SEED}-{pool}-{phaser}")
    checked = check_fast(rows, rng)
    out = [f"  {pool} {phaser}: {len(rows)} genomes, fast verdicts equal learn_then_test_by_group "
           f"in {checked} of {checked} checks"]
    planner, cells = [], []
    full_ltt = {b: learn_then_test_by_group([(g, (s, u), k, w) for s, rs in rows.items()
                                             for g, u, k, w in rs], b, cs.DELTA) for b in cs.BETAS}
    for group in GROUPS:
        genomes = arrays(rows, group)
        _, n_full, k_full = first_threshold(genomes)
        N, K = sum(n_full), sum(k_full)
        r = K / max(1, N)
        u = N / len(genomes)
        for beta in cs.BETAS:
            assert certifies(N, K, beta) == (full_ltt[beta].get(group) is not None)
            need = [units_needed(r, beta, cs.DELTA, p) for p in POWERS]
            planner.append((pool, phaser, group, beta, len(genomes), N, K, r, u,
                            certifies(N, K, beta), need,
                            [None if x is None or u == 0 else math.ceil(x / u) for x in need]))
        for m in SIZES:
            obs = {b: 0 for b in cs.BETAS}
            pred = {b: 0.0 for b in cs.BETAS}
            for _ in range(draws):
                pick = [genomes[rng.randrange(len(genomes))] for _ in range(m)]
                _, n, k = first_threshold(pick)
                for beta in cs.BETAS:
                    obs[beta] += certifies(sum(n), sum(k), beta)
                    pred[beta] += certification_probability(sum(n), r, beta, cs.DELTA)
            for beta in cs.BETAS:
                cells.append((pool, phaser, group, beta, m, obs[beta], pred[beta] / draws))
    print(f"  {pool} {phaser} done", flush=True)
    return out, planner, cells


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draws", type=int, default=400)
    args = parser.parse_args()
    started = time.time()
    print(f"L4 planner check (docs/design/l4-calibration-planner.md), development: both pools "
          f"seen before. Depth {DEPTH}x; phasers {cs.PHASERS}; groups {GROUPS}; betas {cs.BETAS}; "
          f"delta {cs.DELTA}; unit {cs.UNIT // 1000} kb; K1b pool HG002 + 5 Pedigree HiFi, "
          f"agreeing pairs; W1b pool HG002 + 20 HPRC, all pairs, consensus truth; bootstrap "
          f"pools of {list(SIZES)} genomes with replacement, B = {args.draws}, seed {SEED}; "
          f"prediction (1) at the full pool's unit rate r and each draw's n; agreement = "
          f"prediction inside the 95% Wilson interval of the observed frequency; informative = "
          f"prediction in [0.02, 0.98] or observed strictly in (0, 1); rule: >= 80% of "
          f"informative cells agree; powers {POWERS}", flush=True)
    items = [(p, ph, args.draws) for p in ("K1b", "W1b") for ph in cs.PHASERS]
    with ProcessPoolExecutor(len(items)) as ex:
        results = list(ex.map(run, items))
    for out, _, _ in results:
        print("\n".join(out))

    print("\nPlanner: pool phaser group beta | genomes, units at t0 (N), wrong units (K), r = K/N, "
          "units per genome u, full pool certifies | units for power 0.8 / 0.95 | genomes")
    for _, planner, _ in results:
        for pool, ph, g, b, gen, N, K, r, u, cert, need, genomes in planner:
            fmt = lambda xs: " / ".join("never" if x is None else f"{x:,}" for x in xs)  # noqa: E731
            print(f"  {pool} {ph:<9} {g:<9} {b:<6g} | {gen:>2} {N:>8,} {K:>6} {r:.2e} "
                  f"{u:>9,.0f} {'cert' if cert else 'none':<4} | {fmt(need):>19} | {fmt(genomes)}")

    print("\nCells: pool phaser group beta m | observed certifying / B, 95% interval, predicted, "
          "agree, informative")
    total = informative = agree = 0
    worst = (0.0, None)
    for _, _, cells in results:
        for pool, ph, g, b, m, x, p in cells:
            lo, hi = wilson(x, args.draws)
            ok = lo <= p <= hi
            info = 0.02 <= p <= 0.98 or 0 < x < args.draws
            total += 1
            informative += info
            agree += ok and info
            diff = abs(x / args.draws - p)
            if diff > worst[0]:
                worst = (diff, (pool, ph, g, b, m, x / args.draws, p))
            print(f"  {pool} {ph:<9} {g:<9} {b:<6g} {m} | {x:>4}/{args.draws} [{lo:.3f}, {hi:.3f}] "
                  f"{p:.3f} {'agree' if ok else 'DISAGREE'} {'info' if info else 'trivial'}")
    share = agree / max(1, informative)
    print(f"\ncells {total}; informative {informative}; informative agreeing {agree} ({share:.3f}); "
          f"trivial {total - informative}; worst |observed - predicted| {worst[0]:.3f} at "
          f"{worst[1]}")
    print(f"L4 rule (>= 0.80 of informative cells agree): {'holds' if share >= 0.8 else 'fails'}")
    print(f"runtime {time.time() - started:.0f} s")


if __name__ == "__main__":
    main()
