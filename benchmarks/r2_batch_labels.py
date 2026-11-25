"""Round 5, R2 (development): certifying a held-out batch from k of its own labelled units.

`docs/design/r2-batch-labels.md` states the floor, the methods and the gate.
The new batch is HG002 (`cal_v5`); the pool is the other six genomes of
`benchmarks/w3_pooled.py`. Frozen score, read from the cached rows.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/r2_batch_labels.py`
"""

from __future__ import annotations

import math
import random
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import confidence_score as cs  # noqa: E402
import genome_common as common  # noqa: E402
import stratified_score as ss  # noqa: E402

from phasebudget.confidence import Units, coverage_grid, learn_then_test  # noqa: E402

TARGET = "cal_v5"
POOL = ("hg001_v5", "ped_NA12879", "ped_NA12881", "ped_NA12882", "ped_NA12885", "ped_NA12886")
DEPTHS = ("30", "20")
KS = (300, 1000, 3000)
SEEDS = 20
METHODS = ("pool", "target", "checked", "merged")
VALID = ("target", "checked")


def floor_units(beta: float, delta: float, power: float = 1.0) -> float:
    """Labelled units a valid procedure needs to certify `beta` with probability `power` when
    every label is correct: `ln(power / delta) / -ln(1 - beta)`."""
    return math.log(power / delta) / -math.log1p(-beta)


def by_group(rows, beta: float, delta: float, floors: dict | None = None) -> dict:
    """Per group, Learn then Test on `rows` `(group, unit, score, wrong)`. With `floors`, a
    group's grid stops at its floor threshold (never looser); a group whose floor is None
    certifies nothing."""
    grouped: dict = {}
    for group, unit, key, wrong in rows:
        grouped.setdefault(group, []).append((unit, key, wrong))
    out = {}
    for group, members in grouped.items():
        grid = coverage_grid([key for _, key, _ in members])
        if floors is not None:
            limit = floors.get(group)
            grid = () if limit is None else tuple(t for t in grid if t >= limit)
        out[group] = learn_then_test(Units.from_pairs(members), grid, beta, delta)
    return out


def rows_of(platform: str, depth: str, phaser: str, stratum) -> list:
    """`(group, (platform, chrom, unit index), score, wrong)` on the confirmatory chromosomes."""
    return [(g, (platform,) + unit, key, wrong)
            for g, unit, key, wrong in ss.grouped(platform, depth, phaser,
                                                  list(common.CONFIRMATORY), stratum)]


def cell(item) -> list:
    """`[(method, k, seed, beta, group, cells tuple)]` for one depth and phaser."""
    from stratified_gate import strata

    depth, phaser = item
    stratum = strata(list(common.CONFIRMATORY))
    pool = [r for p in POOL for r in rows_of(p, depth, phaser, stratum)]
    target = rows_of(TARGET, depth, phaser, stratum)
    odd = [r for r in target if int(r[1][1][3:]) % 2]
    test = [r for r in target if not int(r[1][1][3:]) % 2]
    odd_units = sorted({r[1] for r in odd})
    out = []
    for beta in cs.BETAS:
        t_pool = by_group(pool, beta, cs.DELTA)
        for g, c in ss.cells(t_pool, test, beta).items():
            out.append(("pool", 0, 0, beta, g, c))
        draws = [(k, seed) for k in KS for seed in range(SEEDS)] + [(len(odd_units), 0)]
        for k, seed in draws:
            chosen = set(random.Random(f"{depth}{phaser}{k}{seed}").sample(odd_units, k))
            mine = [r for r in odd if r[1] in chosen]
            methods = [("target", by_group(mine, beta, cs.DELTA)),
                       ("checked", by_group(mine, beta, cs.DELTA, t_pool))]
            if seed == 0:  # merging recalibrates on the whole pool: one seed per k is enough
                methods.append(("merged", by_group(pool + mine, beta, cs.DELTA)))
            for method, t in methods:
                for g, c in ss.cells(t, test, beta).items():
                    out.append((method, k, seed, beta, g, c))
    return [(depth, phaser) + o for o in out]


def main() -> None:
    print(f"R2 (development, docs/design/r2-batch-labels.md): new batch {TARGET} (HG002), pool "
          f"{', '.join(POOL)}; frozen score; depths {DEPTHS}; phasers {cs.PHASERS}; betas "
          f"{cs.BETAS}; delta {cs.DELTA}; k in {KS} ({SEEDS} seeds) and all odd units of "
          f"{TARGET}; scored on all pairs of {TARGET}'s even chromosomes; {cs.UNIT // 1000} kb "
          "units", flush=True)
    for beta in cs.BETAS:
        print(f"  floor at beta {beta:g}: {floor_units(beta, cs.DELTA):.0f} labelled units "
              f"(power 1), {floor_units(beta, cs.DELTA, 0.5):.0f} (power 0.5)")
    with ProcessPoolExecutor(4) as pool:
        results = [r for part in pool.map(cell, [(d, p) for d in DEPTHS for p in cs.PHASERS])
                   for r in part]
    valid = [r for r in results if r[2] in VALID]
    level = 0.05 / max(1, len(valid))
    print(f"\nviolation level {level:.2e} (Bonferroni over {len(valid)} cells of the valid "
          "methods); per method, k, beta, group, pooled over depths, phasers and seeds: share "
          "of runs certifying, violating runs, mean trusted share, wrong units/units")
    groups = sorted({r[6] for r in results})
    all_k = sorted({r[3] for r in results if r[2] != "pool" and r[3] not in KS})
    for method in METHODS:
        for k in ((0,) if method == "pool" else KS + tuple(all_k)):
            for beta in cs.BETAS:
                for g in groups:
                    mine = [r[7] for r in results if r[2] == method and r[3] == k
                            and r[5] == beta and r[6] == g]
                    if not mine:
                        continue
                    cert = [c for c in mine if c[7]]
                    bad = [c for c in cert if c[0] and c[2] < level]
                    share = sum(c[4] / max(1, c[3]) for c in mine) / len(mine)
                    n, kk = sum(c[0] for c in mine), sum(c[1] for c in mine)
                    print(f"  {method:<8} {'all' if k in all_k else k:>5} {beta:<7g} {g:<10} "
                          f"{len(cert) / len(mine):.2f} {len(bad):>4}  {share:.4f}  "
                          f"{kk}/{n} = {kk / max(1, n):.2e}")
    small = [r for r in valid if r[3] in (300, 1000) and r[5] == 1e-3 and r[7][7]]
    big = {(r[0], r[1]) for r in valid if r[3] in all_k and r[5] == 1e-3 and r[6] == "easy"
           and r[7][7]}
    bad = [r for r in valid if r[7][7] and r[7][0] and r[7][2] < level]
    for r in bad:
        print(f"  violating: {r[0]}x {r[1]} {r[2]} k {r[3]} seed {r[4]} beta {r[5]:g} {r[6]} "
              f"p {r[7][2]:.2e} ({r[7][1]}/{r[7][0]})")
    a, b, c = not small, len(big) == len(DEPTHS) * len(cs.PHASERS), not bad
    print(f"\nR2 gate: (a) floor holds {'ok' if a else 'FAIL'} ({len(small)} certifying cells at "
          f"k 300 or 1000, beta 1e-3); (b) all odd units certify easy at 1e-3 in {len(big)} of "
          f"{len(DEPTHS) * len(cs.PHASERS)} cells {'ok' if b else 'FAIL'}; (c) valid "
          f"{'ok' if c else 'FAIL'} ({len(bad)} violating)")
    print(f"R2 gate: {'passed' if a and b and c else 'not passed'}")


if __name__ == "__main__":
    main()
