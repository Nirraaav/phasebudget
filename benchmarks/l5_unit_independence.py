"""L5 (`docs/design/l5-unit-independence.md`): do wrong 100 kb units cluster, per group?

Part 1, per (platform, depth, phaser, group): joins of adjacent wrong units against a
within-chromosome permutation null. Part 2: Learn then Test thresholds of the K1b and W1b
pools with units widened to 200 and 300 kb. Reads only existing row caches.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/l5_unit_independence.py`
"""

from __future__ import annotations

import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import confidence_score as cs  # noqa: E402

GROUPS = ("indel", "segdup", "lowmap", "tandem", "difficult", "easy")
ALL3 = ("30", "20", "10")
PEDS = ("NA12879", "NA12881", "NA12882", "NA12885", "NA12886")
PLATFORMS = {"hprc_HG002": ALL3, "w4_grch38": ("30",), "cal_v5": ALL3, "cal_r10": ALL3,
             "hg001_v5": ALL3, **{f"ped_{s}": ALL3 for s in PEDS},
             **{f"phifi_{s}": ALL3 for s in PEDS}, "phifi_NA12877": ("30",),
             "k1_hg005": ("30",), "hg005_ont": ALL3}
TESTS = sum(len(d) for d in PLATFORMS.values()) * len(cs.PHASERS) * len(GROUPS)
LEVEL = 0.05 / TESTS
B_FIRST, B_MORE, SEED = 2000, 40000, 0
MERGES = (2, 3)


def unit_flags(rows) -> dict:
    """`{chrom: {unit index: wrong}}` from `(unit=(chrom, index), wrong)` rows."""
    out: dict = {}
    for (chrom, idx), wrong in rows:
        units = out.setdefault(chrom, {})
        units[idx] = units.get(idx, False) or bool(wrong)
    return out


def joins(flags: dict, b: int, rng) -> tuple:
    """`(observed joins, null joins array of length b, adjacent pairs, rho)`."""
    observed, null, adjacent, xs, ys = 0, np.zeros(b, dtype=np.int64), 0, [], []
    for units in flags.values():
        idx = sorted(units)
        w = np.array([units[i] for i in idx], dtype=bool)
        pos = np.array(idx)
        left = np.nonzero(np.diff(pos) == 1)[0]
        adjacent += len(left)
        xs.append(w[left])
        ys.append(w[left + 1])
        if not len(left) or w.sum() < 2:
            continue
        observed += int((w[left] & w[left + 1]).sum())
        for start in range(0, b, 1000):
            m = min(1000, b - start)
            perm = rng.permuted(np.tile(w, (m, 1)), axis=1)
            null[start:start + m] += (perm[:, left] & perm[:, left + 1]).sum(axis=1)
    x = np.concatenate(xs) if xs else np.zeros(0, bool)
    y = np.concatenate(ys) if ys else np.zeros(0, bool)
    rho = float(np.corrcoef(x, y)[0, 1]) if len(x) > 1 and x.std() > 0 and y.std() > 0 else 0.0
    return observed, null, adjacent, rho


def test(flags: dict) -> dict:
    """The join-count test of the design note on one group's units."""
    n = sum(len(u) for u in flags.values())
    k = sum(sum(u.values()) for u in flags.values())
    obs, null, adj, rho = joins(flags, B_FIRST, np.random.default_rng(SEED))
    p = (1 + int((null >= obs).sum())) / (1 + B_FIRST)
    if p < 0.01:
        obs, null, adj, rho = joins(flags, B_MORE, np.random.default_rng(SEED))
        p = (1 + int((null >= obs).sum())) / (1 + B_MORE)
    neff = n * (1 - rho) / (1 + rho) if rho > 0 else n
    return {"units": n, "wrong": k, "adjacent": adj, "obs": obs, "mean": float(null.mean()),
            "p": p, "rho": rho, "neff": neff}


def part1(job) -> list:
    from stratified_gate import strata
    from stratified_score import grouped

    import genome_common as common

    platform, depth, phaser = job
    chroms = list(common.CONFIRMATORY)
    rows = grouped(platform, depth, phaser, chroms, strata(chroms))
    out = []
    for g in GROUPS:
        r = test(unit_flags((u, w) for grp, u, _, w in rows if grp == g))
        out.append((platform, depth, phaser, g, r))
    print(f"  part 1 {platform} {depth}x {phaser} done", flush=True)
    return out


def widen(rows, m: int) -> list:
    """`(group, unit, score, wrong)` rows with the unit's 100 kb index `// m`.

    The index is the unit's last field, or the last field of a nested last field: W1b's
    pooled units are `(genome, (chrom, index))`, K1b's are `(genome, chrom, index)`."""

    def coarse(u):
        return (*u[:-1], coarse(u[-1]) if isinstance(u[-1], tuple) else u[-1] // m)

    return [(g, coarse(u), k, w) for g, u, k, w in rows]


def compare(rows) -> list:
    """`[(beta, group, merge, base threshold, merged threshold, conservative)]`."""
    from phasebudget.confidence import learn_then_test_by_group

    out = []
    for beta in cs.BETAS:
        base = learn_then_test_by_group(rows, beta, cs.DELTA)
        for m in MERGES:
            merged = learn_then_test_by_group(widen(rows, m), beta, cs.DELTA)
            for g in GROUPS:
                t1, tm = base.get(g), merged.get(g)
                ok = tm is None or (t1 is not None and tm >= t1)
                out.append((beta, g, m, t1, tm, ok))
    return out


def k1b_part(phaser: str) -> list:
    import k1b_score as kb
    from stratified_gate import strata

    chroms = list(kb.CHROMS)
    stratum = strata(chroms)
    cal = []
    for label, platform in kb.POOL.items():
        tables = {c: kb.load_table((platform, c)) for c in chroms}
        rows = cs.load_rows(platform, kb.DEPTH, phaser, chroms)
        cal += kb.pair_rows(label, rows, tables, phaser, stratum)
    out = []
    for mode in ("agree", "all"):
        sel = [r[:4] for r in cal if mode == "all" or r[4]]
        out += [("K1b", "30", phaser, mode, *c) for c in compare(sel)]
    print(f"  part 2 K1b {phaser} done", flush=True)
    return out


def w1b_part(job) -> list:
    import w1b_score as wb
    from stratified_gate import strata
    from stratified_score import grouped
    from w1_score import CAL

    import genome_common as common

    depth, phaser = job
    chroms = list(common.CONFIRMATORY)
    stratum = strata(chroms)
    rows = {"HG002": grouped(CAL, depth, phaser, chroms, stratum)}
    rows.update({s: grouped(f"hprc_{s}", depth, phaser, chroms, stratum) for s in wb.POOL})
    out = [("W1b", depth, phaser, "all", *c) for c in compare(wb.pooled(rows))]
    print(f"  part 2 W1b {depth}x {phaser} done", flush=True)
    return out


def fmt(t) -> str:
    return "none" if t is None else f"({t[0]:.4f},{t[1]})"


def main() -> None:
    started = time.time()
    print(f"L5 unit independence (docs/design/l5-unit-independence.md, development: data seen); "
          f"unit {cs.UNIT // 1000} kb per group, wrong = any wrong pair; statistic = joins of "
          f"adjacent wrong units; null = within-chromosome permutation, B {B_FIRST} then {B_MORE} "
          f"when p < 0.01, seed {SEED}; {TESTS} tests, level 0.05/{TESTS} = {LEVEL:.3e}; "
          f"platforms {', '.join(p + ' (' + '/'.join(d) + ')' for p, d in PLATFORMS.items())}; phasers "
          f"{cs.PHASERS}; groups {GROUPS}; merge check widen x{MERGES} on K1b pool (30x, agree and "
          f"all) and W1b pool (depths {cs.DEPTHS}), betas {cs.BETAS}, delta {cs.DELTA}", flush=True)
    jobs = [(p, d, ph) for p, ds in PLATFORMS.items() for d in ds for ph in cs.PHASERS]
    with ProcessPoolExecutor(16) as ex:
        p1 = [r for part in ex.map(part1, jobs) for r in part]
        k1 = ex.map(k1b_part, cs.PHASERS)
        w1 = ex.map(w1b_part, [(d, ph) for d in cs.DEPTHS for ph in cs.PHASERS])
        p2 = [r for part in k1 for r in part] + [r for part in w1 for r in part]

    print("\nPart 1: platform depth phaser group units wrong adjacent joins null_mean ratio p rho "
          "neff/n clustered")
    clustered = []
    for platform, depth, phaser, g, r in p1:
        ratio = r["obs"] / r["mean"] if r["mean"] > 0 else float("nan")
        bad = r["p"] <= LEVEL
        if bad:
            clustered.append((platform, depth, phaser, g))
        print(f"  {platform:<14} {depth:>2} {phaser:<9} {g:<9} {r['units']:>6} {r['wrong']:>5} "
              f"{r['adjacent']:>6} {r['obs']:>5} {r['mean']:>8.2f} {ratio:>6.2f} {r['p']:.2e} "
              f"{r['rho']:+.4f} {r['neff'] / max(1, r['units']):.4f}{'  CLUSTERED' if bad else ''}")
    worst = min(p1, key=lambda x: x[4]["neff"] / max(1, x[4]["units"]))
    print(f"\nclustered tests: {len(clustered)} of {len(p1)}")
    for g in GROUPS:
        print(f"  {g:<9} clustered {sum(c[3] == g for c in clustered)} of "
              f"{sum(x[3] == g for x in p1)}")
    print(f"worst neff/n: {worst[4]['neff'] / max(1, worst[4]['units']):.4f} "
          f"({' '.join(worst[:4])}, rho {worst[4]['rho']:+.4f})")
    ratios = [x[4]["obs"] / x[4]["mean"] for x in p1 if x[4]["mean"] > 0]
    print(f"join ratio observed/null over tests with a null mean > 0: mean "
          f"{sum(ratios) / max(1, len(ratios)):.3f}, max {max(ratios, default=float('nan')):.3f}")

    print("\nPart 2: pool depth phaser mode beta group widen base merged conservative")
    for pool, depth, phaser, mode, beta, g, m, t1, tm, ok in p2:
        print(f"  {pool} {depth:>2} {phaser:<9} {mode:<5} {beta:<6g} {g:<9} x{m} {fmt(t1):<18} "
              f"{fmt(tm):<18} {'yes' if ok else 'NO'}")
    flagged = sorted({c[3] for c in clustered})
    looser = [r for r in p2 if not r[-1]]
    print(f"\nmerged cells less conservative, all groups: {len(looser)} of {len(p2)}")
    looser_flagged = [r for r in looser if r[5] in flagged]
    if not clustered:
        verdict = "supported by (a): no clustered test"
    elif not looser_flagged:
        verdict = f"supported by (b): clustering in {flagged}, merged cells never looser there"
    else:
        verdict = (f"not supported: clustering in {flagged}, {len(looser_flagged)} merged cells "
                   "looser there")
    print(f"L5: {verdict}")
    print(f"runtime {time.time() - started:.0f} s")


if __name__ == "__main__":
    main()
