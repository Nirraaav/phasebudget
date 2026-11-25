"""Round 4, W3 with pooled calibration (development): learned against frozen score, genome held out.

W3's gate failed on one genome (`cal_r10`, then `cal_v5`) because `segdup` has
about 2,900 units per genome, below the 2,995 that Learn then Test needs to
certify beta 1e-3 at delta 0.05 whatever the score. This run pools the
calibration over genomes (`docs/design/w1b-pooled-plan.md` states the pooled
procedure and its risk), on ONT data already scored in earlier tests:
HG002 (`cal_v5`), HG001 (`hg001_v5`) and five Platinum Pedigree children
(`ped_NA12879` to `ped_NA12886`; HG001 is their mother).

For each depth and phaser, one model (as in `w3_gate.py`) trains on
chromosomes 1 to 6 of all seven genomes. Each genome in turn is held out:
thresholds are calibrated per group on chromosomes 7 to 22 of the other six,
pooled, with a unit identified by (genome, chromosome, 100 kb index), and
scored on chromosomes 7 to 22 of the held-out genome. The frozen score is
calibrated and scored the same way. The held-out genome shares parents or
children with the pool, so this is easier than an unrelated genome.

Gate (W3's, with pooled calibration): at beta 1e-3 and depths 30 and 20, for
both phasers, pooled over held-out genomes, the learned score trusts at least
1.5 times the frozen score's share of pairs in `segdup` and in `indel` (and
more than none), at most 0.005 less in `easy`, and no cell violates at 0.05
over the cells scored.

Usage: `PYTHONPATH=src:benchmarks .venv-w3/bin/python benchmarks/w3_pooled.py`
"""

from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import confidence_score as cs  # noqa: E402
import genome_common as common  # noqa: E402
import stratified_score as ss  # noqa: E402
import w3_gate as w3  # noqa: E402
from w3_features import NAMES  # noqa: E402

from phasebudget.confidence import learn_then_test_by_group  # noqa: E402

GENOMES = ("cal_v5", "hg001_v5", "ped_NA12879", "ped_NA12881", "ped_NA12882", "ped_NA12885",
           "ped_NA12886")
LATER = tuple(c for c in common.CONFIRMATORY if c not in w3.TRAIN)
SCORES = ("frozen", "learned", "neural")
NEURAL = HERE.parent / "data" / "w3_neural"


def cell(item) -> list:
    """`[(held-out genome, beta, score, group, cells tuple)]` for one depth and phaser."""
    import numpy as np
    from sklearn.ensemble import HistGradientBoostingClassifier
    from stratified_gate import strata

    depth, phaser = item
    stratum = strata(list(common.CONFIRMATORY))
    data = {g: w3.load(depth, phaser, stratum, g) for g in GENOMES}
    train = [np.isin(chrom, w3.TRAIN) for _, chrom, *_ in data.values()]
    model = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05,
                                           l2_regularization=1.0, random_state=0)
    model.fit(np.vstack([d[0][t] for d, t in zip(data.values(), train)]),
              np.concatenate([d[3][t] for d, t in zip(data.values(), train)]))
    w = NAMES.index("wilson"), NAMES.index("reads")
    # W3b: the neural score of `w3_neural.py`, when it has been written for every genome.
    files = [NEURAL / f"{g}_{depth}_{phaser}.npy" for g in GENOMES]
    neural = {g: np.load(f) for g, f in zip(GENOMES, files)} if all(f.exists() for f in files) \
        else None
    rows = {}  # (genome, score) -> [(group, unit, key, wrong)] on chromosomes 7 to 22
    for g, (x, chrom, pos1, wrong, group) in data.items():
        later = np.isin(chrom, LATER)
        p = model.predict_proba(x[later])[:, 1]
        if neural is not None:
            assert len(neural[g]) == len(x), g
            rows[(g, "neural")] = [
                (str(gr), (g, str(c), int(a) // cs.UNIT), (-float(q), float(s), int(n)), bool(k))
                for gr, c, a, q, s, n, k in zip(group[later], chrom[later], pos1[later],
                                                neural[g][later], x[later][:, w[0]],
                                                x[later][:, w[1]], wrong[later])]
        x, chrom, pos1, wrong, group = x[later], chrom[later], pos1[later], wrong[later], group[later]
        units = [(g, str(c), int(a) // cs.UNIT) for c, a in zip(chrom, pos1)]
        rows[(g, "frozen")] = [(str(gr), u, (float(a), int(b)), bool(k)) for gr, u, a, b, k in
                               zip(group, units, x[:, w[0]], x[:, w[1]], wrong)]
        rows[(g, "learned")] = [(str(gr), u, (-float(q), float(a), int(b)), bool(k))
                                for gr, u, q, a, b, k in
                                zip(group, units, p, x[:, w[0]], x[:, w[1]], wrong)]
    out = []
    for held in GENOMES:
        for score in [s for s in SCORES if (held, s) in rows]:
            cal = [r for g in GENOMES if g != held for r in rows[(g, score)]]
            for beta in cs.BETAS:
                t = learn_then_test_by_group(cal, beta, cs.DELTA)
                for group, c in ss.cells(t, rows[(held, score)], beta).items():
                    out.append((depth, phaser, held, beta, score, group, c))
    print(f"  {depth}x {phaser} done", flush=True)
    return out


def main() -> None:
    print(f"W3 pooled (development): genomes {', '.join(GENOMES)}; model trains on {w3.TRAIN} "
          f"of all genomes; per held-out genome, calibrate on {LATER[0]} to {LATER[-1]} of the "
          f"others pooled, score the held-out genome's; {cs.UNIT // 1000} kb units, delta "
          f"{cs.DELTA}, betas {cs.BETAS}", flush=True)
    items = [(d, ph) for d in cs.DEPTHS for ph in cs.PHASERS]
    with ProcessPoolExecutor(len(items)) as pool:
        results = [r for part in pool.map(cell, items) for r in part]
    level = 0.05 / max(1, len(results))
    bad = [r for r in results if r[6][7] and r[6][0] and r[6][2] < level]
    scores = [s for s in SCORES if any(r[4] == s for r in results)]
    print(f"\nper depth, phaser, beta, group, pooled over held-out genomes: certifying cells, "
          f"trusted share, wrong units/units, each for {', '.join(scores)}; violation level "
          f"{level:.2e}")
    gates = {s: [] for s in scores if s != "frozen"}
    for depth, phaser in items:
        for beta in cs.BETAS:
            for group in w3.GROUPS:
                summary = {}
                for score in scores:
                    mine = [r[6] for r in results if r[:2] == (depth, phaser) and r[3] == beta
                            and r[4] == score and r[5] == group]
                    summary[score] = (sum(c[7] for c in mine), len(mine),
                                      sum(c[4] for c in mine) / max(1, sum(c[3] for c in mine)),
                                      sum(c[1] for c in mine), sum(c[0] for c in mine))
                v = [summary[s] for s in scores]
                print(f"  {depth:>2}x {phaser:<9} {beta:<7g} {group:<10} "
                      + " ".join(f"{x[0]}/{x[1]}" for x in v) + "  "
                      + " ".join(f"{x[2]:.4f}" for x in v) + "  "
                      + " ".join(f"{x[3]}/{x[4]}" for x in v))
                f = summary["frozen"]
                for s, gate in gates.items():
                    ln = summary[s]
                    if beta == w3.GATE_BETA and depth in w3.GATE_DEPTHS:
                        if group in ("segdup", "indel"):
                            gate.append((depth, phaser, group,
                                         ln[2] >= w3.GATE_RATIO * f[2] and ln[2] > 0))
                        elif group == "easy":
                            gate.append((depth, phaser, group, ln[2] >= f[2] - w3.EASY_LOSS))
    print(f"\nviolating cells: {len(bad)}"
          + "".join(f"\n  {d}x {ph} held out {h} beta {b:g} {s} {g} p {c[2]:.2e} ({c[1]}/{c[0]})"
                    for d, ph, h, b, s, g, c in bad))
    for s, gate in gates.items():
        mine = [r for r in bad if r[4] == s]
        print(f"{s} gate conditions: " + "; ".join(f"{d}x {ph} {g} {'ok' if ok else 'FAIL'}"
                                                  for d, ph, g, ok in gate))
        print(f"{s} gate: " + ("passed" if all(ok for *_, ok in gate) and not mine
                               else "not passed"))


if __name__ == "__main__":
    main()
