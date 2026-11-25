"""Round 4, W3 gate: does a learned per-pair score certify more than the frozen score?

The guarantee of Learn then Test holds for any score fixed before
calibration, so the score can be learned, provided its training data are
neither calibration nor test data. Chromosomes are split three ways: the
model trains on chromosomes 1 to 6, thresholds are calibrated per group on
the odd chromosomes from 7 on, and both are scored on the even ones from 8 on
against the truth. The frozen score (Wilson bound, reads) is calibrated and
scored on the same split, so the two are compared on the same pairs.

The model is scikit-learn's `HistGradientBoostingClassifier` on the
features of `w3_features.py` plus the group, one model per depth and phaser,
predicting a wrong pair; the learned score is (minus the predicted
probability, Wilson bound, reads), compared in that order.

Gate (round 4): at beta 1e-3 and depths 30 and 20, for both phasers, the
learned score trusts at least 1.5 times the frozen score's share of pairs in
`segdup` and in `indel`, at most 0.005 less in `easy`, and no cell violates
at 0.05 over the cells scored.

Usage: `PHASEBUDGET_PLATFORM=cal_r10 PYTHONPATH=src:benchmarks .venv-w3/bin/python
benchmarks/w3_gate.py`
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import confidence_score as cs  # noqa: E402
import genome_common as common  # noqa: E402
import stratified_score as ss  # noqa: E402
from w3_features import NAMES, OUT  # noqa: E402

from phasebudget.confidence import learn_then_test_by_group  # noqa: E402

GROUPS = ("easy", "difficult", "lowmap", "segdup", "tandem", "indel")
TRAIN = tuple(f"chr{i}" for i in range(1, 7))
CAL = tuple(c for c in common.CONFIRMATORY if c not in TRAIN and int(c[3:]) % 2)
TEST = tuple(c for c in common.CONFIRMATORY if c not in TRAIN and int(c[3:]) % 2 == 0)
GATE_BETA, GATE_DEPTHS, GATE_RATIO, EASY_LOSS = 1e-3, ("30", "20"), 1.5, 0.005


def load(depth: str, phaser: str, stratum, platform: str = ""):
    import numpy as np

    data = np.load(OUT / f"{platform or common.PLATFORM}_{depth}_{phaser}.npz")
    assert tuple(data["names"]) == NAMES
    x = data["x"]
    chrom, pos1, pos2, wrong = data["chrom"], data["pos1"], data["pos2"], data["wrong"]
    ia, ib = NAMES.index("indel_a"), NAMES.index("indel_b")
    group = np.array(["indel" if x[i, ia] or x[i, ib] else stratum(chrom[i], int(pos1[i]),
                                                                   int(pos2[i]))
                      for i in range(len(x))])
    onehot = np.stack([group == g for g in GROUPS], axis=1).astype(float)
    return np.hstack([x, onehot]), chrom, pos1, wrong, group


def main() -> None:
    import numpy as np
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score
    from stratified_gate import strata

    stratum = strata(list(common.CONFIRMATORY))
    print(f"W3 gate, platform {common.PLATFORM} (development unless stated): train on {TRAIN}, "
          f"calibrate on {CAL}, score on {TEST}; {cs.UNIT // 1000} kb units, delta {cs.DELTA}, "
          f"betas {cs.BETAS}; HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, "
          "l2_regularization=1.0, random_state=0)")
    w = NAMES.index("wilson"), NAMES.index("reads")
    results = []  # (depth, phaser, beta, score, group, cell)
    for depth in cs.DEPTHS:
        for phaser in cs.PHASERS:
            x, chrom, pos1, wrong, group = load(depth, phaser, stratum)
            train = np.isin(chrom, TRAIN)
            model = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05,
                                                   l2_regularization=1.0, random_state=0)
            model.fit(x[train], wrong[train])
            p = model.predict_proba(x)[:, 1]
            keys = {"frozen": [(float(a), int(b)) for a, b in zip(x[:, w[0]], x[:, w[1]])],
                    "learned": [(-float(q), float(a), int(b))
                                for q, a, b in zip(p, x[:, w[0]], x[:, w[1]])]}
            test = np.isin(chrom, TEST)
            print(f"\n{depth}x {phaser}: train {train.sum()} pairs ({wrong[train].sum()} wrong), "
                  f"calibrate {np.isin(chrom, CAL).sum()}, test {test.sum()} "
                  f"({wrong[test].sum()} wrong)")
            print("  group       test_pairs  wrong  AUC_frozen  AUC_learned")
            for g in GROUPS:
                m = test & (group == g)
                if wrong[m].sum() and (~wrong[m]).sum():
                    frozen = roc_auc_score(wrong[m], -x[m, w[0]] - 1e-6 * x[m, w[1]])
                    learned = roc_auc_score(wrong[m], p[m])
                    print(f"  {g:<10} {m.sum():>11} {wrong[m].sum():>6} {frozen:>11.3f}"
                          f" {learned:>12.3f}")
            for name, key in keys.items():
                rows = [(str(group[i]), (str(chrom[i]), int(pos1[i]) // cs.UNIT), key[i],
                         bool(wrong[i])) for i in range(len(key))]
                cal = [r for r, c in zip(rows, chrom) if c in CAL]
                scored = [r for r, c in zip(rows, chrom) if c in TEST]
                for beta in cs.BETAS:
                    t = learn_then_test_by_group(cal, beta, cs.DELTA)
                    for g, cell in ss.cells(t, scored, beta).items():
                        results.append((depth, phaser, beta, name, g, cell))
    level = 0.05 / max(1, len(results))
    print(f"\nper cell: trusted share of the group's test pairs, wrong units/units; violation "
          f"level {level:.2e}")
    print("depth phaser     beta    group       frozen_share  learned_share  frozen_wrong_units"
          "  learned_wrong_units")
    index = {(d, ph, b, n, g): c for d, ph, b, n, g, c in results}
    bad = [r for r in results if r[5][7] and r[5][0] and r[5][2] < level]
    gate = []
    for depth in cs.DEPTHS:
        for phaser in cs.PHASERS:
            for beta in cs.BETAS:
                for g in GROUPS:
                    f, ln = index.get((depth, phaser, beta, "frozen", g)), \
                        index.get((depth, phaser, beta, "learned", g))
                    if f is None or ln is None:
                        continue
                    fs, ls = f[4] / max(1, f[3]), ln[4] / max(1, ln[3])
                    print(f"{depth:>5} {phaser:<10} {beta:<7g} {g:<10} {fs:>12.4f} {ls:>14.4f}"
                          f"  {f[1]:>8}/{f[0]:<9} {ln[1]:>8}/{ln[0]:<9}")
                    if beta == GATE_BETA and depth in GATE_DEPTHS:
                        if g in ("segdup", "indel"):
                            gate.append((depth, phaser, g, ls >= GATE_RATIO * fs and ls > 0))
                        elif g == "easy":
                            gate.append((depth, phaser, g, ls >= fs - EASY_LOSS))
    print(f"\nviolating cells: {len(bad)}"
          + "".join(f"\n  {d}x {ph} beta {b:g} {n} {g} p {c[2]:.2e}" for d, ph, b, n, g, c in bad))
    print("gate conditions: " + "; ".join(f"{d}x {ph} {g} {'ok' if ok else 'FAIL'}"
                                           for d, ph, g, ok in gate))
    print("gate: " + ("passed" if all(ok for *_, ok in gate) and not bad else "not passed"))


if __name__ == "__main__":
    main()
