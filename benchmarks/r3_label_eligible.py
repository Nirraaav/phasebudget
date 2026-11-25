"""Round 5, R3 (development): label-calibrated certificates scored on the pairs they can reach.

`docs/design/r3-label-eligible.md` states the claim and the gate. HG002
(`cal_v5`), labels from both flowcell replicates of HG003's and HG004's
Clair3 calls (`benchmarks/w5_noise.py`). Odd chromosomes calibrate on
eligible pairs against the labels; even chromosomes are scored against the
platform's primary truth (v5.0q for HG002) on eligible pairs and on all pairs.

With `PHASEBUDGET_PLATFORM=hg005_ont` it is R3b's scorer
(`docs/design/r3b-parent-labels-plan.md`): child HG005, labels from HG006 and
HG007, truth GIAB's trio transmission, and the registered verdicts at the end.

Usage: `PHASEBUDGET_PLATFORM=cal_v5 PYTHONPATH=src:benchmarks .venv/bin/python
benchmarks/r3_label_eligible.py`
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import confidence_score as cs  # noqa: E402
import genome_common as common  # noqa: E402
import stratified_score as ss  # noqa: E402
import w5_called as wc  # noqa: E402
import w5_noise as wn  # noqa: E402

from phasebudget.confidence import learn_then_test_by_group  # noqa: E402

FILTERS = ("agree", "all", "double")
R3B_CELLS = 3 * 2 * 3 * 6  # depths, phasers, betas, groups of the agree filter
R3B_LEVEL = 0.05 / R3B_CELLS


def site_reason(chrom: str, pos1: int, ref: str, alt: str, parents: dict) -> str:
    """Why a site has no `agree` label: `uninformative` when, in both replicates, both parents
    are heterozygous or not PASS there (no label is possible from these calls); else `other`."""
    for father, mother in parents.values():
        key = (chrom, pos1, ref, alt)
        if father.get(key, 0) != 1 or mother.get(key, 0) != 1:
            return "other"
    return "uninformative"


def pair_class(a: str, b: str) -> str:
    """A pair's class from its sites' classes: eligible only when both are."""
    if a == b == "eligible":
        return "eligible"
    return "uninformative" if "uninformative" in (a, b) else "other"


def main() -> None:
    from stratified_gate import strata

    sc = common._module("score_hg002.py")
    chroms = list(common.CONFIRMATORY)
    odd = {c for c in chroms if int(c[3:]) % 2}
    stratum = strata(chroms)
    ind = ss.indels(common.PLATFORM, chroms)
    truth_name = "v5.0q" if common.SAMPLE == "HG002" else "GIAB trio transmission"
    print(f"R3 (development, docs/design/r3-label-eligible.md): child {common.SAMPLE} ({common.PLATFORM}); "
          f"labels from both replicates of the parents' Clair3 calls; filters {FILTERS}; "
          "calibrate on eligible odd-chromosome pairs against the labels, score even-chromosome "
          f"pairs against {truth_name}, eligible only and all; {cs.UNIT // 1000} kb units, delta "
          f"{cs.DELTA}; betas {cs.BETAS}; depths {cs.DEPTHS}; phasers {cs.PHASERS}", flush=True)
    parents = {name: (wc.calls(fa, set(chroms)), wc.calls(mo, set(chroms)))
               for name, (fa, mo) in wc.REPLICATES.items()}
    truth, labels, reason = {}, {}, {}
    for chrom in chroms:
        variants = sc.variant_table(common.hetsnv_path(chrom))
        for v, t in zip(variants, common.primary_truth(chrom, variants)):
            if t is not None:
                truth[(chrom, v[0])] = t
        per = {name: wn.labelled(chrom, variants, *parents[name]) for name in wc.REPLICATES}
        for pos1, ref, alt in variants:
            a, b = per["A"].get(pos1), per["B"].get(pos1)
            labels[(chrom, pos1)] = (a, b)
            agree = wn.FILTERS["agree"](a, b)
            reason[(chrom, pos1)] = ("eligible" if agree else
                                     site_reason(chrom, pos1, ref, alt, parents))

    def relation(alleles, c, a, b):
        x, y = alleles.get((c, a)), alleles.get((c, b))
        return None if x is None or y is None else x == y

    summary = []
    for depth in cs.DEPTHS:
        for phaser in cs.PHASERS:
            rows = [("indel" if (c, a) in ind or (c, b) in ind else stratum(c, a, b),
                     (c, a // cs.UNIT), key, wrong, c, a, b)
                    for c, a, key, wrong, b, _ in cs.load_rows(common.PLATFORM, depth, phaser,
                                                               chroms)]
            print(f"\n{depth}x {phaser}: {len(rows)} pairs; per group, class: share of pairs, "
                  f"error against {truth_name}", flush=True)
            for g in wc.GROUPS:
                mine = [r for r in rows if r[0] == g]
                parts = []
                for cls in ("eligible", "uninformative", "other"):
                    have = [r for r in mine
                            if pair_class(reason[(r[4], r[5])], reason[(r[4], r[6])]) == cls]
                    parts.append(f"{cls} {len(have) / max(1, len(mine)):.3f} "
                                 f"{sum(r[3] for r in have) / max(1, len(have)):.2e}")
                print(f"  {g:<10} " + "; ".join(parts))
            for name in FILTERS:
                lab = {}
                for site, (a, b) in labels.items():
                    kept = wn.FILTERS[name](a, b)
                    if kept:
                        lab[site] = kept[0]
                eligible = [relation(lab, *r[4:7]) is not None for r in rows]
                cal = [(r[0], r[1], r[2], r[3] != (relation(lab, *r[4:7])
                                                   != relation(truth, *r[4:7])))
                       for r, e in zip(rows, eligible) if e and r[4] in odd]
                test_all = [r[:4] for r in rows if r[4] not in odd]
                test_eligible = [r[:4] for r, e in zip(rows, eligible) if e and r[4] not in odd]
                for beta in cs.BETAS:
                    t = learn_then_test_by_group(cal, beta, cs.DELTA)
                    for scope, test in (("eligible", test_eligible), ("all", test_all)):
                        for g, c in ss.cells(t, test, beta).items():
                            summary.append((name, scope, depth, phaser, beta, g, c))
    for scope in ("eligible", "all"):
        mine_cells = [s for s in summary if s[1] == scope]
        level = 0.05 / max(1, len(mine_cells))
        print(f"\nscored on {scope} even-chromosome pairs; violation level {level:.2e} over "
              f"{len(mine_cells)} cells; per filter, beta, group, pooled over depths and phasers: "
              "cells certifying, violating, min p, trusted share of the scored pairs, wrong "
              "units/units")
        for name in FILTERS:
            for beta in cs.BETAS:
                for g in wc.GROUPS:
                    mine = [s[6] for s in mine_cells if s[0] == name and s[4] == beta
                            and s[5] == g]
                    if not mine:
                        continue
                    cert = [c for c in mine if c[7]]
                    bad = [c for c in cert if c[0] and c[2] < level]
                    n, k = sum(c[0] for c in mine), sum(c[1] for c in mine)
                    share = sum(c[4] for c in mine) / max(1, sum(c[3] for c in mine))
                    minp = min((c[2] for c in cert), default=1.0)
                    print(f"  {name:<7} {beta:<7g} {g:<10} {len(cert):>2} of {len(mine):<2} "
                          f"{len(bad):>3}  {minp:.2e}  {share:.4f}  {k}/{n} = "
                          f"{k / max(1, n):.2e}")
    # Gate, `agree` filter.
    level = 0.05 / max(1, sum(1 for s in summary if s[1] == "eligible"))
    eligible_cells = [s for s in summary if s[0] == "agree" and s[1] == "eligible"]
    bad = [s for s in eligible_cells if s[6][7] and s[6][0] and s[6][2] < level]
    easy = [s for s in summary if s[0] == "agree" and s[4] == 1e-3 and s[5] == "easy"]
    easy_cert = {(s[2], s[3]) for s in easy if s[1] == "eligible" and s[6][7]}
    # The procedure abstains on ineligible pairs: trusted eligible pairs over all pairs.
    share_all = (sum(s[6][4] for s in easy if s[1] == "eligible")
                 / max(1, sum(s[6][3] for s in easy if s[1] == "all")))
    a_ok, b_ok = not bad, len(easy_cert) == len(cs.DEPTHS) * len(cs.PHASERS) and share_all >= 0.70
    for s in bad:
        print(f"  violating: agree {s[2]}x {s[3]} beta {s[4]:g} {s[5]} p {s[6][2]:.2e} "
              f"({s[6][1]}/{s[6][0]})")
    print(f"\nR3 gate (agree): (a) valid on eligible pairs {'ok' if a_ok else 'FAIL'} "
          f"({len(bad)} violating of {len(eligible_cells)}); (b) easy at 1e-3 certifies in "
          f"{len(easy_cert)} of {len(cs.DEPTHS) * len(cs.PHASERS)} cells, trusted share of all "
          f"easy pairs {share_all:.4f} {'ok' if b_ok else 'FAIL'}")
    print(f"R3 gate: {'passed' if a_ok and b_ok else 'not passed'}")
    if common.SAMPLE == "HG005":
        bad = [s for s in eligible_cells if s[6][7] and s[6][0] and s[6][2] < R3B_LEVEL]
        print(f"\nR3b (registered, level {R3B_LEVEL:.3e} over {R3B_CELLS} cells): "
              f"{len(eligible_cells)} agree cells have pairs, min p "
              f"{min((s[6][2] for s in eligible_cells if s[6][7] and s[6][0]), default=1.0):.3e}")
        for s in bad:
            print(f"  violating: {s[2]}x {s[3]} beta {s[4]:g} {s[5]} p {s[6][2]:.2e} "
                  f"({s[6][1]}/{s[6][0]})")
        print(f"R3b-1 validity: {'supported' if not bad else 'not supported'}")
        print(f"R3b-2 usefulness: easy at 1e-3 certifies in {len(easy_cert)} of "
              f"{len(cs.DEPTHS) * len(cs.PHASERS)} cells, trusted share of all even easy pairs "
              f"{share_all:.4f}: {'supported' if b_ok else 'not supported'}")


if __name__ == "__main__":
    main()
