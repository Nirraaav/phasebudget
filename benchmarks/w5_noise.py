"""Round 4, W5 part 3 (development): where called-genotype labels go wrong, and which labels to keep.

Part 2 (`benchmarks/w5_called.py`) labelled HG002's pairs from ONT Clair3
calls on HG003 and HG004 and found the labels about ten times noisier than
whatshap in `easy`, so calibration on them abstained. This script splits a
site's label by what it rests on (`docs/design/w5-label-noise.md`):

- `double`: both parents homozygous, opposite. The label holds if either
  parent's call is right, so it is wrong only when both are;
- `called`: one parent has a PASS `1/1` call, the other is heterozygous or
  uninformative;
- `absent`: one parent has no record (read as `0/0`), the other is
  heterozygous or uninformative. A missed heterozygous call in that parent
  makes this label wrong half the time.

and by whether the two flowcell replicates of each parent give the same
label. It reports each kind's share and its error against v5.0q, then
calibrates on the labelled odd-chromosome pairs kept by each filter and
scores every even-chromosome pair against v5.0q, as part 2 does.

Usage: `PHASEBUDGET_PLATFORM=cal_v5 PYTHONPATH=src:benchmarks .venv/bin/python
benchmarks/w5_noise.py`
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import confidence_score as cs  # noqa: E402
import genome_common as common  # noqa: E402
import stratified_score as ss  # noqa: E402
import w5_called as wc  # noqa: E402

from phasebudget.confidence import learn_then_test_by_group  # noqa: E402

FILTERS = {
    "all": lambda a, b: a,  # part 2, replicate A
    "agree": lambda a, b: a if a and b and a[0] == b[0] else None,
    "positive": lambda a, b: a if a and a[1] != "absent" else None,
    "agree_positive": lambda a, b: (a if a and b and a[0] == b[0] and "absent" not in (a[1], b[1])
                                    else None),
    "double": lambda a, b: a if a and a[1] == "double" else None,
}


def labelled(chrom: str, variants, father: dict, mother: dict) -> dict:
    """`{pos1: (paternal allele, kind)}`, by `wc.paternal`'s rule, with the kind it rests on."""
    out = {}
    for pos1, ref, alt in variants:
        f = father.get((chrom, pos1, ref, alt), 0)
        m = mother.get((chrom, pos1, ref, alt), 0)
        if f != 1 and m != 1:
            if f != m:
                out[pos1] = (f // 2, "double")
        elif f != 1:
            out[pos1] = (f // 2, "called" if f else "absent")
        elif m != 1:
            out[pos1] = (1 - m // 2, "called" if m else "absent")
    assert {p: a for p, (a, _) in out.items()} == wc.paternal(chrom, variants, father, mother)
    return out


def main() -> None:
    from stratified_gate import strata

    sc = common._module("score_hg002.py")
    chroms = list(common.CONFIRMATORY)
    odd = {c for c in chroms if int(c[3:]) % 2}
    stratum = strata(chroms)
    ind = ss.indels(common.PLATFORM, chroms)
    print(f"W5 part 3 (development): child HG002 ({common.PLATFORM}); labels from both "
          "replicates of the parents' Clair3 calls, split by kind; filters calibrate on "
          "labelled odd-chromosome pairs and score all even-chromosome pairs against v5.0q; "
          f"{cs.UNIT // 1000} kb units, delta {cs.DELTA}", flush=True)
    parents = {name: (wc.calls(fa, set(chroms)), wc.calls(mo, set(chroms)))
               for name, (fa, mo) in wc.REPLICATES.items()}
    truth, labels, kinds = {}, {}, Counter()
    for chrom in chroms:
        variants = sc.variant_table(common.hetsnv_path(chrom))
        for v, t in zip(variants, common.v5_truth(chrom, variants)):
            if t is not None:
                truth[(chrom, v[0])] = t
        per = {name: labelled(chrom, variants, *parents[name]) for name in wc.REPLICATES}
        for pos1, _, _ in variants:
            a, b = per["A"].get(pos1), per["B"].get(pos1)
            labels[(chrom, pos1)] = (a, b)
            t = truth.get((chrom, pos1))
            if a and t is not None:
                group = "indel" if (chrom, pos1) in ind else stratum(chrom, pos1, pos1)
                agree = "agree" if b and b[0] == a[0] else ("disagree" if b else "A only")
                kinds[(group, a[1], agree, a[0] != t)] += 1

    print("\nsites labelled by replicate A with a v5.0q phase: group, kind, replicate B, sites, "
          "wrong, error")
    for group in wc.GROUPS:
        for kind in ("double", "called", "absent"):
            for agree in ("agree", "disagree", "A only"):
                n = kinds[(group, kind, agree, False)] + kinds[(group, kind, agree, True)]
                k = kinds[(group, kind, agree, True)]
                if n:
                    print(f"  {group:<10} {kind:<7} {agree:<9} {n:>9} {k:>6} {k / n:.2e}")

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
            test = [r[:4] for r in rows if r[4] not in odd]
            print(f"\n{depth}x {phaser}: {len(rows)} pairs; filter, group, labelled share, "
                  "label errors, wrong by truth among labelled, among unlabelled", flush=True)
            for name, keep in FILTERS.items():
                lab = {}
                for site, (a, b) in labels.items():
                    kept = keep(a, b)
                    if kept:
                        lab[site] = kept[0]
                cal = []
                for g in wc.GROUPS:
                    mine = [r for r in rows if r[0] == g]
                    have = [r for r in mine if relation(lab, *r[4:7]) is not None]
                    errs = sum(relation(lab, *r[4:7]) != relation(truth, *r[4:7]) for r in have)
                    rest = len(mine) - len(have)
                    print(f"  {name:<15} {g:<10} {len(have) / max(1, len(mine)):.3f} "
                          f"{errs / max(1, len(have)):.2e} "
                          f"{sum(r[3] for r in have) / max(1, len(have)):.2e} "
                          f"{(sum(r[3] for r in mine) - sum(r[3] for r in have)) / max(1, rest):.2e}")
                for r in rows:
                    rel = relation(lab, *r[4:7])
                    if r[4] in odd and rel is not None:
                        cal.append((r[0], r[1], r[2], r[3] != (rel != relation(truth, *r[4:7]))))
                for beta in cs.BETAS:
                    t = learn_then_test_by_group(cal, beta, cs.DELTA)
                    for g, c in ss.cells(t, test, beta).items():
                        summary.append((name, beta, g, c))
            oracle = [r[:4] for r in rows if r[4] in odd]
            for beta in cs.BETAS:
                t = learn_then_test_by_group(oracle, beta, cs.DELTA)
                for g, c in ss.cells(t, test, beta).items():
                    summary.append(("oracle", beta, g, c))
    level = 0.05 / max(1, sum(1 for s in summary if s[0] != "oracle"))
    print(f"\nviolation level {level:.2e}; per calibration, beta and group, pooled over depths "
          "and phasers: cells certifying, violating, trusted share, wrong units/units")
    for name in ("oracle", *FILTERS):
        for beta in cs.BETAS:
            for g in wc.GROUPS:
                mine = [s[3] for s in summary if s[0] == name and s[1] == beta and s[2] == g]
                if not mine:
                    continue
                cert = [c for c in mine if c[7]]
                bad = [c for c in cert if c[0] and c[2] < level]
                n, k = sum(c[0] for c in mine), sum(c[1] for c in mine)
                share = sum(c[4] for c in mine) / max(1, sum(c[3] for c in mine))
                print(f"  {name:<15} {beta:<7g} {g:<10} {len(cert):>2} of {len(mine):<2} "
                      f"{len(bad):>3}  {share:.4f}  {k}/{n} = {k / max(1, n):.2e}")


if __name__ == "__main__":
    main()
