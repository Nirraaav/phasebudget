"""Gate A of the extension plan (development repository): per-stratum Learn then Test on HG001, development.

Reads the cached HG001 ONT sup rows of the registered confidence run and
labels each pair against the Platinum Pedigree truth, on the pairs track C
calls adjudicable (both SNVs heterozygous with one genotype in both truths,
inside the Pedigree high-quality regions). The phaser's relation is
recovered from the row's loss against the 2017 truth, as in track C.

Each pair gets the stratum of `docs/design/stratified-confidence.md`.
Thresholds come from `learn_then_test_by_group` on the odd chromosomes and are
scored on the even ones. Gate A holds when no stratum's realized unit rate
exceeds beta at p < 0.05 in any cell, and the easy stratum keeps at least
0.99 of its pairs at beta 1e-3 and 30-fold, for both phasers.

With `--pooled`, every pair is in one group "all": the single threshold of
the registered test, calibrated and scored the same way, as the baseline.

With `--cross`, the calibration is instead the registered one: HG002 ONT sup
rows on all confirmatory chromosomes, labelled by v5.0q, and every HG001
chromosome is test. This is the setting where the pooled threshold failed
(track C, C2). Violations are also counted at the registered level 0.05/36.

The data were seen before: this is development, not a test.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/stratified_gate.py`
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import adjudicate_truth as adj  # noqa: E402
import confidence_score as cs  # noqa: E402
import genome_common as common  # noqa: E402

from phasebudget.confidence import Units, binomial_cdf, learn_then_test_by_group  # noqa: E402

ORDER = ("segdup", "lowmap", "tandem", "difficult", "easy")
GATE_P = 0.05
EASY_SHARE = 0.99


def strata(chroms):
    """`stratum(chrom, pos1, pos2)` of `docs/design/stratified-confidence.md`, 1-based positions."""
    stratify = common._module("stratify.py")
    beds = {name: stratify.load_bed(path, chroms) for name, path in stratify.fetch().items()}
    named = {"segdup": beds["segdups"], "lowmap": beds["lowmap"], "tandem": beds["TR+homopol"]}
    easy = beds["easy"]

    def stratum(chrom, p1, p2):
        for name, bed in named.items():
            if any(stratify.inside(bed[chrom], [p1 - 1, p2 - 1])):
                return name
        return "easy" if all(stratify.inside(easy[chrom], [p1 - 1, p2 - 1])) else "difficult"

    return stratum


def snvs(chroms) -> set:
    """`(chrom, pos1)` of the SNVs in track B's `cal_r10` input."""
    import gzip

    out = set()
    for chrom in chroms:
        with gzip.open(common.ROOT / "data" / "ext_cal_r10" / chrom / "hetsnv.vcf.gz", "rt") as h:
            for line in h:
                f = line.split("\t", 5)
                if not line.startswith("#") and len(f[3]) == 1 and len(f[4]) == 1:
                    out.add((chrom, int(f[1])))
    return out


def hg001_pedigree(chroms):
    """`rows(depth, phaser)`: HG001 rows on track C's adjudicable pairs as `(chrom, p1, score,
    wrong, p2)`, wrong against the Platinum Pedigree truth; the phaser's relation is recovered
    from the row's loss against the 2017 truth."""
    old = adj.phased_hets(adj.PLATINUM, adj.SAMPLE, chroms)
    new = adj.phased_hets(adj.TRUTHSET, adj.SAMPLE, chroms)
    hq = adj.load_bed(adj.HQ, chroms)

    def judged(chrom, pos):
        a, b = old.get((chrom, pos)), new.get((chrom, pos))
        return a is not None and b is not None and a[:2] == b[:2] and adj.inside(hq, chrom, pos)

    def rows(depth, phaser):
        out = []
        for chrom, p1, key, wrong, p2, _ in cs.load_rows("hg001_ont", depth, phaser, chroms):
            if judged(chrom, p1) and judged(chrom, p2):
                same_phaser = (old[(chrom, p1)][2] == old[(chrom, p2)][2]) != wrong
                bad = same_phaser != (new[(chrom, p1)][2] == new[(chrom, p2)][2])
                out.append((chrom, p1, key, bad, p2))
        return out

    return rows


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pooled", action="store_true")
    parser.add_argument("--cross", action="store_true")
    parser.add_argument("--calibration", default="ont", choices=("ont", "cal_r10"),
                        help="with --cross: the registered rows, or track B's SNV pairs")
    args = parser.parse_args()
    pooled, cross = args.pooled, args.cross
    order = ("all",) if pooled else ORDER
    chroms = common.CONFIRMATORY
    cal_chroms = [] if cross else [c for c in chroms if int(c[3:]) % 2]
    setting = (f"calibrate on HG002 ONT sup rows ({args.calibration}, v5.0q) on all of them, "
               "test on HG001 on all of them" if cross else
               f"calibrate on HG001 {' '.join(cal_chroms)}, test on HG001 on the rest of them")
    print("gate A, the extension plan (development repository) section 1 (development; data seen): HG001 ONT "
          f"sup rows of the registered confidence run, truth Platinum Pedigree {adj.SAMPLE} on "
          f"track C's adjudicable pairs; {len(chroms)} confirmatory chromosomes, {setting}; groups {', '.join(order)} (GIAB v3.6, first "
          f"that holds for either SNV); units {cs.UNIT // 1000} kb per stratum; betas "
          f"{cs.BETAS}; delta {cs.DELTA}; grid from 0.5 in steps of 0.0025; gate p < {GATE_P} "
          f"per group and cell, {order[-1]} keeps >= {EASY_SHARE} at beta 1e-3 and 30x")
    stratum = strata(chroms)
    hg001 = hg001_pedigree(chroms)
    failures, registered, easy_shares = [], [], []
    for depth in cs.DEPTHS:
        for phaser in cs.PHASERS:
            cal, test = [], []
            for chrom, p1, key, bad, p2 in hg001(depth, phaser):
                s = "all" if pooled else stratum(chrom, p1, p2)
                (cal if chrom in cal_chroms else test).append((s, (chrom, p1 // cs.UNIT), key, bad))
            if cross:
                snv = snvs(chroms) if args.calibration == "cal_r10" else None
                cal = [("all" if pooled else stratum(c, p1, p2), (c, p1 // cs.UNIT), key, wrong)
                       for c, p1, key, wrong, p2, _ in cs.load_rows(args.calibration, depth,
                                                                     phaser, chroms)
                       if snv is None or ((c, p1) in snv and (c, p2) in snv)]
            print(f"\n{depth}x {phaser}: calibration pairs {len(cal)}, test pairs {len(test)}")
            print("  stratum     beta   cal pairs (wrong)  test pairs (wrong)  threshold"
                  "             test units  wrong  rate      p_exceed  trusted share  wrong trusted")
            for beta in cs.BETAS:
                chosen = learn_then_test_by_group(cal, beta, cs.DELTA)
                for s in order:
                    c = [r for r in cal if r[0] == s]
                    t = [r for r in test if r[0] == s]
                    head = (f"  {s:<10}  {beta:<6g} {len(c):>9} ({sum(r[3] for r in c):>4})"
                            f"  {len(t):>10} ({sum(r[3] for r in t):>4})")
                    threshold = chosen.get(s)
                    if threshold is None or not t:
                        print(f"{head}  {'none' if threshold is None else threshold}")
                        continue
                    n, k = Units.from_pairs([(u, key, bad) for _, u, key, bad in t]).at(threshold)
                    p = 1.0 if k == 0 else 1.0 - binomial_cdf(k - 1, n, beta)
                    trusted = [r for r in t if r[2] >= threshold]
                    share = len(trusted) / len(t)
                    print(f"{head}  ({threshold[0]:.6f}, {threshold[1]:>3})  {n:>10}  {k:>5}  "
                          f"{(k / n if n else 0):<8.2e}  {p:<8.2e}  {share:>13.4f}  "
                          f"{sum(r[3] for r in trusted):>13}")
                    if p < 0.05 / 36:
                        registered.append(f"{depth}x {phaser} {s} beta {beta:g}")
                    if p < GATE_P:
                        failures.append(f"{depth}x {phaser} {s} beta {beta:g} (p {p:.2e})")
                    if s == order[-1] and beta == 1e-3 and depth == "30":
                        easy_shares.append(share)
    print(f"\nviolations at p < {GATE_P}: {len(failures)}" + "".join(f"\n  {f}" for f in failures))
    print(f"violations at the registered level 0.05/36: {len(registered)}"
          + "".join(f"\n  {f}" for f in registered))
    keeps = len(easy_shares) == len(cs.PHASERS) and min(easy_shares) >= EASY_SHARE
    print(f"{order[-1]} share trusted at beta 1e-3, 30x: {easy_shares} -> "
          f"{'kept' if keeps else 'not kept'}")
    print(f"gate A: {'passed' if not failures and keeps else 'failed'}")


if __name__ == "__main__":
    main()
