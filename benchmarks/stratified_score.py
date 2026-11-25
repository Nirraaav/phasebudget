"""Decide S1 to S4 of `docs/design/stratified-run-plan.md` and report S5.

Per arm, depth and phaser, thresholds come from `learn_then_test_by_group` on
the calibration genome's pairs, grouped as the plan states (`indel` for a pair
with an indel at either end, else the stratum of
`docs/design/stratified-confidence.md`); each test genome's pairs are scored
against its truth. A cell is (test genome, depth, phaser, beta, group).

With `--dry-run` no test genome is read: each arm is calibrated on its
calibration genome's odd chromosomes and scored on its even ones, standing in
for every test genome. That run checks the code, not the hypotheses.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/stratified_score.py [--dry-run]`
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import confidence_score as cs  # noqa: E402
import genome_common as common  # noqa: E402

from phasebudget.confidence import Units, binomial_cdf, learn_then_test_by_group  # noqa: E402

PEDIGREE_TESTS = ("NA12877", "NA12879", "NA12881", "NA12882", "NA12885", "NA12886")
ARMS = {"R10": ("cal_r10", ("hg005_ont",)),
        "Pedigree": ("ped_NA12878", tuple(f"ped_{s}" for s in PEDIGREE_TESTS))}
DIFFICULT = ("segdup", "lowmap", "tandem", "difficult", "indel")
S4_DEPTHS, S4_BETA, S4_SHARE = ("30", "20"), 1e-3, 0.99
S6_MIN_RELATIONS = 100
S8_DEPTHS, S8_MARGIN = ("30", "20"), 0.0


def indels(platform: str, chroms) -> set:
    """`(chrom, pos1)` of the indels in `platform`'s phaser input."""
    out = set()
    for chrom in chroms:
        path = common.ROOT / "data" / f"ext_{platform}" / chrom / "hetsnv.vcf.gz"
        with gzip.open(path, "rt") as handle:
            for line in handle:
                f = line.split("\t", 5)
                if not line.startswith("#") and (len(f[3]) != 1 or len(f[4]) != 1):
                    out.add((chrom, int(f[1])))
    return out


def grouped(platform, depth, phaser, chroms, stratum) -> list:
    """`[(group, unit, score, wrong)]` for `platform`'s consecutive pairs on `chroms`."""
    ind = indels(platform, chroms)
    return [("indel" if (c, a) in ind or (c, b) in ind else stratum(c, a, b), (c, a // cs.UNIT),
             key, wrong)
            for c, a, key, wrong, b, _ in cs.load_rows(platform, depth, phaser, chroms)]


def cells(thresholds: dict, rows: list, beta: float) -> dict:
    """Per group: `(units, wrong units, p_exceed, pairs, trusted, wrong, wrong trusted,
    certified)`."""
    out = {}
    for group in {r[0] for r in rows}:
        mine = [r for r in rows if r[0] == group]
        t = thresholds.get(group)
        wrong = sum(r[3] for r in mine)
        if t is None:
            out[group] = (0, 0, 1.0, len(mine), 0, wrong, 0, False)
            continue
        n, k = Units.from_pairs([(u, key, bad) for _, u, key, bad in mine]).at(t)
        p = 1.0 if k == 0 else 1.0 - binomial_cdf(k - 1, n, beta)
        kept = [r for r in mine if r[2] >= t]
        out[group] = (n, k, p, len(mine), len(kept), wrong, sum(r[3] for r in kept), True)
    return out


def fewer_errors(c) -> bool:
    """S4's error part: half the wrong pairs removed, or, with under 20 wrong pairs to remove,
    under 1 wrong trusted pair per 10,000 trusted (the dry run of 8 October found groups with
    one or two wrong pairs, which halving cannot judge)."""
    if c[5] >= 20:
        return c[6] <= c[5] / 2
    return c[6] < 1e-4 * max(1, c[4])


def s6(thresholds: dict, tests: dict, chroms, test_chroms, stratum) -> list:
    """S6: per (test genome, phaser) at 30x and beta 1e-3, the one-sided Mann-Whitney p that
    duplicated clinical genes have lower trusted shares than the other clinical genes."""
    import bisect

    import chain_relations as cr
    from cmrg_withholding import duplicated
    from scipy.stats import mannwhitneyu

    clinical = cr.clinical_symbols()
    dup = duplicated()
    gene_map = cr.genes(test_chroms, clinical)
    out = []
    for (arm, phaser), t in thresholds.items():
        for test in tests[arm]:
            chains = cr.load_chains(test, phaser, chroms)
            ind = indels(test, test_chroms)
            shares = {}
            for chrom in test_chroms:
                rows = chains[chrom]

                def trusted(a, b, chrom=chrom):
                    g = ("indel" if (chrom, a[0]) in ind or (chrom, b[0]) in ind
                         else stratum(chrom, a[0], b[0]))
                    return t.get(g) is not None and tuple(a[3]) >= t[g]

                pos = [r[0] for r in rows]
                for start, end, name in gene_map[chrom]:
                    lo, hi = bisect.bisect_left(pos, start), bisect.bisect_right(pos, end)
                    if hi - lo >= 2:
                        every, _ = cr.count(rows[lo:hi], trusted)
                        n, tr = shares.get(name, (0, 0))
                        shares[name] = (n + every[0], tr + every[2])
            share = {g: tr / n for g, (n, tr) in shares.items() if n >= S6_MIN_RELATIONS}
            inside = [v for g, v in share.items() if g in dup]
            outside = [v for g, v in share.items() if g not in dup]
            p = mannwhitneyu(inside, outside, alternative="less").pvalue
            out.append((test, phaser, len(inside), len(outside), p))
    return out


def s7(tests, test_chroms) -> list:
    """S7: per Pedigree test genome at 30x, on indel pairs both phasers phase, the one-sided
    binomial p that pairs only whatshap gets wrong outnumber pairs only LongPhase gets wrong."""
    from scipy.stats import binomtest

    out = []
    for test in tests:
        ind = indels(test, test_chroms)
        rows = {ph: {(c, a, b): w for c, a, k, w, b, _ in cs.load_rows(test, "30", ph, test_chroms)
                     if (c, a) in ind or (c, b) in ind} for ph in cs.PHASERS}
        shared = rows["whatshap"].keys() & rows["longphase"].keys()
        only_w = sum(rows["whatshap"][k] and not rows["longphase"][k] for k in shared)
        only_l = sum(rows["longphase"][k] and not rows["whatshap"][k] for k in shared)
        p = binomtest(only_w, only_w + only_l, 0.5, alternative="greater").pvalue \
            if only_w + only_l else 1.0
        out.append((test, len(shared), only_w, only_l, p))
    return out


def main() -> None:
    import argparse

    from stratified_gate import strata

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--arms", nargs="+", default=list(ARMS), choices=list(ARMS))
    args = parser.parse_args()
    chroms = list(common.CONFIRMATORY)
    odd = [c for c in chroms if int(c[3:]) % 2]
    even = [c for c in chroms if int(c[3:]) % 2 == 0]
    stratum = strata(chroms)
    print(f"stratified run scoring{' (DRY RUN: calibration genomes only)' if args.dry_run else ''}"
          f"; units {cs.UNIT // 1000} kb per group; betas {cs.BETAS}; delta {cs.DELTA}; depths "
          f"{cs.DEPTHS}; phasers {cs.PHASERS}")
    results = []  # (arm, test, depth, phaser, beta, group, cell tuple)
    at30 = {}  # (arm, phaser): per-group thresholds at 30x and beta 1e-3, for S6
    tested = {}
    s8 = []  # (arm, test, depth, phaser, stratified share, pooled share, stratified violations)
    for arm, (cal, tests) in ((a, ARMS[a]) for a in args.arms):
        cal_chroms, test_chroms = (odd, even) if args.dry_run else (chroms, chroms)
        tests = (cal,) if args.dry_run else tests
        for depth in cs.DEPTHS:
            for phaser in cs.PHASERS:
                cal_rows = grouped(cal, depth, phaser, cal_chroms, stratum)
                test_rows = {t: grouped(t, depth, phaser, test_chroms, stratum) for t in tests}
                for beta in cs.BETAS:
                    thresholds = learn_then_test_by_group(cal_rows, beta, cs.DELTA)
                    if depth == "30" and beta == S4_BETA:
                        at30[(arm, phaser)] = thresholds
                        tested[arm] = tests
                    if depth in S8_DEPTHS and beta == S4_BETA:
                        pooled = learn_then_test_by_group([("all",) + r[1:] for r in cal_rows],
                                                          beta, cs.DELTA)
                        for t, rows in test_rows.items():
                            per = cells(thresholds, rows, beta)
                            one = cells(pooled, [("all",) + r[1:] for r in rows], beta)["all"]
                            share = sum(c[4] for c in per.values()) / sum(c[3] for c in per.values())
                            bad = sum(1 for c in per.values() if c[0] and c[2] < 0.05)
                            s8.append((arm, t, depth, phaser, share, one[4] / max(1, one[3]), bad))
                    for t, rows in test_rows.items():
                        for group, cell in cells(thresholds, rows, beta).items():
                            results.append((arm, t, depth, phaser, beta, group, cell))
    print("\narm       test          depth phaser     beta    group      units  wrong  p_exceed"
          "  pairs     trusted  wrong  wrong_trusted")
    for arm, t, depth, phaser, beta, group, c in results:
        print(f"{arm:<9} {t:<13} {depth:>5} {phaser:<10} {beta:<7g} {group:<10} {c[0]:>6} {c[1]:>6}"
              f"  {c[2]:<8.2e}  {c[3]:>8} {c[4]:>8} {c[5]:>6} {c[6]:>14}")

    def decide(name, chosen, level):
        bad = [r for r in chosen if r[6][0] and r[6][2] < level]
        print(f"{name}: {len(chosen)} cells, level {level:.2e}, violating {len(bad)}"
              + "".join(f"\n  {r[0]} {r[1]} {r[2]}x {r[3]} beta {r[4]:g} {r[5]} p {r[6][2]:.2e}"
                        for r in bad))
        return bad

    s1 = [r for r in results if r[5] == "easy"]
    s2 = [r for r in results if r[0] == "Pedigree" and r[5] in DIFFICULT]
    s3 = [r for r in results if r[0] == "R10" and r[5] in DIFFICULT]
    print()
    bad1 = decide("S1 (easy, all test genomes)", s1, 0.05 / max(1, len(s1)))
    bad2 = decide("S2 (difficult groups, Pedigree arm)", s2, 0.05 / max(1, len(s2)))
    bad3 = decide("S3 (difficult groups, R10 arm, across a basecaller change)", s3,
                  0.05 / max(1, len(s3)))
    s4 = [r for r in s1 if r[2] in S4_DEPTHS and r[4] == S4_BETA]
    ok4 = [r for r in s4 if r[6][3] and r[6][4] / r[6][3] >= S4_SHARE and fewer_errors(r[6])]
    print(f"S4 (easy, beta {S4_BETA:g}, depths {S4_DEPTHS}): {len(ok4)} of {len(s4)} cells succeed")
    abstain = {}
    for arm, t, depth, phaser, beta, group, c in results:
        key = (arm, t, depth, phaser, beta)
        a, n = abstain.get(key, (0, 0))
        abstain[key] = (a + (0 if c[7] else c[3]), n + c[3])
    print("S5 abstention (share of pairs in groups certifying nothing), beta 1e-3: "
          + "; ".join(f"{k[1]} {k[2]}x {k[3]} {a / n:.3f}" for k, (a, n) in abstain.items()
                      if k[4] == 1e-3))
    r6 = s6(at30, tested, chroms, even if args.dry_run else chroms, stratum)
    level6 = 0.05 / max(1, len(r6))
    print(f"S6 (duplicated clinical genes withhold more; level {level6:.2e}): "
          + "; ".join(f"{t} {ph} {a} vs {b} genes p {p:.2e}" for t, ph, a, b, p in r6))
    ped = tested.get("Pedigree", ())
    r7 = s7(ped, even if args.dry_run else chroms)
    level7 = 0.05 / max(1, len(r7))
    print(f"S7 (whatshap worse on shared indel pairs; level {level7:.2e}): "
          + "; ".join(f"{t} shared {n} only-whatshap {w} only-LongPhase {lp} p {p:.2e}"
                      for t, n, w, lp, p in r7))
    ok8 = [r for r in s8 if r[4] > r[5] + S8_MARGIN]
    print("S8 (per-group trusts more than pooled, beta 1e-3): "
          + "; ".join(f"{a} {t} {d}x {ph} {x:.4f} vs {y:.4f}" for a, t, d, ph, x, y, _ in s8)
          + f"\nS8 {len(ok8)} of {len(s8)} cells: "
          + ("supported" if s8 and len(ok8) == len(s8) else
             "partly supported" if s8 and len(ok8) >= 0.75 * len(s8) else "not supported"))
    sig6 = sum(p < level6 for *_, p in r6)
    sig7 = sum(p < level7 for *_, p in r7)
    print(f"S6 {sig6} of {len(r6)} significant: "
          + ("supported" if r6 and sig6 == len(r6) else
             "partly supported" if sig6 >= len(r6) - 4 and r6 else "not supported")
          + f"; S7 {sig7} of {len(r7)} significant: "
          + ("supported" if r7 and sig7 == len(r7) else "not supported"))
    print(f"\nS1 {'supported' if not bad1 else 'not supported'}; S2 "
          f"{'supported' if not bad2 else 'not supported'}; S3 "
          f"{'supported' if not bad3 else 'not supported'}; S4 "
          + ("supported" if len(ok4) == len(s4) else
             "partly supported" if len(ok4) >= len(s4) - 3 else "not supported"))


if __name__ == "__main__":
    main()
