"""Decide E1 and E2 of `docs/design/confidence-run-plan.md`, and report E3 and E4.

It reads only stored `realign` fragments, phaser VCFs and truths:
- calibration: HG002 ONT sup (`PHASEBUDGET_PLATFORM=ont`), 21 confirmatory
  chromosomes, truth v5.0q;
- test: HG001 ONT sup (`hg001_ont`, Platinum truth) and HG002 ONT hac
  (`ont_hac`, v5.0q).

Each (dataset, depth, phaser) is reduced to rows `(chromosome, position,
score, wrong)` for pairs whose two SNVs the truth phases, in a subprocess
with that dataset's environment. The frozen design then applies:
- 100 kb units;
- `coverage_grid` from 50 per cent;
- `learn_then_test` at delta 0.05 for beta in {3e-4, 1e-3, 3e-3};
- the E1 violation test at 0.05/36;
- the E2 rule.

`--self-check` instead calibrates on the odd and tests on the even HG002 ONT
sup chromosomes. It reproduces the 100 kb development tables
(`2026-10-02_gate_d_ont_*_100kb.txt`) and decides nothing.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/confidence_score.py`
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

UNIT = 100_000
BETAS = (3e-4, 1e-3, 3e-3)
DELTA = 0.05
DATASETS = ("hg001_ont", "ont_hac")
DEPTHS = ("30", "20", "10")
PHASERS = ("whatshap", "longphase")
CELLS = len(DATASETS) * len(DEPTHS) * len(PHASERS) * len(BETAS)
COVERAGES = (0.95, 0.98, 0.99, 1.0)
CACHE = HERE.parent / "data" / "confidence_rows"


def rows_for(depth: str, phaser: str, chromosomes) -> list:
    """`[(chrom, pos1, score, wrong)]` in the current platform's environment."""
    import genome_common as common
    from phase_confidence import phasing_from_vcf
    from phasebudget.confidence import pair_support, score

    extractor = common._module("extract_fragments.py")
    sc = common._module("score_hg002.py")
    out = []
    for chrom in chromosomes:
        variants = sc.variant_table(common.hetsnv_path(chrom))
        matrix, positions, _, _ = extractor.load_fragments(
            str(common.fragments_path(chrom, depth, "realign")))
        truth = common.truth_at(chrom, variants, positions)
        quality = None
        if phaser == "hapcut2":
            from compare_hapcut2_quality import hapcut2_blocks

            phasing, quality = hapcut2_blocks(
                common.chrom_dir(chrom) / "external" / f"hapcut2_{depth}.blocks", variants)
        else:
            vcf = common.chrom_dir(chrom) / "external" / f"{phaser}_{depth}.vcf"
            vcf = vcf if vcf.exists() else vcf.with_name(vcf.name + ".gz")
            phasing = phasing_from_vcf(str(vcf), common.VCF_SAMPLE, chrom, positions)
        for s in pair_support(matrix, phasing):
            a, b = s.first, s.second
            if truth[a] is None or truth[b] is None:
                continue
            wrong = (phasing.alleles[a] == phasing.alleles[b]) != (truth[a] == truth[b])
            hq = None if quality is None else min(quality[a], quality[b])
            out.append((chrom, positions[a] + 1, list(score(s)), bool(wrong), positions[b] + 1, hq))
    return out


def load_rows(platform: str, depth: str, phaser: str, chromosomes) -> list:
    """Rows from the cache, or computed in a subprocess with `platform`'s environment."""
    CACHE.mkdir(parents=True, exist_ok=True)
    truth = os.environ.get("PHASEBUDGET_TRUTH")  # an alternative truth gets its own cache
    path = CACHE / f"{platform}{'_' + truth if truth else ''}_{depth}_{phaser}_{'-'.join(chromosomes)}.json"
    if not path.exists():
        env = dict(os.environ, PHASEBUDGET_PLATFORM=platform,
                   PYTHONPATH=f"{HERE.parent / 'src'}:{HERE}")
        subprocess.run([sys.executable, __file__, "--rows", depth, phaser, str(path),
                        *chromosomes], env=env, check=True)
    return [(c, p, tuple(k), w, q, h) for c, p, k, w, q, h in json.loads(path.read_text())]


def evaluate(cal: list, test: list) -> dict:
    """Thresholds from `cal`, applied to `test`, for every beta."""
    from phasebudget.confidence import Units, binomial_cdf, coverage_grid, learn_then_test

    def unit_rows(rows):
        return [((c, p // UNIT), k, w) for c, p, k, w, *_ in rows]

    cal_units, test_units = Units.from_pairs(unit_rows(cal)), Units.from_pairs(unit_rows(test))
    grid = coverage_grid([r[2] for r in cal], start=0.5, step=0.0025)
    all_wrong = sum(r[3] for r in test)
    out = {"pairs": len(test), "wrong_pairs": all_wrong, "units": len(test_units.best),
           "wrong_units": len(test_units.best_wrong), "betas": {}}
    for beta in BETAS:
        t = learn_then_test(cal_units, grid, beta, DELTA)
        if t is None:
            out["betas"][beta] = None
            continue
        n, k = test_units.at(t)
        trusted = [r for r in test if r[2] >= t]
        p = 1.0 if k == 0 else 1.0 - binomial_cdf(k - 1, n, beta)
        untrusted_wrong = [r for r in test if r[3] and r[2] < t]
        out["betas"][beta] = {"threshold": t, "units": n, "wrong_units": k, "rate": k / n,
                              "strata": strata_shares(untrusted_wrong),
                              "p_exceed": p, "trusted_share": len(trusted) / len(test),
                              "wrong_trusted": sum(r[3] for r in trusted)}
    ranked = sorted(test, key=lambda r: r[2], reverse=True)
    out["risk_coverage"] = {c: sum(r[3] for r in ranked[: int(round(c * len(ranked)))])
                            / max(1, int(round(c * len(ranked)))) for c in COVERAGES}
    return out


_BEDS = {}


def strata_shares(rows) -> dict:
    """E4: per GIAB v3.6 stratum, the wrong pairs among `rows` with both SNVs inside it."""
    import genome_common as common

    stratify = common._module("stratify.py")
    if not _BEDS:
        _BEDS.update({name: stratify.load_bed(path, common.CONFIRMATORY)
                      for name, path in stratify.fetch().items()})
    counts = {name: 0 for name in _BEDS}
    for chrom, pos1, _, _, pos2, _ in rows:
        for name, bed in _BEDS.items():
            inside = stratify.inside(bed[chrom], [pos1 - 1, pos2 - 1])
            counts[name] += inside[0] and inside[1]
    counts["total"] = len(rows)
    return counts


def report(label: str, result: dict) -> list[str]:
    lines = [f"{label}: {result['pairs']} pairs, {result['wrong_pairs']} wrong; {result['units']} "
             f"units, {result['wrong_units']} wrong (base rate "
             f"{result['wrong_units'] / max(1, result['units']):.2e})"]
    for beta, r in result["betas"].items():
        if r is None:
            lines.append(f"    beta {beta:.0e}: none certified")
            continue
        lines.append(f"    beta {beta:.0e}: trusted {r['trusted_share']:.4f} of pairs; units "
                     f"wrong {r['wrong_units']}/{r['units']} = {r['rate']:.2e}; P(X >= k) "
                     f"{r['p_exceed']:.3g}; wrong trusted pairs {r['wrong_trusted']} of "
                     f"{result['wrong_pairs']}; untrusted wrong pairs by stratum "
                     + ", ".join(f"{k} {v}" for k, v in r["strata"].items()))
    rc = ", ".join(f"{c:.2f}: {1e4 * v:.2f}" for c, v in result["risk_coverage"].items())
    lines.append(f"    risk-coverage, error per 10,000 trusted pairs at coverage {rc}")
    return lines


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", nargs="+", help=argparse.SUPPRESS)
    parser.add_argument("--self-check", action="store_true")
    parser.add_argument("--hifi-transfer", action="store_true",
                        help="exploratory, not registered: calibrate on HG002 HiFi, test on HG005")
    parser.add_argument("--hg001-giab", action="store_true",
                        help="exploratory, after scoring: HG001 restricted to GIAB v4.2.1")
    args = parser.parse_args()
    if args.rows:
        depth, phaser, path, *chromosomes = args.rows
        Path(path).write_text(json.dumps(rows_for(depth, phaser, chromosomes)))
        return
    import genome_common as common

    if args.self_check:
        print("self-check: calibrate on odd, test on even HG002 ONT sup chromosomes (development)")
        for depth in DEPTHS:
            for phaser in PHASERS:
                cal = load_rows("ont", depth, phaser, common.ODD)
                test = load_rows("ont", depth, phaser, common.EVEN)
                print("\n".join(report(f"  ont {depth}x {phaser}", evaluate(cal, test))))
        return
    if args.hifi_transfer:
        print("exploratory, not registered: calibrate on HG002 HiFi (v5.0q), test on HG005 HiFi "
              "(transmission truth); whatshap phase at default settings; both seen in PhaseBudget")
        for depth in ("full", "10"):
            cal = load_rows("hifi", depth, "whatshap", common.CONFIRMATORY)
            test = load_rows("hg005", depth, "whatshap", common.CONFIRMATORY)
            print("\n".join(report(f"  hg005 {depth}x whatshap", evaluate(cal, test))))
        return
    if args.hg001_giab:
        hg001_giab(common)
        return
    print("error-controlled phasing, registered test (docs/design/confidence-run-plan.md); "
          f"units {UNIT // 1000} kb, delta {DELTA}, E1 level 0.05/{CELLS}")
    violations, cells, e2 = 0, 0, []
    for dataset in DATASETS:
        for depth in DEPTHS:
            for phaser in PHASERS:
                cal = load_rows("ont", depth, phaser, common.CONFIRMATORY)
                test = load_rows(dataset, depth, phaser, common.CONFIRMATORY)
                result = evaluate(cal, test)
                print("\n".join(report(f"  {dataset} {depth}x {phaser}", result)))
                for beta, r in result["betas"].items():
                    if r is None:
                        continue
                    cells += 1
                    violations += r["p_exceed"] < 0.05 / CELLS
                    if beta == 1e-3 and depth in ("20", "10"):
                        ok = (r["trusted_share"] >= 0.90
                              and r["wrong_trusted"] <= 0.5 * result["wrong_pairs"])
                        e2.append((f"{dataset} {depth}x {phaser}", ok))
    print(f"E1 validity: {violations} violations in {cells} certified cells of {CELLS} -> "
          f"{'supported' if violations == 0 else 'not supported'} (rule: no cell with "
          f"P(X >= k | n, beta) below 0.05/{CELLS})")
    good = sum(ok for _, ok in e2)
    verdict = "supported" if good == 8 else "partly supported" if good >= 6 else "not supported"
    print(f"E2 usefulness: {good} of {len(e2)} cells succeed -> {verdict} (rule: beta 1e-3, "
          "depths 20 and 10, trusted share >= 0.90 and wrong trusted pairs <= half)")
    for label, ok in e2:
        print(f"    {label}: {'succeeds' if ok else 'fails'}")
    e5(common)


def e5(common) -> None:
    """Amendment 2: on HG002 hac, the frozen score against HapCUT2's own quality."""
    wins = 0
    print("E5 (amendment 2): HapCUT2 phasings of HG002 ONT hac, wrong pairs left at a share "
          "trusted, ranked by the frozen score and by HapCUT2's quality")
    for depth in DEPTHS:
        test = load_rows("ont_hac", depth, "hapcut2", common.CONFIRMATORY)
        n = len(test)
        left = {}
        for name, key in (("frozen", lambda r: r[2]), ("hapcut2", lambda r: (r[5], 0))):
            ranked = sorted(test, key=key, reverse=True)
            left[name] = {c: sum(r[3] for r in ranked[: int(round(c * n))]) for c in (0.99, 0.995)}
        wins += left["frozen"][0.99] < left["hapcut2"][0.99]
        print(f"    {depth}x: {n} pairs, {sum(r[3] for r in test)} wrong; at 0.99 frozen "
              f"{left['frozen'][0.99]} against HapCUT2 {left['hapcut2'][0.99]}; at 0.995 frozen "
              f"{left['frozen'][0.995]} against HapCUT2 {left['hapcut2'][0.995]}")
        cal = load_rows("ont", depth, "hapcut2", common.CONFIRMATORY)
        print("\n".join("    " + line for line in
                        report(f"Learn then Test on HapCUT2, {depth}x (descriptive)",
                               evaluate(cal, test))))
    print(f"E5: frozen score fewer wrong pairs at 0.99 in {wins} of {len(DEPTHS)} depths -> "
          f"{'supported' if wins == len(DEPTHS) else 'not supported'} (rule: all three)")


GIAB_HG001 = HERE.parent / "data" / "giab_hg001" / "HG001_GRCh38_1_22_v4.2.1_benchmark"


def giab_hg001_sites(chromosomes) -> dict:
    """Per chromosome, GIAB HG001 v4.2.1 BED intervals and its het SNVs (1-based) inside them."""
    import bisect
    import pysam

    beds: dict = {c: ([], []) for c in chromosomes}
    with open(f"{GIAB_HG001}.bed") as handle:
        for line in handle:
            f = line.split("\t")
            if f[0] in beds:
                beds[f[0]][0].append(int(f[1]))
                beds[f[0]][1].append(int(f[2]))
    sites = {}
    with pysam.TabixFile(f"{GIAB_HG001}.vcf.gz") as tabix:
        for chrom in chromosomes:
            starts, ends = beds[chrom]
            keep = set()
            for line in tabix.fetch(chrom):
                f = line.split("\t", 10)
                if len(f[3]) != 1 or len(f[4]) != 1 or f[9].split(":")[0] not in ("0/1", "1/0"):
                    continue
                pos = int(f[1])
                k = bisect.bisect_right(starts, pos - 1) - 1
                if k >= 0 and pos - 1 < ends[k]:
                    keep.add(pos)
            sites[chrom] = keep
    return sites, beds


def single_site_flips(rows) -> int:
    """Variants that end one wrong pair and start the next: a lone flipped site."""
    ends = {(r[0], r[4]) for r in rows if r[3]}
    return sum((r[0], r[1]) in ends for r in rows if r[3])


def hg001_giab(common) -> None:
    """Exploratory: the registered HG001 cells with both SNVs GIAB v4.2.1 het SNVs in its BED."""
    print("exploratory, after the registered test was scored: HG001 ONT sup against Platinum, "
          "only pairs whose two SNVs are GIAB HG001 v4.2.1 het SNVs inside its benchmark BED "
          "(how the HG002 calibration variants were chosen); same calibration, units, grid, "
          f"delta {DELTA}, E1 level 0.05/{CELLS}")
    import bisect

    sites, beds = giab_hg001_sites(common.CONFIRMATORY)

    def in_bed(chrom, pos1):
        starts, ends = beds[chrom]
        k = bisect.bisect_right(starts, pos1 - 1) - 1
        return k >= 0 and pos1 - 1 < ends[k]

    for depth in DEPTHS:
        for phaser in PHASERS:
            cal = load_rows("ont", depth, phaser, common.CONFIRMATORY)
            full = load_rows("hg001_ont", depth, phaser, common.CONFIRMATORY)
            test = [r for r in full if r[1] in sites[r[0]] and r[4] in sites[r[0]]]
            result = evaluate(cal, test)
            dropped = [r for r in full if r[3] and not (r[1] in sites[r[0]] and r[4] in sites[r[0]])]
            outside = sum(not (in_bed(r[0], r[1]) and in_bed(r[0], r[4])) for r in dropped)
            print(f"  hg001_ont {depth}x {phaser}: dropped wrong pairs {len(dropped)}, {outside} with "
                  f"an SNV outside the GIAB BED, {len(dropped) - outside} inside it but not a GIAB "
                  "het SNV")
            print(f"  hg001_ont {depth}x {phaser}: kept {len(test)} of {len(full)} pairs; wrong "
                  f"pairs on lone flipped sites {single_site_flips(test)} of "
                  f"{result['wrong_pairs']} (all pairs: {single_site_flips(full)} of "
                  f"{sum(r[3] for r in full)})")
            print("\n".join(report(f"  hg001_ont {depth}x {phaser} in GIAB", result)))
            for beta, r in result["betas"].items():
                if r is not None:
                    print(f"    beta {beta:.0e}: "
                          f"{'violates' if r['p_exceed'] < 0.05 / CELLS else 'no violation'}")


if __name__ == "__main__":
    main()
