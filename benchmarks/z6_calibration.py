"""Z6 development comparison, `docs/design/z6-grandparent-genotypes.md`: NA12878 calibrated on labels.

Transmission labels for NA12878's heterozygous SNVs come from NA12891 and
NA12892 genotyped by `grandparent_genotypes.py`. The `ped_NA12878` SNV pairs
are calibrated on the odd chromosomes against the Platinum Pedigree truth and
against the labels, and both thresholds are scored on the even chromosomes
against the truth. No test genome is read.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/z6_calibration.py`
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import adjudicate_truth as adj  # noqa: E402
import confidence_score as cs  # noqa: E402
import genome_common as common  # noqa: E402
from family_calibration import score  # noqa: E402
from grandparent_genotypes import OUT, genotype  # noqa: E402

from phasebudget.confidence import learn_then_test_by_group  # noqa: E402
from phasebudget.truth import Call, Regions, phase_by_transmission  # noqa: E402

PARENTS = ("NA12891", "NA12892")


def labels(chrom: str, rule: str = "R0") -> dict:
    """`{pos1: first-parent allele}` for NA12878 on `chrom`, by transmission from read genotypes."""
    sites = {}
    with gzip.open(OUT / "sites" / f"{chrom}.vcf.gz", "rt") as handle:
        for line in handle:
            if not line.startswith("#"):
                f = line.split("\t", 6)
                sites[int(f[1])] = (f[3], f[4])
    child = {p: Call(p, r, (a,), (0, 1), False, None) for p, (r, a) in sites.items()}
    parents, holes = [], []
    for parent in PARENTS:
        calls, missing = {}, set()
        for path in sorted((OUT / parent).glob(f"{chrom}_*.tsv.gz")):
            with gzip.open(path, "rt") as handle:
                for line in handle:
                    p, ref, alt = map(int, line.split())
                    g = genotype(ref, alt, rule)
                    if g is None:
                        missing.add(p)
                    elif g != (0, 0):
                        calls[p] = Call(p, sites[p][0], (sites[p][1],), g, False, None)
        parents.append(calls)
        edges = [0] + [x for p in sorted(missing) for x in (p - 1, p)] + [10 ** 10]
        holes.append(Regions([(a, b) for a, b in zip(edges[::2], edges[1::2]) if a < b]))
    phased, _ = phase_by_transmission(child, parents[0], parents[1], holes[0], holes[1])
    return {p: pat for p, (pat, _) in phased.items()}


def disagreements(lab: dict, truth: dict, chroms) -> tuple[int, int, int]:
    """`(labelled SNVs, consecutive labelled relations, those disagreeing with the truth)`."""
    n = disagree = 0
    for c in chroms:
        sites = sorted(p for p in lab[c] if (c, p) in truth)
        for a, b in zip(sites, sites[1:]):
            n += 1
            disagree += (lab[c][a] == lab[c][b]) != (truth[(c, a)][2] == truth[(c, b)][2])
    return sum(len(lab[c]) for c in chroms), n, disagree


def main() -> None:
    import argparse

    from grandparent_genotypes import RULES
    from stratified_gate import strata

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--choose", action="store_true",
                        help="the repair: choose the rule on the odd chromosomes first")
    args = parser.parse_args()
    chroms = list(common.CONFIRMATORY)
    odd = {c for c in chroms if int(c[3:]) % 2}
    truth = adj.phased_hets(adj.TRUTHSET, adj.SAMPLE, chroms)
    rule = "R0"
    if args.choose:
        odd_list = sorted(odd, key=lambda c: int(c[3:]))
        found = {r: disagreements({c: labels(c, r) for c in odd_list}, truth, odd_list)
                 for r in RULES}
        for r, (sites, n, k) in found.items():
            print(f"rule {r} {RULES[r]} on odd chromosomes: labelled {sites}, relations {n}, "
                  f"disagree {k} ({k / max(1, n):.2e})")
        allowed = [r for r in RULES if found[r][0] >= 0.8 * found["R0"][0]]
        rule = min(allowed, key=lambda r: (found[r][2], list(RULES).index(r)))
        print(f"chosen rule: {rule}")
    print(f"rule {rule} {RULES[rule]}")
    print("Z6 development comparison (NA12878 only): labels by transmission from NA12891 and "
          "NA12892 read genotypes (grandparent_genotypes.py rule); ped_NA12878 SNV pairs; "
          "calibrate on odd chromosomes (pairs with both SNVs labelled) against truth or labels, "
          "score on all even-chromosome SNV pairs against the Platinum Pedigree truth; pooled and "
          f"per stratum; betas {cs.BETAS}; delta {cs.DELTA}; violation p < 0.05 per group")
    lab = {c: labels(c, rule) for c in chroms}
    sites, n, disagree = disagreements(lab, truth, chroms)
    print(f"labelled SNVs {sites}; consecutive labelled relations "
          f"{n}, disagreeing with the truth {disagree} ({disagree / max(1, n):.2e})")
    stratum = strata(chroms)
    failures, shares = [], {}
    for depth in cs.DEPTHS:
        for phaser in cs.PHASERS:
            cal = {"truth": [], "labels-called": []}
            test = []
            hidden = invented = 0
            for c, p1, key, wrong, p2, _ in cs.load_rows("ped_NA12878", depth, phaser, chroms):
                if (c, p1) not in truth or (c, p2) not in truth:
                    continue  # an indel at either end, or outside the truth's phased SNVs
                row = (c, p1, p2, key)
                if c not in odd:
                    test.append(row + (wrong,))
                    continue
                if p1 not in lab[c] or p2 not in lab[c]:
                    continue
                same_phaser = (truth[(c, p1)][2] == truth[(c, p2)][2]) != wrong
                label_loss = same_phaser != (lab[c][p1] == lab[c][p2])
                hidden += wrong and not label_loss
                invented += label_loss and not wrong
                cal["truth"].append(row + (wrong,))
                cal["labels-called"].append(row + (label_loss,))
            print(f"\n{depth}x {phaser}: calibration pairs {len(cal['truth'])}, test pairs "
                  f"{len(test)}; labels hide {hidden}, invent {invented}")
            print("  mode        beta    source          viol  share   wrong trusted")
            for mode in ("pooled", "stratified"):
                def grouped(rows):
                    return [("all" if mode == "pooled" else stratum(c, a, b), (c, a // cs.UNIT),
                             k, loss) for c, a, b, k, loss in rows]

                for beta in cs.BETAS:
                    res = {}
                    for src, rows in cal.items():
                        t = learn_then_test_by_group(grouped(rows), beta, cs.DELTA)
                        res[src] = score(t, grouped(test), beta)
                        print(f"  {mode:<10}  {beta:<6g}  {src:<14}  {res[src][0]:>4}  "
                              f"{res[src][1]:.4f}  {res[src][2]:>6}")
                    if res["labels-called"][0] > res["truth"][0]:
                        failures.append(f"{depth}x {phaser} {mode} {beta:g}")
                    if beta == 1e-3 and depth == "30":
                        shares[(phaser, mode)] = (res["labels-called"][1], res["truth"][1])
    print(f"\n(a) cells where the labels violate in more groups than the truth: {len(failures)}"
          + "".join(f"\n  {f}" for f in failures))
    close = all(abs(a - b) <= 0.01 for a, b in shares.values())
    print("(b) trusted share at beta 1e-3, 30x, labels against truth: "
          + ", ".join(f"{p} {m} {a:.4f} vs {b:.4f}" for (p, m), (a, b) in shares.items()))
    print(f"Z6 development comparison: {'pass' if not failures and close else 'fail'}")


if __name__ == "__main__":
    main()
