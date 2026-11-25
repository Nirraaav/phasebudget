"""Round 4, W5 part 1 (development): does calibrating on family-labelled pairs keep the bound?

A child's heterozygous site is labelled when a parent is homozygous there, so
that parent's transmitted allele is known; a pair is labelled when both of
its sites are. Only labelled pairs can be scored without a benchmark truth.
Here the labels come from the Platinum Pedigree truthset genotypes of the
parents (NA12877 and NA12878), so they carry no genotype error: this asks
only whether the labelled subset is a biased sample of all pairs.

For each child of NA12877 and NA12878 (second-test genomes, scored 13
October), depth and phaser, it reports per group the share of pairs
labelled and the error rates of labelled and unlabelled pairs, then
calibrates per-group thresholds on the odd chromosomes two ways (all pairs,
the oracle; labelled pairs only) and scores both on all pairs of the even
chromosomes against the truth.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/w5_labels.py`
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import confidence_score as cs  # noqa: E402
import genome_common as common  # noqa: E402
import stratified_score as ss  # noqa: E402

from phasebudget.confidence import learn_then_test_by_group  # noqa: E402

CHILDREN = ("NA12879", "NA12881", "NA12882", "NA12885", "NA12886")
PARENTS = ("NA12877", "NA12878")
GROUPS = ("easy", "difficult", "lowmap", "segdup", "tandem", "indel")


def labelled_sites(child: str, chroms) -> set:
    """`(chrom, pos1)` of the child's heterozygous sites where a parent is homozygous."""
    import pysam

    out = set()
    with pysam.VariantFile(str(common.PEDIGREE / "CEPH1463.GRCh38.family-truthset.ov.vcf.gz")) as v:
        for chrom in chroms:
            for r in v.fetch(chrom):
                gt = r.samples[child]["GT"]
                if None in gt or gt[0] == gt[1]:
                    continue
                for parent in PARENTS:
                    p = r.samples[parent]["GT"]
                    if None not in p and p[0] == p[1]:
                        out.add((chrom, r.pos))
                        break
    return out


def main() -> None:
    from stratified_gate import strata

    chroms = list(common.CONFIRMATORY)
    odd = {c for c in chroms if int(c[3:]) % 2}
    stratum = strata(chroms)
    print("W5 part 1 (development, not registered): calibration on family-labelled pairs; "
          f"children {CHILDREN}; labels from the truthset genotypes of {PARENTS}; calibrate on "
          f"odd, score on all pairs of even chromosomes; {cs.UNIT // 1000} kb units; delta "
          f"{cs.DELTA}")
    summary = []
    for child in CHILDREN:
        platform = f"ped_{child}"
        labels = labelled_sites(child, chroms)
        ind = ss.indels(platform, chroms)
        for depth in cs.DEPTHS:
            for phaser in cs.PHASERS:
                rows = []  # (group, unit, score, wrong, labelled, odd)
                for c, a, key, wrong, b, _ in cs.load_rows(platform, depth, phaser, chroms):
                    g = "indel" if (c, a) in ind or (c, b) in ind else stratum(c, a, b)
                    rows.append((g, (c, a // cs.UNIT), key, wrong,
                                 (c, a) in labels and (c, b) in labels, c in odd))
                print(f"\n{child} {depth}x {phaser}: {len(rows)} pairs")
                print("  group       pairs  labelled_share  err_labelled  err_unlabelled")
                for g in GROUPS:
                    mine = [r for r in rows if r[0] == g]
                    lab = [r for r in mine if r[4]]
                    unl = [r for r in mine if not r[4]]
                    print(f"  {g:<10} {len(mine):>7} {len(lab) / max(1, len(mine)):>15.3f}"
                          f" {sum(r[3] for r in lab) / max(1, len(lab)):>13.2e}"
                          f" {sum(r[3] for r in unl) / max(1, len(unl)):>15.2e}")
                cal_all = [r[:4] for r in rows if r[5]]
                cal_lab = [r[:4] for r in rows if r[5] and r[4]]
                test = [r[:4] for r in rows if not r[5]]
                test_lab = [r[:4] for r in rows if not r[5] and r[4]]
                for beta in cs.BETAS:
                    for name, cal, scored in (("all", cal_all, test), ("labelled", cal_lab, test),
                                              ("self-check", cal_lab, test_lab)):
                        cells = ss.cells(learn_then_test_by_group(cal, beta, cs.DELTA), scored,
                                         beta)
                        for g, c in sorted(cells.items()):
                            summary.append((child, depth, phaser, beta, name, g, c))
    scored = [s for s in summary if s[4] != "self-check"]
    level = 0.05 / max(1, len(scored))
    print(f"\ncells: calibration (all pairs, the oracle; labelled pairs) x child x depth x "
          f"phaser x beta x group, scored on all even-chromosome pairs; violation level "
          f"{level:.2e} (Bonferroni over {len(scored)} cells). `self-check` scores the labelled "
          "thresholds on labelled even-chromosome pairs only, as a user without a truth would")
    print("calibration beta    group       certifying  violating  min_p     trusted_share  "
          "wrong_units/units (pooled over children, depths, phasers)")
    for name in ("all", "labelled", "self-check"):
        for beta in cs.BETAS:
            for g in GROUPS:
                mine = [s[6] for s in summary if s[4] == name and s[3] == beta and s[5] == g]
                cert = [c for c in mine if c[7]]
                bad = [c for c in cert if c[0] and c[2] < level]
                pairs = sum(c[3] for c in mine)
                trusted = sum(c[4] for c in mine)
                n, k = sum(c[0] for c in mine), sum(c[1] for c in mine)
                minp = min((c[2] for c in cert), default=1.0)
                print(f"{name:<11} {beta:<7g} {g:<10} {len(cert):>4} of {len(mine):<4}"
                      f" {len(bad):>8}  {minp:<9.2e} {trusted / max(1, pairs):>13.4f}  "
                      f"{k}/{n} = {k / max(1, n):.2e}")

if __name__ == "__main__":
    main()
