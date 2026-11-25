"""Round 4, W5 part 2 (development): calibration on labels from called parental genotypes.

The child is HG002 (`cal_v5`, ONT sup v5, the third test's calibration
reads). Its parents' genotypes come from ONT's own `wf-human-variation`
calls (Clair3 1.0.8) on the same 2025 sup data, one flowcell per parent
(`data/w5/`, fetched 14 October; md5 in `fetched.md5`); nothing about them
uses a benchmark. A child's heterozygous site is labelled when a parent's
call is homozygous there: `0/0` when the parent has no PASS record at those
alleles, `1/1` when it has one. A pair is labelled when both its sites are,
and its label is the relation of the two paternal alleles.

For each depth and phaser it reports, per group, the share of pairs
labelled and the share of labelled pairs whose label disagrees with v5.0q,
then calibrates per-group thresholds on the labelled pairs of odd
chromosomes against the labels (no truth used) and scores them on all pairs
of the even chromosomes against v5.0q. The oracle calibrates on all odd
pairs against v5.0q. Both flowcell pairs of the parents are run, as
replicates.

Usage: `PHASEBUDGET_PLATFORM=cal_v5 PYTHONPATH=src:benchmarks .venv/bin/python
benchmarks/w5_called.py`
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import confidence_score as cs  # noqa: E402
import genome_common as common  # noqa: E402
import stratified_score as ss  # noqa: E402

from phasebudget.confidence import learn_then_test_by_group  # noqa: E402

CALLS = common.ROOT / "data" / "w5"
# (father, mother) per flowcell replicate; HG005's are R3b's (`docs/design/r3b-parent-labels-plan.md`).
REPLICATES = ({"A": ("HG006_PAY77227", "HG007_PAY12990"), "B": ("HG006_PBA16846", "HG007_PBA20413")}
              if common.SAMPLE == "HG005" else
              {"A": ("HG003_PAY87794", "HG004_PAY87778"), "B": ("HG003_PAY87954", "HG004_PAY88428")})
GROUPS = ("easy", "difficult", "lowmap", "segdup", "tandem", "indel")


def calls(name: str, chroms) -> dict:
    """`{(chrom, pos1, ref, alt): alt alleles}`: 1 or 2 for PASS calls, 1 (uninformative) for any
    other record there, so only a site with no record at all reads as `0/0`."""
    out = {}
    with gzip.open(CALLS / f"{name}.wf_snp.vcf.gz", "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.split("\t", 10)
            if f[0] not in chroms:
                continue
            gt = f[9].split(":")[0].replace("|", "/")
            for alt in f[4].split(","):
                key = (f[0], int(f[1]), f[3], alt)
                ok = f[6] == "PASS" and "," not in f[4] and gt in ("0/1", "1/0", "1/1")
                out[key] = gt.count("1") if ok else 1
    return out


def paternal(chrom: str, variants, father: dict, mother: dict) -> dict:
    """`{pos1: paternal allele}` the parents' calls imply, where they imply one."""
    out = {}
    for pos1, ref, alt in variants:
        f = father.get((chrom, pos1, ref, alt), 0)
        m = mother.get((chrom, pos1, ref, alt), 0)
        if f != 1 and m != 1:  # both homozygous: only a Mendelian-consistent pair labels
            if f != m:
                out[pos1] = f // 2
        elif f != 1:
            out[pos1] = f // 2
        elif m != 1:
            out[pos1] = 1 - m // 2
    return out


def main() -> None:
    from stratified_gate import strata

    sc = common._module("score_hg002.py")
    chroms = list(common.CONFIRMATORY)
    odd = {c for c in chroms if int(c[3:]) % 2}
    stratum = strata(chroms)
    ind = ss.indels(common.PLATFORM, chroms)
    print(f"W5 part 2 (development): child HG002 ({common.PLATFORM}); parents' labels from ONT "
          "wf-human-variation sup calls (Clair3 1.0.8), PASS only, a missing record read as 0/0; "
          f"calibrate on labelled odd-chromosome pairs against the labels, score all "
          f"even-chromosome pairs against v5.0q; {cs.UNIT // 1000} kb units, delta {cs.DELTA}",
          flush=True)
    truth_pat, label_pat = {}, {name: {} for name in REPLICATES}
    parents = {name: (calls(fa, set(chroms)), calls(mo, set(chroms)))
               for name, (fa, mo) in REPLICATES.items()}
    for chrom in chroms:
        variants = sc.variant_table(common.hetsnv_path(chrom))
        for v, t in zip(variants, common.v5_truth(chrom, variants)):
            if t is not None:
                truth_pat[(chrom, v[0])] = t
        for name, (father, mother) in parents.items():
            for pos1, a in paternal(chrom, variants, father, mother).items():
                label_pat[name][(chrom, pos1)] = a

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
            oracle = learn_then_test_by_group([r[:4] for r in rows if r[4] in odd], 1e-3,
                                              cs.DELTA)
            print(f"\n{depth}x {phaser}: {len(rows)} pairs")
            for name in REPLICATES:
                lab = label_pat[name]
                print(f"  replicate {name}: group  labelled_share  label_errors  "
                      "wrong_by_truth_labelled")
                cal = []
                for g in GROUPS:
                    mine = [r for r in rows if r[0] == g]
                    have = [r for r in mine if relation(lab, r[4], r[5], r[6]) is not None]
                    errs = sum(relation(lab, *r[4:7]) != relation(truth_pat, *r[4:7])
                               for r in have)
                    print(f"    {g:<10} {len(have) / max(1, len(mine)):>14.3f} "
                          f"{errs / max(1, len(have)):>13.2e} "
                          f"{sum(r[3] for r in have) / max(1, len(have)):>12.2e}")
                for r in rows:
                    rel = relation(lab, *r[4:7])
                    if r[4] in odd and rel is not None:
                        # wrong by the label: the phaser's relation differs from the label's
                        wrong_label = r[3] != (rel != relation(truth_pat, *r[4:7]))
                        cal.append((r[0], r[1], r[2], wrong_label))
                for beta in cs.BETAS:
                    t = learn_then_test_by_group(cal, beta, cs.DELTA)
                    for g, c in ss.cells(t, test, beta).items():
                        summary.append((depth, phaser, name, beta, g, c))
                    if beta == 1e-3:
                        for g, c in ss.cells(oracle, test, beta).items():
                            summary.append((depth, phaser, "oracle", beta, g, c))
    level = 0.05 / max(1, sum(1 for s in summary if s[2] != "oracle"))
    print(f"\nviolation level {level:.2e}; per calibration, beta and group, pooled over depths "
          "and phasers: cells certifying, violating, trusted share, wrong units/units")
    for name in ("oracle", *REPLICATES):
        for beta in cs.BETAS:
            for g in GROUPS:
                mine = [s[5] for s in summary if s[2] == name and s[3] == beta and s[4] == g]
                if not mine:
                    continue
                cert = [c for c in mine if c[7]]
                bad = [c for c in cert if c[0] and c[2] < level]
                n, k = sum(c[0] for c in mine), sum(c[1] for c in mine)
                share = sum(c[4] for c in mine) / max(1, sum(c[3] for c in mine))
                print(f"  {name:<7} {beta:<7g} {g:<10} {len(cert):>2} of {len(mine):<2} "
                      f"{len(bad):>3}  {share:.4f}  {k}/{n} = {k / max(1, n):.2e}")


if __name__ == "__main__":
    main()
