"""Track C, `docs/design/truth-adjudication-plan.md`: the HG001 wrong pairs against a second truth.

Reads only the cached rows of the registered confidence run
(`data/confidence_rows/hg001_ont_*`), the 2017 Platinum Genomes NA12878 calls,
the Platinum Pedigree truthset's NA12878 column with its high-quality regions,
and the GIAB HG001 v4.2.1 BED. Prints C1 to C3.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/adjudicate_truth.py`
"""

from __future__ import annotations

import bisect
import gzip
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import confidence_score as cs  # noqa: E402

ROOT = HERE.parent
PEDIGREE = ROOT / "data" / "pedigree_truth"
TRUTHSET = PEDIGREE / "CEPH1463.GRCh38.family-truthset.ov.vcf.gz"
HQ = PEDIGREE / "hq_regions_final.bed.gz"
PLATINUM = ROOT / "data" / "platinum" / "NA12878.vcf.gz"
SAMPLE = "NA12878"
LEVEL = 0.05 / 36
MIN_ADJUDICABLE = 0.30


def clopper_pearson(x: int, m: int) -> tuple[float, float]:
    from scipy.stats import beta

    lo = 0.0 if x == 0 else float(beta.ppf(0.025, x, m - x + 1))
    hi = 1.0 if x == m else float(beta.ppf(0.975, x + 1, m - x))
    return lo, hi


def load_bed(path: Path, chromosomes) -> dict:
    raw: dict = {c: [] for c in chromosomes}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            f = line.split("\t", 3)
            if f[0] in raw:
                raw[f[0]].append((int(f[1]), int(f[2])))
    out = {}
    for chrom, intervals in raw.items():
        starts, ends = [], []
        for s, e in sorted(intervals):
            if ends and s <= ends[-1]:
                ends[-1] = max(ends[-1], e)
            else:
                starts.append(s)
                ends.append(e)
        out[chrom] = (starts, ends)
    return out


def inside(bed, chrom: str, pos1: int) -> bool:
    starts, ends = bed[chrom]
    k = bisect.bisect_right(starts, pos1 - 1) - 1
    return k >= 0 and pos1 - 1 < ends[k]


def phased_hets(path: Path, sample: str, chromosomes) -> dict:
    """`{(chrom, pos1): (ref, alt, first allele)}` for biallelic SNV `0|1` or `1|0` calls."""
    import pysam

    out = {}
    with pysam.TabixFile(str(path)) as tabix:
        column = 9 + tabix.header[-1].split("\t")[9:].index(sample)
        for chrom in chromosomes:
            for line in tabix.fetch(chrom):
                f = line.split("\t")
                gt = f[column].split(":")[0]
                if gt in ("0|1", "1|0") and len(f[3]) == 1 and len(f[4]) == 1:
                    out[(chrom, int(f[1]))] = (f[3], f[4], int(gt[0]))
    return out


def breakdown(chroms) -> None:
    """Exploratory, added after the registered output: why outside-GIAB SNVs are not adjudicable."""
    import collections

    import pysam

    old = phased_hets(PLATINUM, SAMPLE, chroms)
    new = phased_hets(TRUTHSET, SAMPLE, chroms)
    hq = load_bed(HQ, chroms)
    _, giab = cs.giab_hg001_sites(chroms)
    records = {}
    with pysam.TabixFile(str(TRUTHSET)) as tabix:
        column = 9 + tabix.header[-1].split("\t")[9:].index(SAMPLE)
        for chrom in chroms:
            for line in tabix.fetch(chrom):
                f = line.split("\t")
                records[(chrom, int(f[1]))] = f[column].split(":")[0]
    print("exploratory, after the registered output: each SNV of the wrong whatshap 30x pairs "
          "with an SNV outside GIAB, by its state in the Pedigree truth")
    counts = collections.Counter()
    for chrom, p1, _, wrong, p2, _ in cs.load_rows("hg001_ont", "30", "whatshap", chroms):
        if not wrong or (inside(giab, chrom, p1) and inside(giab, chrom, p2)):
            continue
        for pos in (p1, p2):
            if not inside(hq, chrom, pos):
                counts["outside the high-quality regions"] += 1
            elif (chrom, pos) in new and new[(chrom, pos)][:2] == old[(chrom, pos)][:2]:
                counts["adjudicable"] += 1
            elif (chrom, pos) not in records:
                counts["inside, no Pedigree record"] += 1
            else:
                counts[f"inside, Pedigree GT {records[(chrom, pos)]}"] += 1
    for key, n in counts.most_common():
        print(f"  {n} {key}")


def main() -> None:
    import argparse

    import genome_common as common

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--breakdown", action="store_true")
    chroms = common.CONFIRMATORY
    if parser.parse_args().breakdown:
        breakdown(chroms)
        return
    print("track C, docs/design/truth-adjudication-plan.md: HG001 ONT sup rows of the registered "
          f"confidence run; truths Platinum 2017 ({PLATINUM.name}) and Platinum Pedigree "
          f"({TRUTHSET.name}, {SAMPLE}, inside {HQ.name}); outside GIAB = an SNV outside "
          f"{cs.GIAB_HG001.name}.bed; C2 level 0.05/36; adjudicable floor {MIN_ADJUDICABLE}")
    old = phased_hets(PLATINUM, SAMPLE, chroms)
    new = phased_hets(TRUTHSET, SAMPLE, chroms)
    hq = load_bed(HQ, chroms)
    _, giab = cs.giab_hg001_sites(chroms)

    def adjudicable(chrom, pos):
        a, b = old.get((chrom, pos)), new.get((chrom, pos))
        return a is not None and b is not None and a[:2] == b[:2] and inside(hq, chrom, pos)

    c1, violations = None, 0
    for depth in cs.DEPTHS:
        for phaser in cs.PHASERS:
            rows = cs.load_rows("hg001_ont", depth, phaser, chroms)
            judged, test = [], []
            outside_wrong = 0
            for row in rows:
                chrom, p1, score, wrong, p2, hqual = row
                outside = not (inside(giab, chrom, p1) and inside(giab, chrom, p2))
                outside_wrong += wrong and outside
                if not (adjudicable(chrom, p1) and adjudicable(chrom, p2)):
                    continue
                same_old = old[(chrom, p1)][2] == old[(chrom, p2)][2]
                same_phaser = same_old != wrong
                same_new = new[(chrom, p1)][2] == new[(chrom, p2)][2]
                judged.append((outside, wrong, same_phaser == same_new, same_old == same_new))
                test.append((chrom, p1, score, same_phaser != same_new, p2, hqual))
            w = [j for j in judged if j[0] and j[1]]
            x, m = sum(j[2] for j in w), len(w)
            share = m / max(1, outside_wrong)
            lo, hi = clopper_pearson(x, m) if m else (0.0, 1.0)
            if share < MIN_ADJUDICABLE:
                verdict = "not decidable"
            elif lo > 0.5:
                verdict = "truth errors dominate"
            elif hi < 0.5:
                verdict = "phaser errors dominate"
            else:
                verdict = "undecided"
            print(f"  {depth}x {phaser}: wrong pairs outside GIAB {outside_wrong}, adjudicable {m} "
                  f"({share:.3f}); Pedigree agrees with the phaser at {x} of {m}"
                  + (f" = {x / m:.3f}, 95% CI {lo:.3f} to {hi:.3f}" if m else "")
                  + f" -> {verdict}")
            inside_w = [j for j in judged if not j[0] and j[1]]
            print(f"    wrong pairs inside GIAB, adjudicable {len(inside_w)}: Pedigree agrees with "
                  f"the phaser at {sum(j[2] for j in inside_w)}")
            if depth == "30" and phaser == "whatshap":
                c1 = verdict
                for side, flag in (("inside GIAB", False), ("outside GIAB", True)):
                    part = [j for j in judged if j[0] == flag]
                    print(f"    C3 {side}: truths disagree at {sum(not j[3] for j in part)} of "
                          f"{len(part)} adjudicable pairs")
            cal = cs.load_rows("ont", depth, phaser, chroms)
            result = cs.evaluate(cal, test)
            print("\n".join("    " + line for line in
                            cs.report(f"C2 {depth}x {phaser} against Pedigree", result)))
            for beta, r in result["betas"].items():
                if r is not None:
                    violations += r["p_exceed"] < LEVEL
    print(f"C1 (whatshap 30x): {c1}")
    rule = ("2017 truth outside GIAB" if violations == 0 else "the score" if violations == 18
            else "neither")
    print(f"C2: {violations} of 18 HG001 cells violate against the Pedigree truth -> "
          f"attributed to {rule}")


if __name__ == "__main__":
    main()
