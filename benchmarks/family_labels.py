"""Y1's first gate, `docs/design/family-and-clinical-gates.md`: transmission labels against v5.0q.

Phases HG002's heterozygous SNVs by transmission from HG003 and HG004
(`phasebudget.truth.phase_by_transmission`; child v5.0q, parents v4.2.1 with
their BEDs) and compares the relation of consecutive labelled sites with
v5.0q's phase, overall and per stratum.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/family_labels.py`
"""

from __future__ import annotations

import gzip
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import genome_common as common  # noqa: E402
from adjudicate_truth import clopper_pearson  # noqa: E402

from phasebudget.truth import Regions, phase_by_transmission, read_calls  # noqa: E402

WG = common.ROOT / "data" / "giab_wg"
V5 = common.ROOT / "data" / "giab_v5"
CHILD = V5 / "HG002_GRCh38_v5.0q_smvar.vcf.gz"
CHILD_BED = V5 / "HG002_GRCh38_v5.0q_smvar.benchmark.bed"
FATHER = WG / "HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
MOTHER = WG / "HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
FATHER_BED = WG / "HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
MOTHER_BED = WG / "HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
RULE = 1e-4


CALLS = common.ROOT / "data" / "trio_calls"
CALLED = (CALLS / "hg002.pass.vcf.gz", CALLS / "hg003.pass.vcf.gz", CALLS / "hg004.pass.vcf.gz")


def low_quality(path: Path, chrom: str, min_gq: int) -> set:
    """1-based positions on `chrom` whose call has GQ below `min_gq`."""
    out = set()
    if min_gq <= 0:
        return out
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(chrom + "\t"):
                f = line.split("\t", 10)
                gq = dict(zip(f[8].split(":"), f[9].rstrip("\n").split(":"))).get("GQ", ".")
                if gq == "." or int(float(gq)) < min_gq:
                    out.add(int(f[1]))
    return out


def transmission(chrom: str, called: bool = False, min_gq: int = 0):
    """`(labels, reasons, child calls)`: HG002's transmission labels on one chromosome.

    From the curated genotypes of step 1, or with `called` from the Clair3 PASS
    calls, where a site missing from a parent's VCF counts as homozygous
    reference. With `min_gq` (step 3), a child call below it is dropped and a
    parent call below it makes the site unusable for that parent.
    """
    if called:
        low = [low_quality(path, chrom, min_gq) for path in CALLED]
        child = {p: c for p, c in read_calls(CALLED[0], chrom).items() if p not in low[0]}

        def usable(holes):  # everywhere except the low-quality positions, as 0-based intervals
            edges = [0] + [x for p in sorted(holes) for x in (p - 1, p)] + [10 ** 10]
            return Regions([(a, b) for a, b in zip(edges[::2], edges[1::2]) if a < b])

        labels, reasons = phase_by_transmission(
            child, read_calls(CALLED[1], chrom), read_calls(CALLED[2], chrom),
            usable(low[1]), usable(low[2]))
        return labels, reasons, child
    child_bed = Regions.from_bed(CHILD_BED, chrom)
    child = {p: c for p, c in read_calls(CHILD, chrom).items() if child_bed.contains(p)}
    labels, reasons = phase_by_transmission(
        child, read_calls(FATHER, chrom), read_calls(MOTHER, chrom),
        Regions.from_bed(FATHER_BED, chrom), Regions.from_bed(MOTHER_BED, chrom))
    return labels, reasons, child


def chromosome(chrom: str):
    """`(reasons, [(pos1, pos2, label same, v5.0q same)])` for one chromosome."""
    labels, reasons, child = transmission(chrom)
    sites = sorted(p for p in labels if child[p].phased)
    pairs = [(a, b, labels[a][0] == labels[b][0], child[a].genotype[0] == child[b].genotype[0])
             for a, b in zip(sites, sites[1:])]
    return reasons, pairs


def main() -> None:
    from stratified_gate import strata

    chroms = list(common.CONFIRMATORY)
    print(f"Y1 gate, docs/design/family-and-clinical-gates.md: child {CHILD.name} inside "
          f"{CHILD_BED.name}; father {FATHER.name} in {FATHER_BED.name}; mother {MOTHER.name} in "
          f"{MOTHER_BED.name}; biallelic SNVs; relations of consecutive labelled sites against "
          f"v5.0q phase; chromosomes {' '.join(chroms)}; rule: upper 95% bound <= {RULE:g} "
          "overall and in the easy stratum")
    with ProcessPoolExecutor(6) as pool:
        results = dict(zip(chroms, pool.map(chromosome, chroms)))
    stratum = strata(chroms)
    reasons: Counter = Counter()
    by: dict = {}
    for chrom, (r, pairs) in results.items():
        reasons += r
        for a, b, label_same, truth_same in pairs:
            for key in ("all", stratum(chrom, a, b)):
                n, k = by.get(key, (0, 0))
                by[key] = (n + 1, k + (label_same != truth_same))
    candidates = reasons["candidate"]
    print(f"\nchild heterozygous biallelic SNVs in the v5.0q BED: {candidates}")
    for reason, count in sorted(reasons.items(), key=lambda item: -item[1]):
        if reason != "candidate":
            print(f"  {reason:<40} {count:>9}  {count / candidates:6.1%}")
    print("\nstratum      relations  disagree  rate      95% CI")
    verdict = True
    for key in ("all", "easy", "difficult", "tandem", "lowmap", "segdup"):
        n, k = by.get(key, (0, 0))
        lo, hi = clopper_pearson(k, n) if n else (0.0, 1.0)
        print(f"  {key:<10} {n:>10}  {k:>8}  {k / max(1, n):<8.2e}  {lo:.2e} to {hi:.2e}")
        if key in ("all", "easy") and hi > RULE:
            verdict = False
    print(f"\nY1 gate: {'continue' if verdict else 'stop'}")


if __name__ == "__main__":
    main()
