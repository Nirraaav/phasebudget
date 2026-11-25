"""Y1 step 3, `docs/design/family-and-clinical-gates.md`: choose the GQ filter on odd chromosomes.

For `g` in {0, 10, 20}, counts on the odd confirmatory chromosomes the sites
labelled from Clair3 calls and the disagreements with v5.0q among consecutive
labelled sites, then applies the stated choice rule.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/family_gq.py`
"""

from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import genome_common as common  # noqa: E402
from family_labels import read_calls, transmission  # noqa: E402

GQS = (0, 10, 20)
KEEP = 0.8


def one(job):
    chrom, gq = job
    labels, _, called_child = transmission(chrom, called=True, min_gq=gq)
    v5 = {p: c for p, c in read_calls(common.ROOT / "data" / "giab_v5" /
                                      "HG002_GRCh38_v5.0q_smvar.vcf.gz", chrom).items()
          if c.phased and tuple(c.genotype) in ((0, 1), (1, 0))}
    sites = sorted(p for p in labels if p in v5 and v5[p].alt == called_child[p].alt)
    disagree = sum((labels[a][0] == labels[b][0]) != (v5[a].genotype[0] == v5[b].genotype[0])
                   for a, b in zip(sites, sites[1:]))
    return gq, len(labels), len(sites), disagree


def main() -> None:
    odd = [c for c in common.CONFIRMATORY if int(c[3:]) % 2]
    print(f"Y1 step 3 choice of GQ on {' '.join(odd)}: labelled sites (all, and phased "
          f"heterozygous in v5.0q), disagreements with v5.0q among consecutive labelled sites; "
          f"keep >= {KEEP} of g = 0's labelled sites")
    totals = {g: [0, 0, 0] for g in GQS}
    with ProcessPoolExecutor(6) as pool:
        for gq, n, m, k in pool.map(one, product(odd, GQS)):
            for i, v in enumerate((n, m, k)):
                totals[gq][i] += v
    for g, (n, m, k) in totals.items():
        print(f"  g {g:>2}: labelled {n}, in v5.0q {m}, disagree {k} ({k / max(1, m - 1):.2e})")
    allowed = [g for g in GQS if totals[g][0] >= KEEP * totals[0][0]]
    chosen = min(allowed, key=lambda g: (totals[g][2], g))
    print(f"chosen g: {chosen}")


if __name__ == "__main__":
    main()
