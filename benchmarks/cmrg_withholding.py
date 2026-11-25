"""Y2 step 2, `docs/design/family-and-clinical-gates.md`: withholding against GIAB's CMRG genes.

Reads the per-gene counts written by `chain_relations.py --clinical
--stratified --per-gene` and tests, for each phaser on HG002 ONT hac at beta
1e-3, whether CMRG genes have lower trusted shares (one-sided Mann-Whitney U)
and are over-represented below a trusted share of 0.5 (one-sided Fisher).

With `--duplicated` (step 3) the list is instead the protein-coding genes
with an unprocessed pseudogene named after them in GENCODE v50 basic.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/cmrg_withholding.py PER_GENE.json`
"""

from __future__ import annotations

import json
import statistics
import sys
from pathlib import Path

from scipy.stats import fisher_exact, mannwhitneyu

CMRG = Path(__file__).resolve().parent.parent / "data" / "cmrg" / \
    "GRCh38_CMRG_benchmark_gene_coordinates.bed"
MIN_RELATIONS = 100
LOW = 0.5


GENCODE = CMRG.parent.parent / "gencode" / "gencode.v50.basic.annotation.gtf.gz"


def duplicated() -> set:
    """Protein-coding genes X with an unprocessed pseudogene named X + 'P' + digits."""
    import gzip
    import re

    coding, pseudo = set(), set()
    with gzip.open(GENCODE, "rt") as handle:
        for line in handle:
            f = line.split("\t", 9)
            if line.startswith("#") or f[2] != "gene":
                continue
            name = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
            kind = re.search(r'gene_type "([^"]+)"', f[8]).group(1)
            if kind == "protein_coding":
                coding.add(name)
            elif kind in ("unprocessed_pseudogene", "transcribed_unprocessed_pseudogene"):
                pseudo.add(name)
    return {m.group(1) for p in pseudo if (m := re.fullmatch(r"(.+)P\d+", p))} & coding


def main() -> None:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    per_gene = json.loads(Path(args[0]).read_text())
    if "--duplicated" in sys.argv:
        cmrg, label = duplicated(), "duplicated (GENCODE v50 basic unprocessed pseudogene X+P+digits)"
    else:
        cmrg = {line.split()[3] for line in CMRG.read_text().splitlines() if line.strip()}
        label = f"CMRG ({CMRG.name})"
    print(f"Y2: {len(cmrg)} {label} genes; clinical genes with at least "
          f"{MIN_RELATIONS} relations; beta 1e-3, 30x, per-stratum thresholds; one-sided tests at "
          f"0.05; low = trusted share below {LOW}")
    for key, genes in per_gene.items():
        shares = {g: tr / n for g, (n, tr) in genes.items() if n >= MIN_RELATIONS}
        inside = [v for g, v in shares.items() if g in cmrg]
        outside = [v for g, v in shares.items() if g not in cmrg]
        u = mannwhitneyu(inside, outside, alternative="less")
        table = [[sum(v < LOW for v in inside), sum(v >= LOW for v in inside)],
                 [sum(v < LOW for v in outside), sum(v >= LOW for v in outside)]]
        f = fisher_exact(table, alternative="greater")
        print(f"\n{key}: {len(inside)} listed and {len(outside)} other genes")
        print(f"  median trusted share: listed {statistics.median(inside):.3f}, other "
              f"{statistics.median(outside):.3f}; Mann-Whitney p {u.pvalue:.2e}")
        print(f"  share < {LOW}: listed {table[0][0]} of {len(inside)}, other {table[1][0]} of "
              f"{len(outside)}; odds ratio {f.statistic:.2f}, Fisher p {f.pvalue:.2e}")
        low = sorted((v, g) for g, v in shares.items() if v < LOW)
        print(f"  low genes in the list: {', '.join(g for v, g in low if g in cmrg)}")


if __name__ == "__main__":
    main()
