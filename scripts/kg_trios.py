"""Subset the 1000 Genomes high-coverage phased panel to the trios W1's consensus truth uses.

W1 amendment 1 (`docs/design/w1-ancestry-plan.md`): every test genome in
`docs/design/w1-candidates.tsv` whose father and mother are both in the
3,202-sample pedigree, and NA12878 with her parents for the gate. Writes
`data/kg_panel/trios/<chrom>.vcf.gz` and its index for the 21 confirmatory
autosomes, keeping sites where a kept sample carries an alt allele.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python scripts/kg_trios.py`
"""

from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))

import genome_common as common  # noqa: E402

PANEL = "1kGP_high_coverage_Illumina.{c}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"


def children() -> list[str]:
    with open(common.ROOT / "docs" / "design" / "w1-candidates.tsv") as handle:
        handle.readline()
        samples = [line.split("\t")[0] for line in handle]
    return [s for s in samples + ["NA12878"] if common.kg_parents(s)]


def subset(chrom: str, samples: list[str]) -> str:
    import pysam
    import pysam.bcftools  # a submodule: `import pysam` alone does not load it

    out = common.KG_PANEL / "trios" / f"{chrom}.vcf.gz"
    if not Path(str(out) + ".tbi").exists():
        part = out.with_name(out.name + ".part")
        pysam.bcftools.view("-s", ",".join(samples), "-c", "1", "-Oz", "-o", str(part),
                            str(common.KG_PANEL / PANEL.format(c=chrom)), catch_stdout=False)
        part.rename(out)
        pysam.tabix_index(str(out), preset="vcf", force=True)
    return chrom


def main() -> None:
    kids = children()
    samples = sorted({s for c in kids for s in (c, *common.kg_parents(c))})
    print(f"{len(kids)} children, {len(samples)} samples: {' '.join(kids)}", flush=True)
    (common.KG_PANEL / "trios").mkdir(exist_ok=True)
    with ProcessPoolExecutor(11) as pool:
        for chrom in pool.map(subset, common.CONFIRMATORY, [samples] * len(common.CONFIRMATORY)):
            print(f"{chrom} done", flush=True)


if __name__ == "__main__":
    main()
