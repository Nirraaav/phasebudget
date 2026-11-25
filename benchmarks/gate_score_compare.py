"""Follow-up gates G-A and G-B on HG002 HiFi: compare two support scores for `whatshap phase`.

G-C found that most remaining HiFi errors have a positive read margin, but
with many reads disagreeing. This ranks every consecutive phased pair of
stock `whatshap phase` (cap 15, the HiFi tool-level run, 21 chromosomes) by:
- `margin`: agreeing minus disagreeing reads (`gate_evidence_score.py`);
- `conflict`: the disagreeing share of the reads calling both variants.

It prints the AUC for wrong pairs and the errors left when the top share of
pairs is trusted. The data were seen in PhaseBudget, so this is development
only. Usage: `PYTHONPATH=src:benchmarks .venv/bin/python
benchmarks/gate_score_compare.py`
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402
from gate_evidence_score import COVERAGES, SCORES, auc, pair_scores  # noqa: E402
from gate_truth_check import TOOLS, existing  # noqa: E402

extractor = common._module("extract_fragments.py")
score = common._module("score_hg002.py")


def main() -> None:
    rows = []
    for chrom in common.CONFIRMATORY:
        variants = score.variant_table(common.hetsnv_path(chrom))
        v5 = common.v5_truth(chrom, variants)
        wh = score.whatshap_phasing(existing(TOOLS / chrom / "stock_c15.vcf"), variants)
        matrix, _, _, _ = extractor.load_fragments(
            str(common.fragments_path(chrom, "full", "realign")))
        rows.extend(pair_scores(matrix.fragments, wh.alleles, wh.blocks, v5))
    n, errors = len(rows), sum(r[2] for r in rows)
    print(f"gates G-A and G-B, HG002 HiFi full depth, stock whatshap phase cap 15, 21 chromosomes, "
          f"v5.0q: {n} pairs, {errors} wrong")
    for name, key in SCORES.items():
        ranked = sorted(rows, key=key, reverse=True)
        curve = ", ".join(f"{c:.3f}: {sum(r[2] for r in ranked[: int(round(c * n))])}"
                          for c in COVERAGES)
        as_rank = [(key(r)[0], key(r)[1], r[2]) for r in rows]
        print(f"{name:>8}: AUC {auc(as_rank):.3f}; wrong left at coverage {curve}")


if __name__ == "__main__":
    main()
