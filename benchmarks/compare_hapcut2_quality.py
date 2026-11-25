"""Development: HapCUT2's own per-variant quality against the frozen read-support score.

Both scores rank the same pairs, HapCUT2's consecutive phased variants in
one block, from the ONT run's `hapcut2_{depth}.blocks` (HG002 ONT sup, 21
chromosomes, against v5.0q). The scores are:
- `hapcut2`: the smaller of the two variants' mismatch qualities (column 11
  of the blocks file, the VCF's PQ: a phred-scaled probability that the
  variant is phased wrongly relative to its block);
- `support`: `phasebudget.confidence.score` from the run's `realign`
  fragments.

It prints each score's AUC for wrong pairs and the errors left at matched
shares of pairs trusted. The data were seen in PhaseBudget, so this is
development only.

Usage: `PHASEBUDGET_PLATFORM=ont PYTHONPATH=src:benchmarks .venv/bin/python
benchmarks/compare_hapcut2_quality.py --depths 30 20 10`
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402
from gate_evidence_score import auc  # noqa: E402

from phasebudget._types import Phasing  # noqa: E402
from phasebudget.confidence import pair_support, score  # noqa: E402

extractor = common._module("extract_fragments.py")
sc = common._module("score_hg002.py")
COVERAGES = (0.95, 0.98, 0.99, 0.995, 0.999, 1.0)


def hapcut2_blocks(path: Path, variants):
    """`(Phasing, quality per variant)` from a HapCUT2 blocks file; pruned variants unphased."""
    index = {(p, r, a): j for j, (p, r, a) in enumerate(variants)}
    alleles = [None] * len(variants)
    blocks = [None] * len(variants)
    quality = [None] * len(variants)
    block = -1
    for line in path.read_text().splitlines():
        if line.startswith("BLOCK"):
            block += 1
            continue
        if line.startswith("*"):
            continue
        f = line.split("\t")
        j = index.get((int(f[4]), f[5], f[6]))
        if j is None or f[1] not in ("0", "1") or f[1] == f[2] or f[8] == "1":
            continue
        alleles[j] = int(f[1])
        blocks[j] = block
        quality[j] = float(f[10]) if f[10] not in (".", "") else 0.0
    return Phasing(alleles=tuple(alleles), blocks=tuple(blocks)), quality


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depths", nargs="+", default=["30", "20", "10"])
    args = parser.parse_args()
    print(f"development: HapCUT2 quality against the read-support score, platform "
          f"{common.PLATFORM}, 21 chromosomes, HapCUT2 phasings")
    for depth in args.depths:
        rows = []  # (support score, hapcut2 score, wrong)
        for chrom in common.CONFIRMATORY:
            variants = sc.variant_table(common.hetsnv_path(chrom))
            truth = common.primary_truth(chrom, variants)
            phasing, quality = hapcut2_blocks(
                common.chrom_dir(chrom) / "external" / f"hapcut2_{depth}.blocks", variants)
            matrix, _, _, _ = extractor.load_fragments(
                str(common.fragments_path(chrom, depth, "realign")))
            for s in pair_support(matrix, phasing):
                a, b = s.first, s.second
                if truth[a] is None or truth[b] is None:
                    continue
                wrong = (phasing.alleles[a] == phasing.alleles[b]) != (truth[a] == truth[b])
                rows.append((score(s), (min(quality[a], quality[b]), 0), wrong))
        n, errors = len(rows), sum(r[2] for r in rows)
        print(f"  {depth}x: {n} pairs, {errors} wrong")
        for name, k in (("support", 0), ("hapcut2", 1)):
            ranked = sorted(rows, key=lambda r: r[k], reverse=True)
            curve = ", ".join(f"{c}: {sum(r[2] for r in ranked[: int(round(c * n))])}"
                              for c in COVERAGES)
            value = auc([(r[k][0], r[k][1], r[2]) for r in rows])
            print(f"    {name:>8}: AUC {value:.3f}; wrong left at share trusted {curve}")


if __name__ == "__main__":
    main()
