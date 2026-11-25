"""Follow-up gate G-A (`docs/followup/error-controlled-phasing.md`): does read support flag errors?

For each pair of consecutive phased variants that a tool puts in one block,
the support margin is the number of reads calling both variants whose
relation agrees with the called one, minus the number that disagree; reads
come from the run's `realign` fragments (all reads, not a selection). Pairs
are ranked by margin, then by the number of supporting reads, lowest first.
The script prints the AUC of the rank for a wrong pair (against the truth)
and the risk-coverage curve: the errors left when only the top share of
pairs is trusted.

Development only: chromosome 20, which enters no PhaseBudget test. Usage:
`PHASEBUDGET_PLATFORM=ont PYTHONPATH=src:benchmarks .venv/bin/python
benchmarks/gate_evidence_score.py chr20 --depths full 10`
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402

extractor = common._module("extract_fragments.py")
score = common._module("score_hg002.py")
COVERAGES = (1.0, 0.999, 0.995, 0.99, 0.98, 0.95)
# Higher is more trusted; ties are broken by the number of supporting reads.
SCORES = {
    "margin": lambda r: (r[0], r[1]),
    # The disagreeing share of the reads calling both variants (added after gate G-C).
    "conflict": lambda r: (-(r[1] - r[0]) / 2 / r[1] if r[1] else -1.0, r[1]),
    # Wilson lower bound (z = 2) on the agreeing share: low support and heavy conflict
    # both lower it (added after the first two scores were compared on ONT chromosome 20).
    "wilson": lambda r: (wilson_lower((r[1] + r[0]) // 2, r[1]), r[1]),
}


def wilson_lower(k: int, n: int, z: float = 2.0) -> float:
    """Lower end of the Wilson interval for `k` successes in `n`; -1 when `n` is 0."""
    if n == 0:
        return -1.0
    p = k / n
    centre = p + z * z / (2 * n)
    spread = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (centre - spread) / (1 + z * z / n)


def pair_scores(fragments, alleles, blocks, truth):
    """`[(margin, support, wrong)]` per consecutive phased pair in one block with truth."""
    phased = [j for j, a in enumerate(alleles) if a is not None and truth[j] is not None]
    nxt = {a: b for a, b in zip(phased, phased[1:]) if blocks[a] == blocks[b]}
    agree, disagree = {}, {}
    for f in fragments:
        calls = {f.start + k: a for k, a in enumerate(f.alleles) if a is not None}
        for j, a in calls.items():
            k = nxt.get(j)
            if k is None or k not in calls:
                continue
            same_read = (a == calls[k])
            same_call = (alleles[j] == alleles[k])
            target = agree if same_read == same_call else disagree
            target[j] = target.get(j, 0) + 1
    out = []
    for j, k in nxt.items():
        wrong = (alleles[j] == alleles[k]) != (truth[j] == truth[k])
        a, d = agree.get(j, 0), disagree.get(j, 0)
        out.append((a - d, a + d, wrong))
    return out


def auc(rows):
    """P(a wrong pair ranks below a right one), ties counted half."""
    ranked = sorted(rows, key=lambda r: (r[0], r[1]))
    wrong_total = sum(r[2] for r in ranked)
    right_total = len(ranked) - wrong_total
    if not wrong_total or not right_total:
        return float("nan")
    seen_right = 0
    total = 0.0
    i = 0
    while i < len(ranked):
        j = i
        while j < len(ranked) and ranked[j][:2] == ranked[i][:2]:
            j += 1
        w = sum(r[2] for r in ranked[i:j])
        r = (j - i) - w
        total += w * (right_total - seen_right - r / 2)
        seen_right += r
        i = j
    return total / (wrong_total * right_total)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chrom")
    parser.add_argument("--depths", nargs="+", default=["full", "10"])
    parser.add_argument("--tools", nargs="+", default=["whatshap", "longphase"])
    args = parser.parse_args()
    variants = score.variant_table(common.hetsnv_path(args.chrom))
    truth = (common.v5_truth(args.chrom, variants) if common.SAMPLE == "HG002"
             else common.load_truths(args.chrom, [p - 1 for p, _, _ in variants])[1])
    print(f"gate G-A, platform {common.PLATFORM}, {args.chrom}; support margin from all realign "
          "reads; errors left among the trusted top share of pairs")
    for depth in args.depths:
        matrix, _, _, _ = extractor.load_fragments(
            str(common.fragments_path(args.chrom, depth, "realign")))
        for tool in args.tools:
            d = common.chrom_dir(args.chrom) / "external"
            vcf = d / f"{tool}_{depth}.vcf"
            vcf = vcf if vcf.exists() else d / f"{tool}_{depth}.vcf.gz"
            if not vcf.exists():
                print(f"{depth:>4} {tool}: no output")
                continue
            p = score.whatshap_phasing(vcf, variants)
            rows = pair_scores(matrix.fragments, p.alleles, p.blocks, truth)
            n, errors = len(rows), sum(r[2] for r in rows)
            low = sum(1 for r in rows if r[0] <= 0)
            print(f"{depth:>4} {tool:>9}: pairs {n}, wrong {errors}, pairs with margin <= 0: {low}")
            for name, key in SCORES.items():
                ranked = sorted(rows, key=key, reverse=True)
                curve = ", ".join(f"{c:.3f}: {sum(r[2] for r in ranked[: int(round(c * n))])}"
                                  for c in COVERAGES)
                as_rank = [(key(r)[0], key(r)[1], r[2]) for r in rows]
                print(f"{'':>15}{name:>8}: AUC {auc(as_rank):.3f}; wrong left at coverage {curve}")


if __name__ == "__main__":
    main()
