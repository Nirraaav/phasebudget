"""Round 4, W3: read-evidence features of every scored pair, for a learned score.

For a platform, depth and phaser, each consecutive phased pair that the
truth phases (the population `confidence_score.rows_for` scores) gets the
features below, computed from the stored `realign` fragments, the phaser VCF
and the phaser input only. `wrong` is the truth's verdict, used as the
training label and by the scorer, never as a feature.

Features: reads covering both variants with a called allele at each, the
agreeing and disagreeing ones, the frozen score's Wilson bound, the same
counts weighted by the smaller of the two allele qualities, reads spanning the
pair with a gap at either variant, per variant the reads covering it, the
share of them calling an allele, the allele balance and mean allele quality,
the distance between the variants, indel flags and lengths, and the number of
phaser-input variants within 10 kb.

Output: `data/w3_features/<platform>_<depth>_<phaser>.npz` with `chrom`,
`pos1`, `pos2`, `wrong` and a float matrix `x` whose columns are `NAMES`.

Usage: `PHASEBUDGET_PLATFORM=cal_r10 PYTHONPATH=src:benchmarks .venv/bin/python
benchmarks/w3_features.py --depth 30 --phaser whatshap`
"""

from __future__ import annotations

import argparse
import bisect
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

NAMES = ("reads", "agree", "disagree", "wilson", "w_agree", "w_disagree", "spanning", "gaps",
         "cov_a", "cov_b", "call_a", "call_b", "balance_a", "balance_b", "qual_a", "qual_b",
         "distance", "indel_a", "indel_b", "len_a", "len_b", "density")
OUT = HERE.parent / "data" / "w3_features"


def chromosome(chrom: str, depth: str, phaser: str) -> list:
    import genome_common as common
    from phase_confidence import phasing_from_vcf
    from phasebudget.confidence import phased_pairs, wilson_lower

    extractor = common._module("extract_fragments.py")
    sc = common._module("score_hg002.py")
    variants = sc.variant_table(common.hetsnv_path(chrom))
    matrix, positions, _, _ = extractor.load_fragments(
        str(common.fragments_path(chrom, depth, "realign")))
    truth = common.truth_at(chrom, variants, positions)
    vcf = common.chrom_dir(chrom) / "external" / f"{phaser}_{depth}.vcf"
    vcf = vcf if vcf.exists() else vcf.with_name(vcf.name + ".gz")
    phasing = phasing_from_vcf(str(vcf), common.VCF_SAMPLE, chrom, positions)
    alleles = {v[0] - 1: (v[1], v[2]) for v in variants}
    n = len(positions)
    pairs = phased_pairs(phasing)
    following = dict(pairs)
    stats = {a: [0, 0, 0.0, 0.0, 0, 0] for a, _ in pairs}  # agree, disagree, w+, w-, spanning, gaps
    cover, calls, ones, quality = [0] * n, [0] * n, [0] * n, [0.0] * n
    for f in matrix:
        called = {}
        for k, (allele, weight) in enumerate(zip(f.alleles, f.weights)):
            v = f.start + k
            cover[v] += 1
            if allele is not None:
                called[v] = (allele, weight)
                calls[v] += 1
                ones[v] += allele
                quality[v] += weight
        for a in range(f.start, f.end - 1):
            b = following.get(a)
            if b is None or b >= f.end:
                continue
            s = stats[a]
            s[4] += 1
            if a not in called or b not in called:
                s[5] += 1
                continue
            (xa, wa), (xb, wb) = called[a], called[b]
            same = (xa == xb) == (phasing.alleles[a] == phasing.alleles[b])
            s[0 if same else 1] += 1
            s[2 if same else 3] += min(wa, wb)
    rows = []
    for a, b in pairs:
        if truth[a] is None or truth[b] is None:
            continue
        agree, disagree, wp, wm, spanning, gaps = stats[a]
        reads = agree + disagree
        ra, rb = alleles[positions[a]], alleles[positions[b]]
        pa, pb = positions[a], positions[b]
        density = bisect.bisect_right(positions, pa + 10_000) - bisect.bisect_left(positions,
                                                                                 pa - 10_000)
        wrong = (phasing.alleles[a] == phasing.alleles[b]) != (truth[a] == truth[b])
        x = (reads, agree, disagree, wilson_lower(agree, reads), wp, wm, spanning, gaps,
             cover[a], cover[b], calls[a] / max(1, cover[a]), calls[b] / max(1, cover[b]),
             abs(ones[a] / max(1, calls[a]) - 0.5), abs(ones[b] / max(1, calls[b]) - 0.5),
             quality[a] / max(1, calls[a]), quality[b] / max(1, calls[b]),
             pb - pa, float(len(ra[0]) != 1 or len(ra[1]) != 1),
             float(len(rb[0]) != 1 or len(rb[1]) != 1), abs(len(ra[0]) - len(ra[1])),
             abs(len(rb[0]) - len(rb[1])), density)
        rows.append((chrom, pa + 1, pb + 1, wrong, x))
    return rows


def main() -> None:
    import numpy as np

    import genome_common as common

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", required=True)
    parser.add_argument("--phaser", required=True)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{common.PLATFORM}_{args.depth}_{args.phaser}.npz"
    rows = []
    for chrom in common.CONFIRMATORY:
        rows += chromosome(chrom, args.depth, args.phaser)
    np.savez_compressed(path, chrom=np.array([r[0] for r in rows]),
                        pos1=np.array([r[1] for r in rows]), pos2=np.array([r[2] for r in rows]),
                        wrong=np.array([r[3] for r in rows]),
                        x=np.array([r[4] for r in rows], dtype=float), names=np.array(NAMES))
    print(f"{path.name}: {len(rows)} pairs, {sum(r[3] for r in rows)} wrong")


if __name__ == "__main__":
    main()
