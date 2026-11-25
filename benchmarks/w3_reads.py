"""Round 4, W3b: the reads around every scored pair, for the neural score.

For a platform, depth and phaser, each consecutive phased pair that the truth
phases (the pairs of `w3_features.py`, in the same order) gets the window of
`docs/design/w3b-neural-score.md`: 8 phased variants of its block around it,
and up to 32 of the reads spanning either variant of the pair. Per read and
window slot: +1 when the read calls the phaser's first-haplotype allele, -1
for the other allele, 0 for no call, and the call's weight. Built from the
stored `realign` fragments and the phaser VCF only; `wrong` is not stored here.

Output: `data/w3_reads/<platform>_<depth>_<phaser>.npz` with `chrom`, `pos1`,
`pos2`, `calls` (int8, pairs x 32 x 8), `weights` (uint8, same shape),
`read_calls` (uint16, pairs x 32: variants each read calls) and `reads`
(uint16: reads spanning the pair before sampling).

Usage: `PHASEBUDGET_PLATFORM=cal_v5 PYTHONPATH=src:benchmarks .venv/bin/python
benchmarks/w3_reads.py --depth 30 --phaser whatshap`
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

BEFORE, AFTER, MAX_READS = 3, 3, 32
SLOTS = BEFORE + 2 + AFTER
OUT = HERE.parent / "data" / "w3_reads"


def window(phased: list, i: int, blocks) -> list:
    """Variant indices of the window slots around the pair `(phased[i], phased[i + 1])`; -1
    where a slot falls past either end of `phased` or in another block."""
    block = blocks[phased[i]]
    out = []
    for k in range(i - BEFORE, i + 2 + AFTER):
        v = phased[k] if 0 <= k < len(phased) else -1
        out.append(v if v >= 0 and blocks[v] == block else -1)
    return out


def chromosome(chrom: str, depth: str, phaser: str) -> tuple:
    import numpy as np

    import genome_common as common
    from phase_confidence import phasing_from_vcf
    from phasebudget.confidence import phased_pairs

    extractor = common._module("extract_fragments.py")
    sc = common._module("score_hg002.py")
    variants = sc.variant_table(common.hetsnv_path(chrom))
    matrix, positions, _, _ = extractor.load_fragments(
        str(common.fragments_path(chrom, depth, "realign")))
    truth = common.truth_at(chrom, variants, positions)
    vcf = common.chrom_dir(chrom) / "external" / f"{phaser}_{depth}.vcf"
    vcf = vcf if vcf.exists() else vcf.with_name(vcf.name + ".gz")
    phasing = phasing_from_vcf(str(vcf), common.VCF_SAMPLE, chrom, positions)
    phased = [j for j, a in enumerate(phasing.alleles) if a is not None]
    index = {v: k for k, v in enumerate(phased)}
    reads, called, cover = [], [], [[] for _ in positions]
    for r, f in enumerate(matrix):
        calls = {f.start + k: (a, w) for k, (a, w) in enumerate(zip(f.alleles, f.weights))
                 if a is not None}
        reads.append(calls)
        called.append(len(calls))
        for v in range(f.start, f.end):
            cover[v].append(r)
    pairs = [(a, b) for a, b in phased_pairs(phasing)
             if truth[a] is not None and truth[b] is not None]
    n = len(pairs)
    c_arr = np.zeros((n, MAX_READS, SLOTS), dtype=np.int8)
    w_arr = np.zeros((n, MAX_READS, SLOTS), dtype=np.uint8)
    rc_arr = np.zeros((n, MAX_READS), dtype=np.uint16)
    nr_arr = np.zeros(n, dtype=np.uint16)
    for row, (a, b) in enumerate(pairs):
        span = sorted(set(cover[a]) | set(cover[b]))
        nr_arr[row] = min(len(span), 65535)
        if len(span) > MAX_READS:
            span = sorted(random.Random(positions[a]).sample(span, MAX_READS))
        slots = window(phased, index[a], phasing.blocks)
        for k, r in enumerate(span):
            calls = reads[r]
            rc_arr[row, k] = min(called[r], 65535)
            for s, v in enumerate(slots):
                c = calls.get(v) if v >= 0 else None
                if c is not None:
                    c_arr[row, k, s] = 1 if c[0] == phasing.alleles[v] else -1
                    w_arr[row, k, s] = min(255, max(0, int(round(c[1]))))
    pos1 = np.array([positions[a] + 1 for a, _ in pairs])
    pos2 = np.array([positions[b] + 1 for _, b in pairs])
    return np.array([chrom] * n), pos1, pos2, c_arr, w_arr, rc_arr, nr_arr


def main() -> None:
    import numpy as np

    import genome_common as common

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", required=True)
    parser.add_argument("--phaser", required=True)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    parts = [chromosome(c, args.depth, args.phaser) for c in common.CONFIRMATORY]
    keys = ("chrom", "pos1", "pos2", "calls", "weights", "read_calls", "reads")
    path = OUT / f"{common.PLATFORM}_{args.depth}_{args.phaser}.npz"
    np.savez_compressed(path, **{k: np.concatenate([p[i] for p in parts])
                                 for i, k in enumerate(keys)})
    print(f"{path.name}: {sum(len(p[0]) for p in parts)} pairs")


if __name__ == "__main__":
    main()
