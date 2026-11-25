"""Keep a random fraction of the reads of a BAM, one draw per read name.

The real-data plan registers a 50 per cent subsample drawn with seed 0 before
extraction. Each read name gets one draw from `random.Random(seed)` in the
order the name first appears in the file, so a read's primary, secondary and
supplementary records are kept or dropped together, and the same file, seed
and fraction always give the same subsample.
"""

from __future__ import annotations

import argparse
import random


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bam")
    parser.add_argument("out")
    parser.add_argument("--fraction", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    import pysam

    rng = random.Random(args.seed)
    keep: dict[str, bool] = {}
    kept = total = 0
    with pysam.AlignmentFile(args.bam) as source, \
            pysam.AlignmentFile(args.out, "wb", template=source) as target:
        for record in source:
            total += 1
            name = record.query_name
            if name not in keep:
                keep[name] = rng.random() < args.fraction
            if keep[name]:
                target.write(record)
                kept += 1
    pysam.index(args.out)
    names = sum(keep.values())
    print(f"seed {args.seed}, fraction {args.fraction:g}: kept {kept} of {total} records, "
          f"{names} of {len(keep)} read names")


if __name__ == "__main__":
    main()
