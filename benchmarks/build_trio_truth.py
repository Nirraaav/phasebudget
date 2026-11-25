"""Build the transmission truth for HG002 on one chromosome and score GIAB's phased files.

Inputs are the GIAB v4.2.1 benchmark VCF and BED of HG002, HG003 and HG004,
and optionally the two phased HG002 files from `SupplementaryFiles/`. The
rules are `docs/design/truth-sets.md`. The script prints where every candidate
site went, then each phased file's per-site discordance with the truth, and
writes the truth as a phased VCF when `--out` is given.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from phasebudget.truth import Regions, compare_phasing, phase_by_transmission, read_calls, write_vcf


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chromosome", default="chr20")
    parser.add_argument("--data", type=Path, default=Path("data/giab"))
    parser.add_argument("--child", default="HG002_chr20_v4.2.1.vcf.gz")
    parser.add_argument("--father", default="HG003_chr20_v4.2.1.vcf.gz")
    parser.add_argument("--mother", default="HG004_chr20_v4.2.1.vcf.gz")
    parser.add_argument("--father-bed", default="HG003_chr20_v4.2.1.bed")
    parser.add_argument("--mother-bed", default="HG004_chr20_v4.2.1.bed")
    parser.add_argument(
        "--phased",
        nargs="*",
        default=[
            "HG002_chr20_StrandSeqANDTrio.vcf.gz:1",
            "HG002_chr20_hifiasm_phasetransfer.vcf.gz:PATMAT",
        ],
        help="FILE:SET[,SET...], the phase sets that are paternal|maternal",
    )
    parser.add_argument("--out", type=Path)
    parser.add_argument("--sample", default="HG002", help="sample name written to --out")
    args = parser.parse_args()

    chromosome = args.chromosome
    data = args.data
    print(f"chromosome {chromosome}; data directory {data}")
    print(f"child {args.child}; father {args.father}; mother {args.mother}")
    print(f"father BED {args.father_bed}; mother BED {args.mother_bed}")
    print("rules: docs/design/truth-sets.md; GIAB v4.2.1 benchmark calls")

    child = read_calls(data / args.child, chromosome)
    father = read_calls(data / args.father, chromosome)
    mother = read_calls(data / args.mother, chromosome)
    phased, reasons = phase_by_transmission(
        child,
        father,
        mother,
        Regions.from_bed(data / args.father_bed, chromosome),
        Regions.from_bed(data / args.mother_bed, chromosome),
    )
    candidates = reasons["candidate"]
    print()
    print(f"child heterozygous biallelic SNVs: {candidates}")
    for reason, count in sorted(reasons.items(), key=lambda item: -item[1]):
        if reason == "candidate":
            continue
        print(f"  {reason:<40} {count:>7}  {count / candidates:6.1%}")

    if args.phased:
        print()
        header = (
            f"{'phased file':<45} {'compared':>9} {'discordant':>11} {'rate':>8} "
            f"{'not oriented':>13} {'unphased/absent':>16}"
        )
        print(header)
        print("-" * len(header))
    for spec in args.phased:
        name, _, sets = spec.partition(":")
        other = read_calls(data / name, chromosome)
        result = compare_phasing(phased, other, sets.split(","))
        print(
            f"{name:<45} {result.compared:>9} {result.discordant:>11} "
            f"{result.discordance:>8.4%} {result.not_oriented:>13} {result.unphased_or_absent:>16}"
        )

    if args.out is not None:
        write_vcf(args.out, chromosome, child, phased, args.sample)
        print()
        print(f"wrote {len(phased)} phased sites to {args.out}")


if __name__ == "__main__":
    main()
