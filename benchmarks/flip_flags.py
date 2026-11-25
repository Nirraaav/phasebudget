"""Z7's gate, the October directions memo (development repository): do low-support flips flag genotype errors?

A Platinum 2017 heterozygous SNV of HG001 inside the Platinum Pedigree
high-quality regions is flagged when both consecutive pairs it belongs to
have a Wilson lower bound below 0.5; it is not heterozygous when the Pedigree
truth has no heterozygous call (phased or not) at it. Data seen before.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/flip_flags.py`
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import adjudicate_truth as adj  # noqa: E402
import confidence_score as cs  # noqa: E402
import genome_common as common  # noqa: E402

LOW = 0.5
ENRICHMENT = 20
RECALL = 0.3


def pedigree_hets(chroms) -> set:
    """`(chrom, pos1)` of NA12878 heterozygous biallelic SNVs in the Pedigree truth, any phase."""
    import pysam

    out = set()
    with pysam.TabixFile(str(adj.TRUTHSET)) as tabix:
        column = 9 + tabix.header[-1].split("\t")[9:].index(adj.SAMPLE)
        for chrom in chroms:
            for line in tabix.fetch(chrom):
                f = line.split("\t")
                gt = f[column].split(":")[0].replace("|", "/")
                if gt in ("0/1", "1/0") and len(f[3]) == 1 and len(f[4]) == 1:
                    out.add((chrom, int(f[1])))
    return out


def main() -> None:
    chroms = list(common.CONFIRMATORY)
    print(f"Z7 gate (development; data seen): HG001 ONT sup rows, 30x; site flagged when both "
          f"adjacent pairs have Wilson lower bound < {LOW}; not heterozygous = no heterozygous "
          f"Pedigree truth call; sites inside {adj.HQ.name}; continue if enrichment >= "
          f"{ENRICHMENT} and recall >= {RECALL} for both phasers")
    hets = pedigree_hets(chroms)
    hq = adj.load_bed(adj.HQ, chroms)
    verdict = True
    for phaser in cs.PHASERS:
        left, right = {}, {}
        for c, p1, key, wrong, p2, _ in cs.load_rows("hg001_ont", "30", phaser, chroms):
            right[(c, p1)] = key[0]
            left[(c, p2)] = key[0]
        sites = [s for s in left.keys() & right.keys() if adj.inside(hq, *s)]
        table = {(f, n): 0 for f in (True, False) for n in (True, False)}
        for s in sites:
            table[(max(left[s], right[s]) < LOW, s not in hets)] += 1
        flagged = table[(True, True)] + table[(True, False)]
        unflagged = table[(False, True)] + table[(False, False)]
        rate_f = table[(True, True)] / max(1, flagged)
        rate_u = table[(False, True)] / max(1, unflagged)
        bad = table[(True, True)] + table[(False, True)]
        recall = table[(True, True)] / max(1, bad)
        enrichment = rate_f / rate_u if rate_u else float("inf")
        print(f"\n{phaser}: sites with two scored pairs {len(sites)}, not heterozygous {bad}")
        print(f"  flagged {flagged}: not heterozygous {table[(True, True)]} ({rate_f:.3f})")
        print(f"  unflagged {unflagged}: not heterozygous {table[(False, True)]} ({rate_u:.2e})")
        print(f"  enrichment {enrichment:.1f}, recall {recall:.3f}")
        verdict &= enrichment >= ENRICHMENT and recall >= RECALL
    print(f"\nZ7 gate: {'continue' if verdict else 'stop'}")


if __name__ == "__main__":
    main()
