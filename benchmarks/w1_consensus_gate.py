"""W1 amendment 1's gate: how often 1000 Genomes transmission is wrong, on NA12878.

The consensus truth (`genome_common.consensus_truth`) phases a pair only where
dipcall and transmission from the 1000 Genomes panel give the same paternal
alleles, so a consensus pair is wrong only where the transmission pair is.
The share of 100 kb units holding a wrong transmission pair therefore bounds
the consensus truth's, with no assumption about how dipcall's and the panel's
errors relate.

NA12878 is the one child here with 1000 Genomes parents and an independent
truth. Her variants are the Platinum Pedigree truthset's heterozygous calls
inside its high-quality regions, as the pedigree arm's phaser input selects
them. A pair is two consecutive variants of that list that transmission
labels and the Pedigree truth phases; it is wrong when the two give different
relative phases. Groups as A1 and A2 (`indel` for a pair with an indel, else
the stratum). A (group, beta) is testable when the one-sided 95 per cent exact
upper bound on the share of wrong units is at most beta/3.

`--describe` instead reports, without deciding anything: for HG002 and HG005,
the consensus of dipcall and the GIAB trio transmission truth against v5.0q
and against GIAB's hifiasm phase transfer; and for every test genome whose
dipcall truth exists, dipcall against 1000 Genomes transmission. No phaser
output is read in either mode.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/w1_consensus_gate.py
[--describe]`
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import genome_common as common  # noqa: E402
import w1_gate as gate  # noqa: E402

GROUPS = gate.GROUPS


def count(table: list, a: list, b: list, chrom: str, stratum, counts: dict, key) -> None:
    """Add `table`'s consecutive pairs that `a` and `b` both phase to `counts[(key, group)]`,
    as `[units, wrong units, pairs, wrong pairs]`."""
    indel = [len(r) != 1 or len(t) != 1 for _, r, t in table]
    units = {}
    for i in range(len(table) - 1):
        if None in (a[i], a[i + 1], b[i], b[i + 1]):
            continue
        p, q = table[i][0], table[i + 1][0]
        group = "indel" if indel[i] or indel[i + 1] else stratum(chrom, p, q)
        bad = (a[i] == a[i + 1]) != (b[i] == b[i + 1])
        c = counts.setdefault((key, group), [0, 0, 0, 0])
        c[2] += 1
        c[3] += bad
        units[(group, p // gate.UNIT)] = units.get((group, p // gate.UNIT), False) or bad
    for (group, _), bad in units.items():
        c = counts[(key, group)]
        c[0] += 1
        c[1] += bad


def print_counts(counts: dict) -> None:
    print("comparison                            group       units  wrong_units  unit_rate  "
          "upper95    pairs  wrong_pairs")
    for (key, group), (n, k, p, w) in sorted(counts.items()):
        print(f"{key:<37} {group:<10} {n:>6} {k:>12}  {k / max(1, n):>9.2e}  "
              f"{gate.upper(k, n):>7.2e} {p:>8} {w:>12}")


def decide(chroms, stratum) -> dict:
    """Print the gate and return `{group: [testable betas]}`."""
    extractor = common._module("extract_fragments.py")
    score = common._module("score_hg002.py")
    counts, labelled = {}, [0, 0]
    with tempfile.TemporaryDirectory() as tmp:
        for chrom in chroms:
            out = Path(tmp) / f"{chrom}.vcf.gz"
            extractor.het_snvs(str(common.PEDIGREE / "CEPH1463.GRCh38.family-truthset.ov.vcf.gz"),
                               "NA12878", chrom, str(common.PEDIGREE / "hq_regions_final.bed"),
                               str(out), indels=True, unphase=True)
            table = score.variant_table(out)
            truth = common.pedigree_truth(chrom, table, "NA12878")
            kg = common.kg_transmission(chrom, table, "NA12878")
            labelled[0] += len(table)
            labelled[1] += sum(k is not None for k in kg)
            count(table, kg, truth, chrom, stratum, counts, "NA12878 transmission vs Pedigree")
            print(f"  NA12878 {chrom}: {len(table)} variants", flush=True)
    print(f"\nNA12878: {labelled[1]} of {labelled[0]} variants labelled by transmission "
          f"({labelled[1] / labelled[0]:.3f})\n")
    print_counts(counts)
    print("\ngate (registered in amendment 1):")
    testable = {}
    for group in GROUPS:
        n, k = counts.get(("NA12878 transmission vs Pedigree", group), [0, 0])[:2]
        u = gate.upper(k, n)
        testable[group] = [b for b in gate.BETAS if u <= b / 3]
        print(f"  {group:<10} units {n:>6} wrong {k:>4} upper95 {u:.2e}  testable at beta "
              + (", ".join(f"{b:g}" for b in testable[group]) or "none"))
    return testable


def describe(chroms, stratum) -> None:
    counts = {}
    with tempfile.TemporaryDirectory() as tmp:
        for sample, names in gate.GENOMES.items():
            calls = gate.dipcall_calls(sample)
            other = [n for n in names if n != "transmission"][0]
            for chrom in chroms:
                table = gate.variants(sample, chrom, Path(tmp))
                dip = [calls[(chrom, *v)] for v in table]
                trans = gate.reference(sample, "transmission", chrom, table)
                consensus = [d if t == d else None for d, t in zip(dip, trans)]
                count(table, consensus, gate.reference(sample, other, chrom, table), chrom,
                      stratum, counts, f"{sample} consensus vs {other}")
            print(f"  {sample} done", flush=True)
        for sample in sorted(p.name for p in common.HPRC.iterdir()):
            if not (common.HPRC / sample / f"{sample}.dip.vcf.gz").exists() \
                    or not common.kg_parents(sample):
                continue
            calls = gate.dipcall_calls(sample)
            for chrom in chroms:
                table = gate.variants(sample, chrom, Path(tmp))
                dip = [calls[(chrom, *v)] for v in table]
                count(table, dip, common.kg_transmission(chrom, table, sample), chrom, stratum,
                      counts, f"{sample} dipcall vs 1kGP transmission")
            print(f"  {sample} done", flush=True)
    print()
    print_counts(counts)


def main() -> None:
    from stratified_gate import strata

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--describe", action="store_true")
    args = parser.parse_args()
    chroms = list(common.CONFIRMATORY)
    stratum = strata(chroms)
    print("W1 amendment 1, " + ("descriptive comparisons" if args.describe else "gate") +
          f": consecutive pairs both truths phase; {gate.UNIT // 1000} kb units; one-sided 95% "
          "exact upper bound; testable when at most beta/3", flush=True)
    (describe if args.describe else decide)(chroms, stratum)


if __name__ == "__main__":
    main()
