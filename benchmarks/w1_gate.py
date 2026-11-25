"""W1's gate (`docs/design/w1-ancestry-plan.md`): dipcall truths against trio transmission truths.

For HG005 (release 2 assemblies) and HG002 (year 1 assemblies), the variants
are dipcall's heterozygous calls as W1's phaser input selects them
(`extract_fragments.het_snvs` inside dipcall's BED). A pair is two
consecutive variants of that list whose two variants the reference truth
phases; it disagrees when the two truths give it different relative phases.
Pairs are grouped as A1 and A2 do and counted in 100 kb units. A (group,
beta) is testable when the one-sided 95 per cent exact upper bound on the
share of disagreeing units, pooled over both genomes, is at most beta/3.

The decision uses the transmission truths. The comparisons with v5.0q (HG002)
and with GIAB's hifiasm phase transfer (HG005) are reported, not used.
No phaser output is read.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/w1_gate.py`
"""

from __future__ import annotations

import gzip
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import genome_common as common  # noqa: E402

UNIT = 100_000
BETAS = (3e-4, 1e-3, 3e-3)
GROUPS = ("easy", "difficult", "lowmap", "segdup", "tandem", "indel")
GENOMES = {
    "HG005": {"transmission": ("giab_hg005", "truth/HG005_{c}_trio_truth.vcf.gz", "PATMAT"),
              "hifiasm phase transfer": ("giab_hg005", "split/HG005_{c}_hifiasm_phasetransfer.vcf.gz",
                                         "PATMAT")},
    "HG002": {"transmission": ("giab_wg", "truth/HG002_{c}_trio_truth.vcf.gz", "PATMAT"),
              "v5.0q": None},
}


def dipcall_calls(sample: str) -> dict:
    """`{(chrom, pos1, ref, alt): first allele}` of dipcall's `0|1` and `1|0` calls."""
    calls = {}
    with gzip.open(common.HPRC / sample / f"{sample}.dip.vcf.gz", "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.split("\t", 10)
            gt = f[9].split(":")[0]
            if gt in ("0|1", "1|0"):
                calls[(f[0], int(f[1]), f[3], f[4])] = int(gt[0])
    return calls


def variants(sample: str, chrom: str, tmp: Path) -> list:
    """W1's phaser input for `sample` on `chrom`: `[(pos1, ref, alt)]`."""
    extractor = common._module("extract_fragments.py")
    score = common._module("score_hg002.py")
    out = tmp / f"{sample}_{chrom}.vcf.gz"
    extractor.het_snvs(str(common.HPRC / sample / f"{sample}.dip.vcf.gz"), "syndip", chrom,
                       str(common.HPRC / sample / f"{sample}.dip.bed"), str(out),
                       indels=True, unphase=True)
    return score.variant_table(out)


def reference(sample: str, name: str, chrom: str, table: list) -> list:
    score = common._module("score_hg002.py")
    spec = GENOMES[sample][name]
    if spec is None:  # v5.0q, through the function the runs use (the default platform is HG002's)
        return common.v5_truth(chrom, table)
    folder, pattern, phase_set = spec
    return score.parental_truth(common.ROOT / "data" / folder / pattern.format(c=chrom), table,
                                phase_set=phase_set)


def upper(k: int, n: int) -> float:
    from scipy.stats import beta

    return 1.0 if n == 0 else (1.0 if k == n else float(beta.ppf(0.95, k + 1, n - k)))


def main() -> None:
    from stratified_gate import strata

    chroms = list(common.CONFIRMATORY)
    stratum = strata(chroms)
    counts = {}  # (sample, reference, group) -> [units, wrong units, pairs, wrong pairs]
    shared = {}  # (sample, reference) -> [dipcall variants, phased by the reference]
    print("W1 gate: dipcall truths against reference truths; consecutive pairs of dipcall's "
          f"heterozygous calls inside its BED; {UNIT // 1000} kb units; testable when the "
          "one-sided 95% upper bound on disagreeing units, pooled, is at most beta/3")
    with tempfile.TemporaryDirectory() as tmp:
        for sample in GENOMES:
            calls = dipcall_calls(sample)
            for chrom in chroms:
                table = variants(sample, chrom, Path(tmp))
                dip = [calls[(chrom, *v)] for v in table]
                indel = [len(r) != 1 or len(a) != 1 for _, r, a in table]
                for name in GENOMES[sample]:
                    ref = reference(sample, name, chrom, table)
                    s = shared.setdefault((sample, name), [0, 0])
                    s[0] += len(table)
                    s[1] += sum(t is not None for t in ref)
                    units = {}
                    for i in range(len(table) - 1):
                        if ref[i] is None or ref[i + 1] is None:
                            continue
                        a, b = table[i][0], table[i + 1][0]
                        group = "indel" if indel[i] or indel[i + 1] else stratum(chrom, a, b)
                        bad = (dip[i] == dip[i + 1]) != (ref[i] == ref[i + 1])
                        c = counts.setdefault((sample, name, group), [0, 0, 0, 0])
                        c[2] += 1
                        c[3] += bad
                        key = (group, a // UNIT)
                        units[key] = units.get(key, False) or bad
                    for (group, _), bad in units.items():
                        c = counts[(sample, name, group)]
                        c[0] += 1
                        c[1] += bad
                print(f"  {sample} {chrom}: {len(table)} dipcall variants", flush=True)
    print("\nsample  reference               phased_share  group       units  wrong_units"
          "  unit_rate  upper95    pairs  wrong_pairs")
    for (sample, name, group), (n, k, p, w) in sorted(counts.items()):
        s = shared[(sample, name)]
        print(f"{sample:<7} {name:<23} {s[1] / s[0]:>12.3f}  {group:<10} {n:>6} {k:>12}"
              f"  {k / max(1, n):>9.2e}  {upper(k, n):>7.2e}  {p:>8} {w:>12}")
    print("\npooled over HG002 and HG005, transmission truth:")
    for group in GROUPS:
        n = sum(counts.get((s, "transmission", group), [0, 0])[0] for s in GENOMES)
        k = sum(counts.get((s, "transmission", group), [0, 0])[1] for s in GENOMES)
        u = upper(k, n)
        ok = [b for b in BETAS if u <= b / 3]
        print(f"  {group:<10} units {n:>6} wrong {k:>4} upper95 {u:.2e}  testable at beta "
              + (", ".join(f"{b:g}" for b in ok) if ok else "none"))


if __name__ == "__main__":
    main()
