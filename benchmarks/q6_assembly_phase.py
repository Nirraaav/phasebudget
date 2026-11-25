"""Round 6, Q6 (development): phase relations of a single-sample HiFi assembly in the hard groups.

`docs/design/q6-assembly-phase.md` states the question, the contig rule and the gate.
The assembly is hifiasm 0.25.0 on HG002's HiFi reads alone (`data/r6_q6/`), called
against GRCh38 by dipcall 0.3 (`data/r6_q6/dip/`). The read-based side is W4's
`w4_grch38` platform at 30x. Each chromosome's table (variants, both truths, both
phasers, the assembly) is read once in a subprocess with that platform and cached.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/q6_assembly_phase.py`
"""

from __future__ import annotations

import argparse
import bisect
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

ROOT = HERE.parent
ASM = ROOT / "data" / "r6_q6"
CACHE = ASM / "tables"
PLATFORM = "w4_grch38"
DEPTH = "30"
CHROMS = tuple(f"chr{i}" for i in range(1, 23) if i != 20)
SOURCES = ("whatshap", "longphase", "assembly")
HARD = ("segdup", "lowmap", "tandem", "indel")
GATE_A = ("segdup", "lowmap", "tandem")
MIN_PAIRS = 20


def segments(intervals) -> tuple[list, list, list]:
    """Disjoint `(starts, ends, names)` from `[(start, end, name)]` (0-based, half open): a
    segment's name is the one alignment covering it, None where two or more do."""
    events = sorted({p for s, e, _ in intervals for p in (s, e)})
    starts, ends, names = [], [], []
    for s, e in zip(events, events[1:]):
        cover = [n for a, b, n in intervals if a <= s and e <= b]
        if cover:
            starts.append(s)
            ends.append(e)
            names.append(cover[0] if len(cover) == 1 else None)
    return starts, ends, names


def contig_at(segs, pos0: int):
    """The one contig covering 0-based `pos0`, None if none or several do."""
    starts, ends, names = segs
    k = bisect.bisect_right(starts, pos0) - 1
    return names[k] if k >= 0 and pos0 < ends[k] else None


def contig_of(segs, pos1: int, ref: str):
    """The one contig covering a variant's whole reference span, None otherwise."""
    first = contig_at(segs, pos1 - 1)
    return first if first is not None and contig_at(segs, pos1 - 2 + len(ref)) == first else None


def assembly_values(variants, calls: dict, inside, hap1, hap2) -> list:
    """`[allele, block]` per variant from dipcall's `calls` `{(pos1, ref, alt): GT}`: only
    `0|1` or `1|0` at the same alleles, inside dipcall's BED (`inside(pos1, ref)`), with the
    block naming the hap1 and hap2 contigs that cover the variant; None otherwise."""
    out = []
    for pos1, ref, alt in variants:
        gt = calls.get((pos1, ref, alt))
        c1, c2 = contig_of(hap1, pos1, ref), contig_of(hap2, pos1, ref)
        ok = gt in ("0|1", "1|0") and inside(pos1, ref) and c1 is not None and c2 is not None
        out.append([int(gt[0]), f"{c1}|{c2}"] if ok else None)
    return out


def relation(value_a, value_b):
    """1 if the source puts the two variants on one haplotype, 0 if not, None if it does not
    relate them. A truth's value is its first allele; a phaser's or the assembly's is
    `[allele, block]`, and two blocks that differ relate nothing."""
    if value_a is None or value_b is None:
        return None
    if isinstance(value_a, list):
        if value_a[1] is None or value_a[1] != value_b[1]:
            return None
        return int(value_a[0] == value_b[0])
    return int(value_a == value_b)


def read_alignments(bam: Path) -> dict:
    """`{chrom: segments}` of a dipcall haplotype BAM's alignments."""
    import pysam

    by_chrom: dict = {}
    with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as handle:
        for r in handle.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary:
                continue
            by_chrom.setdefault(r.reference_name, []).append(
                (r.reference_start, r.reference_end, r.query_name))
    return {c: segments(iv) for c, iv in by_chrom.items()}


def read_table(chrom: str) -> dict:
    """One chromosome's table in `w4_grch38`'s environment."""
    import genome_common as common
    from phase_confidence import phasing_from_vcf

    sc = common._module("score_hg002.py")
    variants = sc.variant_table(common.hetsnv_path(chrom))
    positions = [p - 1 for p, _, _ in variants]
    table: dict = {"variants": [list(v) for v in variants],
                   "v5": common.v5_truth(chrom, variants),
                   # The per-chromosome trio files have no tabix index: read them whole (R1).
                   "trio": sc.parental_truth(common.GIAB / "truth"
                                             / f"HG002_{chrom}_trio_truth.vcf.gz", variants)}
    for phaser in ("whatshap", "longphase"):
        vcf = common.chrom_dir(chrom) / "external" / f"{phaser}_{DEPTH}.vcf"
        vcf = vcf if vcf.exists() else vcf.with_name(vcf.name + ".gz")
        phasing = phasing_from_vcf(str(vcf), common.VCF_SAMPLE, chrom, positions)
        table[phaser] = [None if a is None else [a, b]
                         for a, b in zip(phasing.alleles, phasing.blocks)]
    table["assembly"] = dipcall_values(chrom, variants, ASM / "dip", "HG002")
    return table


def dipcall_values(chrom: str, variants, dip: Path, name: str) -> list:
    """`assembly_values` from dipcall's outputs `<dip>/<name>.{dip.bed,dip.vcf.gz,hap1.bam,
    hap2.bam}` on one chromosome."""
    import pysam

    starts, ends = [], []
    with open(dip / f"{name}.dip.bed") as handle:
        for line in handle:
            f = line.split("\t")
            if f[0] == chrom:
                starts.append(int(f[1]))
                ends.append(int(f[2]))

    def inside(pos1: int, ref: str) -> bool:
        k = bisect.bisect_right(starts, pos1 - 1) - 1
        return k >= 0 and pos1 - 2 + len(ref) < ends[k]

    calls = {}
    with pysam.VariantFile(str(dip / f"{name}.dip.vcf.gz")) as vcf:
        for rec in vcf.fetch(chrom) if vcf.index is not None else vcf:
            if rec.chrom != chrom or rec.alts is None or len(rec.alts) != 1:
                continue
            gt = rec.samples[0]["GT"]
            if None in gt or not rec.samples[0].phased:
                continue
            calls[(rec.pos, rec.ref, rec.alts[0])] = f"{gt[0]}|{gt[1]}"
    hap1 = read_alignments(dip / f"{name}.hap1.bam").get(chrom, ([], [], []))
    hap2 = read_alignments(dip / f"{name}.hap2.bam").get(chrom, ([], [], []))
    return assembly_values(variants, calls, inside, hap1, hap2)


def load_table(chrom: str) -> dict:
    path = CACHE / f"{chrom}.json"
    if not path.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        env = dict(os.environ, PHASEBUDGET_PLATFORM=PLATFORM, PYTHONPATH=f"{ROOT / 'src'}:{HERE}")
        subprocess.run([sys.executable, __file__, "--table", chrom, str(path)], env=env,
                       check=True)
    return json.loads(path.read_text())


def pairs_of(chrom: str, table: dict, stratum) -> list:
    """`(chrom, pos a, pos b, group, {source: relation}, v5 relation, trio relation)` for each
    pair of consecutive input variants."""
    v = table["variants"]
    out = []
    for i in range(len(v) - 1):
        (a, ra, aa), (b, rb, ab) = v[i], v[i + 1]
        indel = len(ra) != 1 or len(aa) != 1 or len(rb) != 1 or len(ab) != 1
        group = "indel" if indel else stratum(chrom, a, b)
        rel = {s: relation(table[s][i], table[s][i + 1]) for s in SOURCES}
        out.append((chrom, a, b, group, rel, relation(table["v5"][i], table["v5"][i + 1]),
                    relation(table["trio"][i], table["trio"][i + 1])))
    return out


def gate_a(pairs, truth: int) -> dict:
    """Per group of `GATE_A`: `(n, right)` over pairs whatshap phases wrongly and the
    assembly phases, against truth column `truth` (5 v5.0q, 6 trio)."""
    out = {}
    for g in GATE_A:
        mine = [p for p in pairs if p[3] == g and p[truth] is not None
                and p[4]["whatshap"] is not None and p[4]["whatshap"] != p[truth]
                and p[4]["assembly"] is not None]
        out[g] = (len(mine), sum(p[4]["assembly"] == p[truth] for p in mine))
    return out


def gate_b(pairs, truth: int) -> dict:
    """Per hard group: `(n, whatshap wrong, longphase wrong, assembly wrong)` on pairs all three
    sources phase, against truth column `truth`."""
    out = {}
    for g in HARD:
        mine = [p for p in pairs if p[3] == g and p[truth] is not None
                and all(p[4][s] is not None for s in SOURCES)]
        out[g] = (len(mine),) + tuple(sum(p[4][s] != p[truth] for p in mine) for s in SOURCES)
    return out


def a_holds(result: dict) -> tuple[bool, list]:
    """(a): in every judged group (at least `MIN_PAIRS` pairs) the assembly is right on at least
    0.8, and at least one group is judged."""
    judged = [g for g, (n, _) in result.items() if n >= MIN_PAIRS]
    return bool(judged) and all(result[g][1] >= 0.8 * result[g][0] for g in judged), judged


def b_holds(result: dict) -> tuple[bool, list]:
    """(b): in every judged group (at least `MIN_PAIRS` whatshap errors) the assembly's errors
    are at most a third of whatshap's, and at least one group is judged."""
    judged = [g for g, (_, w, _, _) in result.items() if w >= MIN_PAIRS]
    return bool(judged) and all(3 * result[g][3] <= result[g][1] for g in judged), judged


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", nargs=2, metavar=("CHROM", "OUT"))
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    if args.table:
        Path(args.table[1]).write_text(json.dumps(read_table(args.table[0])))
        return

    import confidence_score as cs
    import stratified_score as ss
    from stratified_gate import strata

    from phasebudget.confidence import learn_then_test_by_group

    print(f"Q6 (development, docs/design/q6-assembly-phase.md): HG002 HiFi; assembly hifiasm "
          f"0.25.0-r726 default options on data/w4/reads.fq.gz, dipcall 0.3 against GRCh38 no-alt "
          f"(data/r6_q6/dip); read-based {PLATFORM} {DEPTH}x whatshap and LongPhase; pairs of "
          f"consecutive v5.0q heterozygous variants inside its BED; chromosomes "
          f"{' '.join(CHROMS)}; truths v5.0q and GIAB trio transmission (PATMAT); assembly relation "
          f"only within one hap1 contig and one hap2 contig; certificates frozen score, per group "
          f"on odd chromosomes, betas {cs.BETAS}, delta {cs.DELTA}, {cs.UNIT // 1000} kb units; "
          f"gate groups (a) {GATE_A}, (b) {HARD}, at least {MIN_PAIRS}", flush=True)
    stratum = strata(list(CHROMS))
    with ThreadPoolExecutor(args.workers) as pool:
        tables = dict(zip(CHROMS, pool.map(load_table, CHROMS)))
    pairs = [p for c in CHROMS for p in pairs_of(c, tables[c], stratum)]
    even = [p for p in pairs if not int(p[0][3:]) % 2]

    # Recheck: the cached W4 rows and this table must agree on every pair both hold.
    rows = {}
    for phaser in ("whatshap", "longphase"):
        rows[phaser] = {(c, a, b): (key, wrong)
                        for c, a, key, wrong, b, _ in cs.load_rows(PLATFORM, DEPTH, phaser,
                                                                    list(CHROMS))}
        both = [p for p in pairs if (p[0], p[1], p[2]) in rows[phaser]
                and p[4][phaser] is not None and p[5] is not None]
        bad = sum(rows[phaser][(p[0], p[1], p[2])][1] != (p[4][phaser] != p[5]) for p in both)
        print(f"check {phaser}: {len(both)} pairs in both the W4 rows and the table, "
              f"{bad} disagree on wrong; rows {len(rows[phaser])}")

    print("\nPairs each source phases, and errors against v5.0q (all 21 chromosomes)")
    print(f"  {'group':<10} {'pairs':>9} {'v5 rel':>9} | " + " | ".join(
        f"{s:>9} phased   wrong" for s in SOURCES))
    for g in ("easy", "difficult") + HARD:
        mine = [p for p in pairs if p[3] == g]
        cols = []
        for s in SOURCES:
            ph = [p for p in mine if p[4][s] is not None and p[5] is not None]
            cols.append(f"{len(ph):>16} {sum(p[4][s] != p[5] for p in ph):>7}")
        print(f"  {g:<10} {len(mine):>9} {sum(p[5] is not None for p in mine):>9} | "
              + " | ".join(cols))

    verdicts = {}
    for label, truth in (("v5.0q", 5), ("trio", 6)):
        ra, rb = gate_a(even, truth), gate_b(pairs, truth)
        a, ja = a_holds(ra)
        b, jb = b_holds(rb)
        verdicts[label] = (a, b)
        print(f"\nAgainst {label}")
        print("  (a) even chromosomes, pairs whatshap phases wrongly and the assembly phases: "
              "assembly right")
        for g, (n, right) in ra.items():
            print(f"    {g:<10} {right}/{n} = {right / max(1, n):.3f}"
                  f"{'' if g in ja else ' (not judged)'}")
        print("  (b) all chromosomes, pairs all three phase: wrong whatshap / LongPhase / assembly")
        for g, (n, w, lp, asm) in rb.items():
            print(f"    {g:<10} n {n:>8}  {w:>6} {lp:>6} {asm:>6}  rates {w / max(1, n):.2e} "
                  f"{lp / max(1, n):.2e} {asm / max(1, n):.2e}{'' if g in jb else ' (not judged)'}")
        print(f"  (a) {'ok' if a else 'FAIL'} (judged {ja}); (b) {'ok' if b else 'FAIL'} "
              f"(judged {jb})")

    print("\nHybrid (descriptive, v5.0q, even chromosomes): per phaser and beta, pairs the "
          "certificate does not trust and the assembly phases, and errors there")
    for phaser in ("whatshap", "longphase"):
        odd = [c for c in CHROMS if int(c[3:]) % 2]
        cal = ss.grouped(PLATFORM, DEPTH, phaser, odd, stratum)
        for beta in cs.BETAS:
            t = learn_then_test_by_group(cal, beta, cs.DELTA)
            for g in ("easy", "difficult") + HARD:
                mine = [p for p in even if p[3] == g and p[5] is not None
                        and (p[0], p[1], p[2]) in rows[phaser] and p[4][phaser] is not None]
                thr = t.get(g)
                untrusted = [p for p in mine
                             if thr is None or rows[phaser][(p[0], p[1], p[2])][0] < thr]
                both = [p for p in untrusted if p[4]["assembly"] is not None]
                print(f"  {phaser:<9} {beta:<7g} {g:<10} pairs {len(mine):>8} untrusted "
                      f"{len(untrusted):>8} wrong {sum(p[4][phaser] != p[5] for p in untrusted):>6}"
                      f"; assembly phases {len(both):>8}: {phaser} wrong "
                      f"{sum(p[4][phaser] != p[5] for p in both):>6}, assembly wrong "
                      f"{sum(p[4]['assembly'] != p[5] for p in both):>6}")

    c = all(verdicts["trio"])
    a, b = verdicts["v5.0q"]
    print(f"\nQ6 gate: (a) {'ok' if a else 'FAIL'}; (b) {'ok' if b else 'FAIL'}; (c) "
          f"{'ok' if c else 'FAIL'} (against the trio truth)")
    print(f"Q6 gate: {'passed' if a and b and c else 'not passed'}")


if __name__ == "__main__":
    main()
