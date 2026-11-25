"""Round 4, W4 gate (development): confidently wrong relations on GRCh38 and on CHM13.

`docs/design/w4-reference-gate.md` fixes the definitions, the region list and
the gate. Inputs are the stored outputs of `confidence_run.py` for
`w4_grch38` and `w4_chm13` at 30-fold: fragments, phaser VCFs, and v5.0q on
each reference. A GRCh38 variant is carried to CHM13 by the T2T chain.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/w4_gate.py`
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

import stratify  # noqa: E402

ROOT = HERE.parent
W4 = ROOT / "data" / "w4"
CHAIN = W4 / "grch38-chm13v2.chain"  # md5 913bb9a3352e88c23010d968d8b3af8c
REGIONS = ("GRCh38_false_duplications_correct_copy", "GRCh38_false_duplications_incorrect_copy",
           "GRCh38_collapsed_duplication_FP_regions")
CACHE = ROOT / "data" / "confidence_rows"
CHROMS = tuple(f"chr{i}" for i in range(1, 23) if i != 20)
PHASERS = ("whatshap", "longphase")
MIN_READS = 10
GROUPS = ("segdup", "lowmap", "tandem", "difficult", "easy")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def load_chain(path: Path, chrom: str) -> list:
    """Ungapped blocks `(t_start, t_end, q_start, strand, q_size)` from `chrom` to the same name,
    sorted by `t_start`; coordinates 0-based, `q_start` on the chain's query strand."""
    blocks, keep = [], False
    with open(path) as handle:
        for line in handle:
            f = line.split()
            if not f:
                continue
            if f[0] == "chain":
                keep = f[2] == chrom and f[7] == chrom
                t, q, strand, q_size = int(f[5]), int(f[10]), f[9], int(f[8])
                continue
            if keep:
                size = int(f[0])
                blocks.append((t, t + size, q, strand, q_size))
                if len(f) == 3:
                    t += size + int(f[1])
                    q += size + int(f[2])
    return sorted(blocks)


def lift(blocks: list, starts: list, pos0: int):
    """`(pos0 on the query, strand)` for a 0-based target position, or None in a gap.
    `starts` is `[b[0] for b in blocks]`."""
    k = bisect.bisect_right(starts, pos0) - 1
    if k < 0 or pos0 >= blocks[k][1]:
        return None
    t, _, q, strand, q_size = blocks[k]
    q_pos = q + pos0 - t
    return (q_pos, "+") if strand == "+" else (q_size - 1 - q_pos, "-")


def match(ref: str, alt: str, strand: str, other: tuple):
    """None if `other` (ref, alt) is not the same SNV, else whether its alleles are swapped."""
    if len(ref) != 1 or len(alt) != 1 or len(other[0]) != 1 or len(other[1]) != 1:
        return None
    if strand == "-":
        ref, alt = ref.translate(COMPLEMENT), alt.translate(COMPLEMENT)
    if (ref, alt) == other:
        return False
    if (alt, ref) == other:
        return True
    return None


def relation(table: dict, i: int, j: int):
    """`(wrong, confident)` for variants `i < j` of one table, or None when the phaser has them
    in different blocks or the truth does not phase both."""
    truth, alleles = table["truth"], table["alleles"]
    if truth[i] is None or truth[j] is None or alleles[i] is None or alleles[j] is None:
        return None
    nxt, support = table["next"], table["support"]
    k, confident = i, True
    while k != j:
        if k not in nxt or nxt[k] > j:
            return None
        agree, disagree = support[k]
        confident = confident and disagree == 0 and agree >= MIN_READS
        k = nxt[k]
    wrong = (alleles[i] == alleles[j]) != (truth[i] == truth[j])
    return wrong, confident


def build_table(chrom: str, phaser: str) -> dict:
    """Per variant: position, alleles, phaser allele, truth; per consecutive phased pair its
    support. In the current platform's environment."""
    import genome_common as common
    from phase_confidence import phasing_from_vcf
    from phasebudget.confidence import pair_support

    extractor = common._module("extract_fragments.py")
    sc = common._module("score_hg002.py")
    variants = sc.variant_table(common.hetsnv_path(chrom))
    by_pos = {p: (r, a) for p, r, a in variants}
    matrix, positions, _, _ = extractor.load_fragments(
        str(common.fragments_path(chrom, "30", "realign")))
    vcf = common.chrom_dir(chrom) / "external" / f"{phaser}_30.vcf"
    vcf = vcf if vcf.exists() else vcf.with_name(vcf.name + ".gz")
    phasing = phasing_from_vcf(str(vcf), common.VCF_SAMPLE, chrom, positions)
    pairs = [(s.first, s.second, s.agree, s.disagree) for s in pair_support(matrix, phasing)]
    return {"pos1": [p + 1 for p in positions],
            "alleles_ra": [list(by_pos[p + 1]) for p in positions],
            "alleles": list(phasing.alleles),
            "truth": common.truth_at(chrom, variants, positions),
            "pairs": pairs}


def load_table(platform: str, chrom: str, phaser: str) -> dict:
    path = CACHE / f"w4_{platform}_30_{phaser}_{chrom}.json"
    if not path.exists():
        env = dict(os.environ, PHASEBUDGET_PLATFORM=platform,
                   PYTHONPATH=f"{ROOT / 'src'}:{HERE}")
        subprocess.run([sys.executable, __file__, "--table", chrom, phaser, str(path)],
                       env=env, check=True)
    table = json.loads(path.read_text())
    table["next"] = {a: b for a, b, _, _ in table["pairs"]}
    table["support"] = {a: (g, d) for a, _, g, d in table["pairs"]}
    return table


def regions(chroms) -> list:
    return [stratify.load_bed(W4 / f"{name}.bed.gz", chroms) for name in REGIONS]


def in_regions(beds: list, chrom: str, positions0) -> bool:
    return any(any(stratify.inside(bed[chrom], positions0)) for bed in beds)


def analyse(phaser: str, stratum, beds) -> tuple[list, dict]:
    lines = []
    totals = {ref: [0, 0, 0] for ref in ("GRCh38", "CHM13")}  # relations, wrong, confident wrong
    groups = {g: [0, 0, 0] for g in GROUPS}
    region = [0, 0, 0]
    matched = {"relations": 0, "truths disagree": 0, "cw38": 0, "cw13": 0,
               "cw38 right13": 0, "cw38 right13 confident": 0, "cw38 wrong13": 0,
               "cw38 not one block13": 0}
    for chrom in CHROMS:
        g38 = load_table("w4_grch38", chrom, phaser)
        c13 = load_table("w4_chm13", chrom, phaser)
        blocks = load_chain(CHAIN, chrom)
        starts = [b[0] for b in blocks]
        index13 = {p: k for k, p in enumerate(c13["pos1"])}
        for table, ref in ((g38, "GRCh38"), (c13, "CHM13")):
            for a, b, _, _ in table["pairs"]:
                r = relation(table, a, b)
                if r is None:
                    continue
                totals[ref][0] += 1
                totals[ref][1] += r[0]
                totals[ref][2] += r[0] and r[1]
        for a, b, _, _ in g38["pairs"]:
            r = relation(g38, a, b)
            if r is None:
                continue
            p1, p2 = g38["pos1"][a], g38["pos1"][b]
            cw = r[0] and r[1]
            rows = [groups[stratum(chrom, p1, p2)]]
            if in_regions(beds, chrom, [p1 - 1, p2 - 1]):
                rows.append(region)
            for row in rows:
                row[0] += 1
                row[1] += r[0]
                row[2] += cw
            mapped = []
            for k in (a, b):
                lifted = lift(blocks, starts, g38["pos1"][k] - 1)
                j = None if lifted is None else index13.get(lifted[0] + 1)
                swap = None if j is None else match(*g38["alleles_ra"][k], lifted[1],
                                                    tuple(c13["alleles_ra"][j]))
                mapped.append(None if swap is None else (j, swap))
            if None in mapped or c13["truth"][mapped[0][0]] is None \
                    or c13["truth"][mapped[1][0]] is None:
                continue
            (i13, s1), (j13, s2) = sorted(mapped)
            t38 = g38["truth"][a] == g38["truth"][b]
            t13 = (c13["truth"][i13] == c13["truth"][j13]) != (s1 != s2)
            matched["relations"] += 1
            if t38 != t13:
                matched["truths disagree"] += 1
                continue
            r13 = relation(c13, i13, j13) if i13 != j13 else None
            matched["cw13"] += bool(r13 and r13[0] and r13[1])
            if cw:
                matched["cw38"] += 1
                if r13 is None:
                    matched["cw38 not one block13"] += 1
                elif r13[0]:
                    matched["cw38 wrong13"] += 1
                else:
                    matched["cw38 right13"] += 1
                    matched["cw38 right13 confident"] += r13[1]
    for ref, (n, w, cw) in totals.items():
        lines.append(f"  {phaser:<9} {ref:<7} relations {n:>9} wrong {w:>6} confident wrong {cw:>5}")
    lines.append(f"  {phaser:<9} GRCh38 by group: relations, wrong, confident wrong")
    for g, (n, w, cw) in groups.items():
        lines.append(f"      {g:<10} {n:>9} {w:>6} {cw:>5}")
    n_all = totals["GRCh38"][0]
    cw_all = totals["GRCh38"][2]
    share_rel = region[0] / max(1, n_all)
    share_cw = region[2] / max(1, cw_all)
    lines.append(f"      region list {region[0]:>6} relations ({share_rel:.2e}), wrong {region[1]}, "
                 f"confident wrong {region[2]} ({share_cw:.2e} of all confident wrong)")
    lines.append(f"  {phaser:<9} matched: " + ", ".join(f"{k} {v}" for k, v in matched.items()))
    removal = (matched["cw38"] > 0 and matched["cw13"] <= matched["cw38"] / 2
               and matched["cw38 right13"] >= matched["cw38"] / 2)
    location = region[2] >= 5 and share_cw >= 5 * share_rel
    return lines, {"removal": removal, "location": location}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", nargs=3, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.table:
        chrom, phaser, path = args.table
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path + ".part").write_text(json.dumps(build_table(chrom, phaser)))
        Path(path + ".part").rename(path)
        return
    import stratified_gate

    print(f"W4 gate (development, docs/design/w4-reference-gate.md): HG002 HPRC HiFi realigned "
          f"with minimap2, 30x, chromosomes {' '.join(CHROMS)}; truth v5.0q per reference; "
          f"confident = every consecutive pair has at least {MIN_READS} reads and none disagrees; "
          f"region list {', '.join(REGIONS)} (GIAB v3.6); chain {CHAIN.name}", flush=True)
    jobs = [(p, c, ph) for p in ("w4_grch38", "w4_chm13") for c in CHROMS for ph in PHASERS]
    with ThreadPoolExecutor(8) as pool:
        list(pool.map(lambda job: load_table(*job), jobs))
    stratum = stratified_gate.strata(CHROMS)
    beds = regions(CHROMS)
    verdicts = {}
    for phaser in PHASERS:
        lines, verdicts[phaser] = analyse(phaser, stratum, beds)
        print("\n".join(lines), flush=True)
    removal = all(v["removal"] for v in verdicts.values())
    location = any(v["location"] for v in verdicts.values())
    print("\n" + "; ".join(f"{p}: removal {'ok' if v['removal'] else 'FAIL'}, location "
                           f"{'ok' if v['location'] else 'FAIL'}" for p, v in verdicts.items()))
    print(f"W4 gate: {'passed' if removal and location else 'not passed'}"
          + ("" if removal == location else
             f" (removal {'holds' if removal else 'fails'}, location "
             f"{'holds' if location else 'fails'})"))


if __name__ == "__main__":
    main()
