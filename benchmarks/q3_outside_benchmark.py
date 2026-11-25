"""Round 6, Q3 (development): how much phasing error lies outside the v5.0q benchmark.

`docs/design/q3-outside-benchmark.md` fixes the definitions, the gate and the
cause split. Inputs are the stored outputs of `confidence_run.py` for
`q3_q100` (HG002 HiFi, W4's minimap2 GRCh38 alignments, 30x, phased on the Q100
dipcall heterozygous variants inside the Q100 dip BED) and, for calibration,
the cached `w4_grch38` rows (the same reads on v5.0q's variants, scored
against v5.0q inside its BED).

Usage: `PHASEBUDGET_PLATFORM=q3_q100 PYTHONPATH=src:benchmarks .venv/bin/python
benchmarks/q3_outside_benchmark.py`
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Optional, Sequence

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from phasebudget._types import Fragment  # noqa: E402
from phasebudget.fragments import FragmentMatrix  # noqa: E402

PLATFORM = "q3_q100"
CALIBRATION = "w4_grch38"
DEPTH = "30"
PHASERS = ("whatshap", "longphase")
CHROMS = tuple(f"chr{i}" for i in range(1, 23) if i != 20)
UNIT = 100_000
FLANK = 10  # variants kept on each side of the pair in the cause window
CAP = 15  # `whatshap phase --internal-downsampling` default (WhatsHap 2.8)
GATE_RATIO = 3.0
GATE_SHARE = 0.20
GROUPS = ("indel", "segdup", "lowmap", "tandem", "difficult", "easy")
CAUSES = ("none", "data", "search", "other")


# ---- pure logic, tested in tests/test_q3_outside_benchmark.py ----------------

def window_matrix(matrix: FragmentMatrix, a: int, b: int, flank: int = FLANK,
                  candidates: Optional[Sequence[int]] = None
                  ) -> tuple[FragmentMatrix, int, list[int]]:
    """The fragments of `matrix` on columns `a-flank..a` and `b..b+flank`, renumbered.

    Columns strictly between `a` and `b` are dropped, so `a` and `b` are
    adjacent and their relation is the boundary returned. Each read keeps its
    called alleles on the kept columns and its id; a read calling none of them
    is left out. `candidates` (indices into `matrix`) limits the reads looked
    at; it must include every read overlapping the window. Returns
    `(window, boundary, kept columns)`.
    """
    if not 0 <= a < b < matrix.n_variants:
        raise ValueError("need 0 <= a < b < n_variants")
    columns = list(range(max(0, a - flank), a + 1)) + \
        list(range(b, min(matrix.n_variants, b + flank + 1)))
    new = {c: k for k, c in enumerate(columns)}
    lo, hi = columns[0], columns[-1]
    fragments = []
    for index in (range(len(matrix)) if candidates is None else candidates):
        f = matrix[index]
        if f.end <= lo or f.start > hi:
            continue
        calls = [(new[j], allele, weight) for j, allele, weight in f.observed() if j in new]
        if not calls:
            continue
        start = calls[0][0]
        alleles: list = [None] * (calls[-1][0] - start + 1)
        weights = [0.0] * len(alleles)
        for k, allele, weight in calls:
            alleles[k - start] = allele
            weights[k - start] = weight
        fragments.append(Fragment(f.read_id, start, tuple(alleles), tuple(weights)))
    return FragmentMatrix(fragments, len(columns)), new[a], columns


def relation_costs(window: FragmentMatrix, boundary: int, cap: Optional[int] = CAP,
                   positions: Optional[Sequence[int]] = None) -> tuple[float, float]:
    """Exact wMEC optimum with `h[boundary] ^ h[boundary+1]` forced to 0 and to 1.

    With `cap`, reads are first chosen by the library's port of WhatsHap 2.8's
    read selection at that coverage, as `whatshap phase` does by default.
    """
    from phasebudget.dp import CoverageDP
    from phasebudget.selection import select_whatshap

    if cap is not None:
        window = window.subset(select_whatshap(window, cap, positions=positions).selected)
    return CoverageDP(window).boundary_costs(boundary)


def majority(agree: int, disagree: int) -> str:
    """Spanning reads against a wrong call: `none`, or whether most back the wrong relation."""
    if agree + disagree == 0:
        return "none"
    return "wrong" if agree > disagree else "tie" if agree == disagree else "right"


def optimum(costs: tuple[float, float], true_relation: int) -> str:
    """`right`, `wrong` or `tie`: which relation the constrained optima prefer."""
    right, wrong = costs[true_relation], costs[1 - true_relation]
    return "right" if right < wrong else "wrong" if right > wrong else "tie"


def cause(vote: str, best: str) -> str:
    """The cause of a wrong relation (design note): `none` (i), `data` (ii), `search` (iii),
    `other` (spanning reads back the truth but the window optimum does not)."""
    if vote == "none":
        return "none"
    if best == "right":
        return "search"
    return "data" if vote in ("wrong", "tie") else "other"


def unit_shares(rows) -> tuple[int, int, int, int]:
    """`(pairs, wrong pairs, units, wrong units)` for rows `(chrom, pos1, wrong)`; a unit is
    `(chrom, pos1 // UNIT)` holding a pair, wrong when it holds a wrong pair."""
    units, bad = set(), set()
    pairs = wrong = 0
    for chrom, pos1, w in rows:
        pairs += 1
        wrong += bool(w)
        units.add((chrom, pos1 // UNIT))
        if w:
            bad.add((chrom, pos1 // UNIT))
    return pairs, wrong, len(units), len(bad)


def gate(inside: tuple, outside: tuple) -> tuple[float, float, bool]:
    """`(per-pair ratio outside/inside, share of wrong outside, holds)` for one phaser, from
    `unit_shares` tuples."""
    rate_in = inside[1] / max(1, inside[0])
    rate_out = outside[1] / max(1, outside[0])
    ratio = rate_out / rate_in if rate_in else float("inf")
    share = outside[1] / max(1, inside[1] + outside[1])
    return ratio, share, ratio >= GATE_RATIO and share >= GATE_SHARE


# ---- rows from the stored outputs ------------------------------------------------

def rows_path(phaser: str, chrom: str) -> Path:
    import genome_common as common

    return common.ROOT / "data" / "r6_q3" / "rows" / f"{phaser}_{DEPTH}_{chrom}.json"


def build_rows(chrom: str, phaser: str, platform: str = PLATFORM) -> list:
    """Per consecutive phased pair both of whose variants the platform's truth (Q100 here)
    phases: `[pos1, pos2, agree, disagree, score, wrong, indel, costs]`, `costs` the window's
    constrained optima `[relation 0, relation 1, true relation]` for wrong pairs, else None.
    Round 7's K2 calls it for other platforms (`benchmarks/k2_k3_error_causes.py`)."""
    import genome_common as common
    from phase_confidence import phasing_from_vcf
    from phasebudget.confidence import pair_support, score

    if common.PLATFORM != platform:
        raise SystemExit(f"set PHASEBUDGET_PLATFORM={platform}")
    extractor = common._module("extract_fragments.py")
    sc = common._module("score_hg002.py")
    variants = sc.variant_table(common.hetsnv_path(chrom))
    by_pos = {p: (r, a) for p, r, a in variants}
    matrix, positions, _, _ = extractor.load_fragments(
        str(common.fragments_path(chrom, DEPTH, "realign")))
    truth = common.truth_at(chrom, variants, positions)
    vcf = common.chrom_dir(chrom) / "external" / f"{phaser}_{DEPTH}.vcf"
    vcf = vcf if vcf.exists() else vcf.with_name(vcf.name + ".gz")
    phasing = phasing_from_vcf(str(vcf), common.VCF_SAMPLE, chrom, positions)
    order = sorted(range(len(matrix)), key=lambda i: matrix[i].start)
    starts = [matrix[i].start for i in order]
    longest = max((len(f) for f in matrix), default=0)
    out = []
    for s in pair_support(matrix, phasing):
        a, b = s.first, s.second
        if truth[a] is None or truth[b] is None:
            continue
        wrong = (phasing.alleles[a] == phasing.alleles[b]) != (truth[a] == truth[b])
        indel = any(len(x) != 1 for k in (a, b) for x in by_pos[positions[k] + 1])
        costs = None
        if wrong:
            lo, hi = max(0, a - FLANK), min(matrix.n_variants - 1, b + FLANK)
            near = order[bisect.bisect_left(starts, lo - longest):bisect.bisect_right(starts, hi)]
            window, boundary, columns = window_matrix(matrix, a, b, FLANK, near)
            costs = [*relation_costs(window, boundary, CAP, [positions[c] for c in columns]),
                     truth[a] ^ truth[b]]
        out.append([positions[a] + 1, positions[b] + 1, s.agree, s.disagree, list(score(s)),
                    bool(wrong), indel, costs])
    return out


def load_rows(job: tuple[str, str]) -> list:
    phaser, chrom = job
    path = rows_path(phaser, chrom)
    if not path.exists():
        started = time.time()
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_name(path.name + f".{os.getpid()}.part")
        part.write_text(json.dumps(build_rows(chrom, phaser)))
        part.replace(path)
        print(f"  rows {phaser} {chrom} in {time.time() - started:.0f} s", flush=True)
    return json.loads(path.read_text())


# ---- report -----------------------------------------------------------------------

def bed_lookup(path: Path, chroms):
    starts: dict = {c: [] for c in chroms}
    ends: dict = {c: [] for c in chroms}
    with open(path) as handle:
        for line in handle:
            f = line.split("\t")
            if f[0] in starts:
                starts[f[0]].append(int(f[1]))
                ends[f[0]].append(int(f[2]))

    def inside(chrom: str, pos1: int) -> bool:
        k = bisect.bisect_right(starts[chrom], pos1 - 1) - 1
        return k >= 0 and pos1 - 1 < ends[chrom][k]

    return inside


def genotypes(path: Path, chroms, het_only: bool) -> dict:
    """`(chrom, pos1, ref, alt) -> GT` of a VCF's records on `chroms`; with `het_only`, only
    `0|1` and `1|0` records."""
    calls = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line[0] == "#":
                continue
            f = line.split("\t", 10)
            if f[0] not in chroms:
                continue
            gt = f[9].split(":")[0].rstrip("\n")
            if not het_only or gt in ("0|1", "1|0"):
                calls[(f[0], int(f[1]), f[3], f[4])] = gt
    return calls


def truth_check(wrongs: list, q100: dict, v5: dict) -> dict:
    """Counts over wrong relations `(chrom, pos1, pos2, ...)`: does v5.0q call both variants
    as the same heterozygous variant as Q100 (`q100` holds Q100's `0|1`/`1|0` records), and
    if so, does v5.0q's phase give Q100's relation?"""
    at = {key[:2]: (key, gt) for key, gt in q100.items()}
    counts = {"genotype differs": 0, "same genotype, same relation": 0,
              "same genotype, other relation": 0, "same genotype, v5.0q unphased": 0}
    for r in wrongs:
        mine = [at[(r[0], p)] for p in (r[1], r[2])]
        theirs = [v5.get(key) for key, _ in mine]
        if any(t is None or t.replace("|", "/") not in ("0/1", "1/0") for t in theirs):
            counts["genotype differs"] += 1
        elif any("|" not in t for t in theirs):
            counts["same genotype, v5.0q unphased"] += 1
        else:
            same = (theirs[0][0] == theirs[1][0]) == (mine[0][1][0] == mine[1][1][0])
            counts["same genotype, " + ("same" if same else "other") + " relation"] += 1
    return counts


def certificate_cells(cal_rows: list, rows: list, beta: float) -> tuple[dict, dict]:
    """Per group, thresholds from `cal_rows` and `stratified_score.cells` of `rows` under them."""
    import confidence_score as cs
    from stratified_score import cells

    from phasebudget.confidence import learn_then_test_by_group

    thresholds = learn_then_test_by_group(cal_rows, beta, cs.DELTA)
    return thresholds, cells(thresholds, rows, beta)


def main() -> None:
    import confidence_score as cs
    import genome_common as common
    import stratified_score
    from stratified_gate import strata

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--chromosomes", nargs="+", default=list(CHROMS))
    args = parser.parse_args()
    chroms = tuple(args.chromosomes)
    descriptive = chroms != CHROMS
    started = time.time()
    print(f"Q3 (development, docs/design/q3-outside-benchmark.md): HG002 HPRC HiFi, minimap2 "
          f"2.31 GRCh38 (data/w4/HG002_grch38_mm2.bam), {DEPTH}x slice drawn as confidence_run "
          f"(seed 0, trim seed 3); variants Q100 dipcall {common.Q100_VCF.name} biallelic "
          f"heterozygous SNVs and indels <= 50 bp inside {common.Q100_BED.name}, phase removed; "
          f"phasers whatshap phase 2.8 (defaults, --reference, --ignore-read-groups) and "
          f"LongPhase 2.0.2 (--pb --indels -t 1); truth Q100 phase; inside = both variants in "
          f"HG002_GRCh38_v5.0q_smvar.benchmark.bed; groups {', '.join(GROUPS)} (indel first, "
          f"then GIAB v3.6 stratum); unit {UNIT // 1000} kb; gate per-pair ratio >= {GATE_RATIO:g}"
          f" and share outside >= {GATE_SHARE:g} for both phasers; certificates calibrated on "
          f"{CALIBRATION} {DEPTH}x rows (v5.0q, 21 chromosomes) with learn_then_test_by_group, "
          f"betas {cs.BETAS}, delta {cs.DELTA}; cause window {FLANK} variants per side, "
          f"WhatsHap selection port at cap {CAP}, exact CoverageDP boundary costs; chromosomes "
          f"{' '.join(chroms)}" + (" (DESCRIPTIVE, not the gate)" if descriptive else ""),
          flush=True)
    jobs = [(p, c) for p in PHASERS for c in chroms]
    with ProcessPoolExecutor(args.workers) as pool:
        loaded = dict(zip(jobs, pool.map(load_rows, jobs)))
    stratum = strata(list(dict.fromkeys(chroms + CHROMS)))  # calibration reads CHROMS too
    v5_inside = bed_lookup(common.GIAB_V5 / "HG002_GRCh38_v5.0q_smvar.benchmark.bed", chroms)
    rows = {p: [] for p in PHASERS}  # (chrom, p1, p2, agree, disagree, score, wrong, group, in, n_in, costs)
    for (phaser, chrom), table in loaded.items():
        for p1, p2, agree, disagree, key, wrong, indel, costs in table:
            n_in = v5_inside(chrom, p1) + v5_inside(chrom, p2)
            group = "indel" if indel else stratum(chrom, p1, p2)
            rows[phaser].append((chrom, p1, p2, agree, disagree, tuple(key), wrong, group,
                                 n_in == 2, n_in, costs))

    # 1. Error inside and outside, and the gate.
    print("\n1. Consecutive relations against Q100: pairs, wrong, error per pair; 100 kb units, "
          "wrong units, wrong unit share")
    verdicts = {}
    for phaser in PHASERS:
        mine = rows[phaser]
        for group in ("all",) + GROUPS:
            sel = [r for r in mine if group in ("all", r[7])]
            cut = {"inside": [r for r in sel if r[8]], "outside": [r for r in sel if not r[8]],
                   "  one in": [r for r in sel if r[9] == 1],
                   "  none in": [r for r in sel if r[9] == 0]}
            stats = {k: unit_shares([(r[0], r[1], r[6]) for r in v]) for k, v in cut.items()}
            for k, (n, w, u, wu) in stats.items():
                print(f"  {phaser:<9} {group:<9} {k:<9} pairs {n:>9} wrong {w:>6} "
                      f"rate {w / max(1, n):.3e}  units {u:>6} wrong {wu:>5} share "
                      f"{wu / max(1, u):.4f}")
            ratio, share, ok = gate(stats["inside"], stats["outside"])
            unit_ratio = ((stats["outside"][3] / max(1, stats["outside"][2]))
                          / max(1e-12, stats["inside"][3] / max(1, stats["inside"][2])))
            print(f"  {phaser:<9} {group:<9} outside/inside per pair {ratio:.2f}, per unit "
                  f"{unit_ratio:.2f}; share of wrong outside {share:.3f}")
            if group == "all":
                verdicts[phaser] = (ratio, share, ok, unit_ratio)

    # 2. Frozen-score certificates outside the benchmark.
    print(f"\n2. Certificates: thresholds from {CALIBRATION} {DEPTH}x (v5.0q) per group, applied "
          "to Q100-scored relations inside and outside; cell = units, wrong trusted units, "
          "P(X >= k | n, beta), pairs, trusted, wrong, wrong trusted")
    outside_cells = []
    for phaser in PHASERS:
        cal = stratified_score.grouped(CALIBRATION, DEPTH, phaser, list(CHROMS), stratum)
        for beta in cs.BETAS:
            for where in ("inside", "outside"):
                sel = [(r[7], (r[0], r[1] // UNIT), r[5], r[6]) for r in rows[phaser]
                       if r[8] == (where == "inside")]
                thresholds, result = certificate_cells(cal, sel, beta)
                for group in GROUPS:
                    if group not in result:
                        continue
                    c = result[group]
                    t = thresholds.get(group)
                    print(f"  {phaser:<9} beta {beta:<6g} {where:<7} {group:<9} threshold "
                          f"{'none' if t is None else f'({t[0]:.4f}, {t[1]})':<16} units {c[0]:>6}"
                          f" wrong {c[1]:>4} p {c[2]:.2e} pairs {c[3]:>8} trusted {c[4]:>8} "
                          f"wrong {c[5]:>6} wrong trusted {c[6]:>5}")
                    if where == "outside" and c[0]:
                        outside_cells.append((phaser, beta, group, c))
    level = 0.05 / max(1, len(outside_cells))
    bad = [x for x in outside_cells if x[3][2] < level]
    nominal = [x for x in outside_cells if x[3][2] < 0.05]
    print(f"  outside cells with a threshold and a unit: {len(outside_cells)}; violating at "
          f"0.05/{len(outside_cells)} = {level:.2e}: {len(bad)}; at 0.05: {len(nominal)}"
          + "".join(f"\n    {p} beta {b:g} {g}: {c[1]} of {c[0]} units, p {c[2]:.2e}"
                    for p, b, g, c in nominal))

    # 3. Causes of wrong relations.
    print(f"\n3. Causes of wrong relations: none = no read calls both variants (i); data = "
          f"spanning majority backs the wrong relation or ties, window optimum not right (ii); "
          f"search = window optimum (cap {CAP}) strictly right (iii); other = majority right, "
          "optimum wrong or tied")
    for phaser in PHASERS:
        for where in ("inside", "outside"):
            wrongs = [r for r in rows[phaser] if r[6] and r[8] == (where == "inside")]
            for group in ("all",) + GROUPS:
                sel = [r for r in wrongs if group in ("all", r[7])]
                if not sel:
                    continue
                counts = dict.fromkeys(CAUSES, 0)
                cross: dict = {}
                for r in sel:
                    vote = majority(r[3], r[4])
                    best = optimum((r[10][0], r[10][1]), r[10][2])
                    counts[cause(vote, best)] += 1
                    cross[(vote, best)] = cross.get((vote, best), 0) + 1
                print(f"  {phaser:<9} {where:<7} {group:<9} wrong {len(sel):>6}: "
                      + ", ".join(f"{k} {v} ({v / len(sel):.3f})" for k, v in counts.items()))
                if group == "all":
                    print("      majority x optimum: " + ", ".join(
                        f"{v}/{o} {n}" for (v, o), n in sorted(cross.items())))

    # 4. Truth check: v5.0q's genotype at the variants of outside-BED wrong relations.
    print("\n4. Outside-BED wrong relations against v5.0q's VCF (all records, BED ignored): "
          "both variants the same heterozygous call in v5.0q, and then whether v5.0q's phase "
          "gives the same relation as Q100")
    v5 = genotypes(common.GIAB_V5 / "HG002_GRCh38_v5.0q_smvar.vcf.gz", chroms, het_only=False)
    q100 = genotypes(common.Q100_VCF, chroms, het_only=True)
    for phaser in PHASERS:
        wrongs = [r for r in rows[phaser] if r[6] and not r[8]]
        counts = truth_check(wrongs, q100, v5)
        print(f"  {phaser:<9} outside wrong {len(wrongs)}: "
              + ", ".join(f"{k} {v}" for k, v in counts.items()))

    print()
    for phaser, (ratio, share, ok, unit_ratio) in verdicts.items():
        print(f"{phaser}: per pair outside/inside {ratio:.2f} (>= {GATE_RATIO:g}), share of wrong "
              f"outside {share:.3f} (>= {GATE_SHARE:g}): {'holds' if ok else 'fails'}; per unit "
              f"{unit_ratio:.2f} (reported)")
    passed = all(v[2] for v in verdicts.values())
    print(f"Q3 gate{' (descriptive chromosomes)' if descriptive else ''}: "
          f"{'finds a gap' if passed else 'does not find a gap'}")
    print(f"runtime {time.time() - started:.0f} s")


if __name__ == "__main__":
    main()
