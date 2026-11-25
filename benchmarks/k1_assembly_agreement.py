"""Round 7, K1a and K4 (development): certificates that require a single-sample assembly to agree.

`docs/design/k1-assembly-agreement.md` states the certificate and both gates. Inputs are Q6's:
the assembly tables of `benchmarks/q6_assembly_phase.py`, W4's cached `w4_grch38` 30x rows and,
for K4, Q4's chain files for `w4_grch38` (`benchmarks/q4_clinical_cis_trans.py`).

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/k1_assembly_agreement.py`
"""

from __future__ import annotations

import argparse
import gzip
import json
import multiprocessing
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import q6_assembly_phase as q6  # noqa: E402

ROOT = HERE.parent
PLATFORM, DEPTH, CHROMS = q6.PLATFORM, q6.DEPTH, q6.CHROMS
ODD = tuple(c for c in CHROMS if int(c[3:]) % 2)
EVEN = tuple(c for c in CHROMS if not int(c[3:]) % 2)
PHASERS = ("whatshap", "longphase")
GROUPS = ("easy", "difficult") + q6.HARD
GATE_GROUPS = ("segdup", "tandem", "indel")
GATE_BETA = 1e-3
MIN_TRUST = 0.5
K4_TARGET, K4_SHARE = 1e-2, 0.5
LOW = (-1.0, -1)  # below every threshold: the score of a link the assembly does not agree with

# Filled in `main` before the K4 pool forks.
ENTRIES: dict = {}
ASSEMBLY: dict = {}


# ---- pure logic, tested in tests/test_k1_assembly_agreement.py ----------------------------

def agrees(phaser_relation, assembly_relation) -> bool:
    """True when the assembly relates the pair and gives the phaser's relation."""
    return (phaser_relation is not None and assembly_relation is not None
            and phaser_relation == assembly_relation)


def entries(pairs, rows: dict, phaser: str) -> list:
    """`(chrom, pos a, group, unit, score, wrong, agree, both wrong, trio wrong)` per Q6 pair that
    `phaser` relates in `rows` (`{(chrom, a, b): (score, wrong)}`) and v5.0q phases; `trio
    wrong` is None where the trio truth does not phase the pair."""
    import confidence_score as cs

    out = []
    for chrom, a, b, group, rel, v5, trio in pairs:
        hit = rows.get((chrom, a, b))
        if hit is None or rel[phaser] is None or v5 is None:
            continue
        score, wrong = hit
        asm = rel["assembly"]
        out.append((chrom, a, group, (chrom, a // cs.UNIT), score, wrong,
                    agrees(rel[phaser], asm), bool(wrong) and asm is not None and asm != v5,
                    None if trio is None else rel[phaser] != trio))
    return out


def rewrite_links(raw: list, assembly: dict) -> list:
    """Q4 chain entries `[pos1, allele, truth, score, indel]` with the score of each link the
    assembly does not agree with set to `LOW`. A link joins a phased variant to the next phased
    one; `assembly` maps pos1 to Q6's `[allele, block]`."""
    out = [list(v) for v in raw]
    phased = [k for k, v in enumerate(raw) if v[1] is not None]
    for k, nxt in zip(phased, phased[1:]):
        if raw[k][3] is None:
            continue
        mine = int(raw[k][1] == raw[nxt][1])
        if not agrees(mine, q6.relation(assembly.get(raw[k][0]), assembly.get(raw[nxt][0]))):
            out[k][3] = list(LOW)
    return out


# ---- K4: Q4's call with agreeing links ------------------------------------------------------

def k4_cell(job: tuple) -> dict:
    """Q4's crossfit tallies for one phaser, with (`agree_only`) or without the agreement rule."""
    phaser, agree_only = job
    import confidence_score as cs
    import q4_clinical_cis_trans as q4
    from stratified_gate import strata

    stratum = strata(list(CHROMS))
    genes = q4.clinical_genes(list(CHROMS))
    nt = len(q4.TARGETS)
    total = {"pairs": 0, "chain": 0, "certified": [0] * nt, "wrong": [0] * nt}
    chosen: dict = {}
    for cal, test in ((ODD, EVEN), (EVEN, ODD)):
        rows = [(e[2], e[3], e[4], e[5]) for e in ENTRIES[phaser]
                if e[0] in cal and (e[6] or not agree_only)]
        thresholds = q4.thresholds_for(rows)
        for chrom in test:
            with gzip.open(q4.chain_path(ROOT, PLATFORM, DEPTH, phaser, chrom), "rt") as handle:
                raw = json.load(handle)
            if agree_only:
                raw = rewrite_links(raw, ASSEMBLY[chrom])
            phased = [k for k, v in enumerate(raw) if v[1] is not None]
            nxt = dict(zip(phased, phased[1:]))
            variants = []
            for k, (pos, allele, truth, key, indel) in enumerate(raw):
                link = None
                if key is not None:
                    other = raw[nxt[k]]
                    group = "indel" if indel or other[4] else stratum(chrom, pos, other[0])
                    t = thresholds.get(group, [None] * len(q4.BETAS))
                    link = (q4.link_level(tuple(key), t), (group, pos // cs.UNIT), group)
                variants.append((pos, allele, truth, link))
            reach, holding = q4.gene_index(variants, genes[chrom])
            for i, j, status, cells, level, _, _ in q4.pairs(variants, reach, q4.WINDOW):
                total["pairs"] += 1
                if status != "chain":
                    continue
                total["chain"] += 1
                wrong = ((variants[i][1] ^ variants[i][2])
                         != (variants[j][1] ^ variants[j][2]))
                pos_j = variants[j][0]
                names = [g for e, g in holding[i] if e >= pos_j]
                for x, r in enumerate(q4.TARGETS):
                    b = q4.b_max(len(cells), r)
                    if b is None or level > b:
                        continue
                    total["certified"][x] += 1
                    total["wrong"][x] += wrong
                    for name in names:
                        crc = zlib.crc32(f"{name}:{variants[i][0]}:{pos_j}".encode())
                        old = chosen.get((name, x))
                        if old is None or crc < old[0]:
                            chosen[(name, x)] = (crc, wrong)
    validity = []
    for x in range(nt):
        picks = [w for (_, y), (_, w) in chosen.items() if y == x]
        validity.append((len(picks), sum(picks)))
    return {"phaser": phaser, "agree_only": agree_only, "total": total, "validity": validity}


# ---- report ---------------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    import confidence_score as cs
    import q4_clinical_cis_trans as q4
    import stratified_score as ss
    from stratified_gate import strata

    from phasebudget.confidence import learn_then_test_by_group

    started = time.time()
    print(f"K1a and K4 (development, docs/design/k1-assembly-agreement.md): HG002 HPRC HiFi, "
          f"{PLATFORM} {DEPTH}x, whatshap and LongPhase; assembly hifiasm 0.25.0-r726 HiFi only, "
          f"dipcall 0.3 (Q6 tables, data/r6_q6/tables); pairs of consecutive v5.0q heterozygous "
          f"variants inside its BED, truth v5.0q (trio transmission descriptive); agreement = the "
          f"assembly relates the pair (Q6 contig rule) with the phaser's relation; K1a: "
          f"learn_then_test_by_group on agreeing pairs of odd chromosomes {' '.join(ODD)}, test on "
          f"even {' '.join(EVEN)}, betas {cs.BETAS}, delta {cs.DELTA}, {cs.UNIT // 1000} kb units; "
          f"gate (a) no violating cell at 0.05/cells, (b) {', '.join(GATE_GROUPS)} each certify at "
          f"{GATE_BETA:g} for one phaser trusting >= {MIN_TRUST} of the group's even pairs; K4: "
          f"Q4 call, betas {q4.BETAS}, targets {q4.TARGETS}, window {q4.WINDOW}, two folds by "
          f"parity, gate (a) no cell below 0.05/4, (b) certified share >= {K4_SHARE} at r "
          f"{K4_TARGET:g} for one phaser", flush=True)

    stratum = strata(list(CHROMS))
    with ThreadPoolExecutor(6) as pool:
        tables = dict(zip(CHROMS, pool.map(q6.load_table, CHROMS)))
    pairs = [p for c in CHROMS for p in q6.pairs_of(c, tables[c], stratum)]
    for c in CHROMS:
        ASSEMBLY[c] = {v[0]: a for v, a in zip(tables[c]["variants"], tables[c]["assembly"])}
    del tables
    for phaser in PHASERS:
        rows = {(c, a, b): (key, wrong) for c, a, key, wrong, b, _ in
                cs.load_rows(PLATFORM, DEPTH, phaser, list(CHROMS))}
        ENTRIES[phaser] = entries(pairs, rows, phaser)
    del pairs

    # 1. Agreement by group.
    print("\n1. Pairs each phaser relates (v5.0q phases), share the assembly agrees with, and "
          "wrong pairs among agreeing and other pairs (all 21 chromosomes)")
    for phaser in PHASERS:
        for g in GROUPS:
            mine = [e for e in ENTRIES[phaser] if e[2] == g]
            agree = [e for e in mine if e[6]]
            rest = [e for e in mine if not e[6]]
            print(f"  {phaser:<9} {g:<9} pairs {len(mine):>8} agree {len(agree) / max(1, len(mine)):.4f}"
                  f"  wrong among agreeing {sum(e[5] for e in agree):>5} of {len(agree):>8}"
                  f"  among others {sum(e[5] for e in rest):>5} of {len(rest):>7}")

    # 2. Both wrong: what no threshold on agreeing pairs can remove.
    print("\n2. Pairs wrong in both the phaser and the assembly (all 21 chromosomes): pairs, "
          "units of the group, units holding one, share (above 1e-3 the group cannot certify "
          "at 1e-3)")
    for phaser in PHASERS:
        for g in GROUPS:
            mine = [e for e in ENTRIES[phaser] if e[2] == g]
            units = {e[3] for e in mine}
            bad = {e[3] for e in mine if e[7]}
            print(f"  {phaser:<9} {g:<9} both wrong {sum(e[7] for e in mine):>5}  units "
                  f"{len(units):>6} holding one {len(bad):>4} share {len(bad) / max(1, len(units)):.2e}")

    # 3. Certificates: frozen score on all pairs against agreement, odd to even.
    print("\n3. Certificates, calibrated on odd and tested on even chromosomes. all = frozen score "
          "on all pairs; agree = frozen score on agreeing pairs only. Cell: threshold, units with "
          "a trusted pair, units with a wrong trusted pair, P(X >= k | n, beta), trusted pairs, "
          "trusted share of the group's even pairs, wrong trusted")
    cells, results = [], {}
    for phaser in PHASERS:
        cal = [e for e in ENTRIES[phaser] if e[0] in ODD]
        test = [e for e in ENTRIES[phaser] if e[0] in EVEN]
        for beta in cs.BETAS:
            for mode in ("all", "agree"):
                keep = (lambda e: True) if mode == "all" else (lambda e: e[6])
                thresholds = learn_then_test_by_group(
                    [(e[2], e[3], e[4], e[5]) for e in cal if keep(e)], beta, cs.DELTA)
                res = ss.cells(thresholds, [(e[2], e[3], e[4], e[5]) for e in test if keep(e)],
                               beta)
                for g in GROUPS:
                    if g not in res:
                        continue
                    n, k, p, _, kept, _, wrong_kept, _ = res[g]
                    total = sum(e[2] == g for e in test)
                    t = thresholds.get(g)
                    share = kept / max(1, total)
                    results[(phaser, beta, mode, g)] = (t is not None, share, n, k, p)
                    if mode == "agree" and t is not None and n:
                        cells.append((phaser, beta, g, n, k, p))
                    print(f"  {phaser:<9} {beta:<6g} {mode:<5} {g:<9} "
                          f"{'none' if t is None else f'({t[0]:.4f}, {t[1]})':<16} units {n:>6} "
                          f"wrong {k:>3} p {p:.2e} trusted {kept:>8} share {share:.4f} "
                          f"wrong trusted {wrong_kept:>4}")

    level = 0.05 / max(1, len(cells))
    violating = [c for c in cells if c[5] < level]
    ok_a = not violating
    print(f"\n  agreement cells with a threshold and a unit: {len(cells)}; level {level:.2e}; "
          f"violating {len(violating)}"
          + "".join(f"\n    {ph} {b:g} {g}: {k} of {n} units, p {p:.2e}"
                    for ph, b, g, n, k, p in violating))
    ok_b = True
    for g in GATE_GROUPS:
        best = [(ph, results.get((ph, GATE_BETA, "agree", g), (False, 0.0))) for ph in PHASERS]
        good = [ph for ph, (cert, share, *_) in best if cert and share >= MIN_TRUST]
        ok_b &= bool(good)
        print(f"  (b) {g:<7} at {GATE_BETA:g}: " + ", ".join(
            f"{ph} {'certifies' if cert else 'none'} share {share:.4f}"
            for ph, (cert, share, *_) in best) + f": {'holds' if good else 'fails'}")

    # 4. Trio truth, descriptive.
    print("\n4. Against GIAB trio transmission (descriptive): agreeing pairs it phases, and phaser "
          "wrong among them (all 21 chromosomes)")
    for phaser in PHASERS:
        for g in GROUPS:
            mine = [e for e in ENTRIES[phaser] if e[2] == g and e[6] and e[8] is not None]
            print(f"  {phaser:<9} {g:<9} agreeing, trio-phased {len(mine):>8} wrong "
                  f"{sum(e[8] for e in mine):>4}")
    print(f"\nK1a gate: (a) {'holds' if ok_a else 'fails'}; (b) {'holds' if ok_b else 'fails'}")
    print(f"K1a gate: {'passed' if ok_a and ok_b else 'not passed'}", flush=True)

    # 5. K4.
    jobs = [(ROOT, PLATFORM, DEPTH, ph, c) for ph in PHASERS for c in CHROMS]
    with ProcessPoolExecutor(args.workers, mp_context=multiprocessing.get_context("fork")) as pool:
        list(pool.map(q4._extract_job, jobs))
        k4 = list(pool.map(k4_cell, [(ph, mode) for ph in PHASERS for mode in (False, True)]))
    print("\n5. K4: within-gene pairs (Q4 population) on w4_grch38 30x, two folds. Per phaser and "
          "rule: pairs, chains, then per target certified, share, wrong, genes picked, wrong "
          "picks, P(X >= k | n, r)")
    k4_cells = []
    for res in k4:
        t = res["total"]
        rule = "agree" if res["agree_only"] else "all"
        for x, r in enumerate(q4.TARGETS):
            n, k = res["validity"][x]
            p = q4.p_exceed(k, n, r)
            share = t["certified"][x] / max(1, t["pairs"])
            if res["agree_only"]:
                k4_cells.append((res["phaser"], r, share, p))
            print(f"  {res['phaser']:<9} {rule:<5} pairs {t['pairs']:>8} chains {t['chain']:>8}"
                  f"  r {r:<6g} certified {t['certified'][x]:>8} share {share:.4f} wrong "
                  f"{t['wrong'][x]:>4} genes {n:>4} wrong {k:>2} p {p:.2e}")
    k4_a = all(p >= 0.05 / 4 for *_, p in k4_cells)
    k4_b = any(share >= K4_SHARE for ph, r, share, _ in k4_cells if r == K4_TARGET)
    print(f"\nK4 gate: (a) {'holds' if k4_a else 'fails'}; (b) {'holds' if k4_b else 'fails'}")
    print(f"K4 gate: {'passed' if k4_a and k4_b else 'not passed'}")
    print(f"runtime {time.time() - started:.0f} s")


if __name__ == "__main__":
    main()
