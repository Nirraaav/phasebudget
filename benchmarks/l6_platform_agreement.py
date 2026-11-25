"""L6 (development, `docs/design/l6-platform-agreement.md`): HiFi and ONT agreement.

For each person phased on both platforms, every adjacent pair one platform's phaser relates
(cached rows of `confidence_score.load_rows`, scored against that platform's truth) is checked
against the other platform's phasing by the same phaser, matched by (chrom, pos, ref, alt).
`agree` is K1b's rule with the second platform in place of the assembly. The gate is decided on
`segdup`; step 2 (K1b-style certificates on agreeing pairs) runs only if the gate passes.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/l6_platform_agreement.py`
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

ROOT = HERE.parent
DEPTH = "30"
CHROMS = tuple(f"chr{i}" for i in range(1, 23) if i != 20)
PHASERS = ("whatshap", "longphase")
GROUPS = ("indel", "segdup", "lowmap", "tandem", "difficult", "easy")
PEDIGREE = ("NA12879", "NA12881", "NA12882", "NA12885", "NA12886")
PEOPLE = {"HG002": ("hprc_HG002", "cal_v5"),
          **{s: (f"phifi_{s}", f"ped_{s}") for s in PEDIGREE},
          "HG005": ("k1_hg005", "hg005_ont")}
POOL = ("HG002", *PEDIGREE)
GATE_SHARE, GATE_MIN_WRONG, GATE_AGREE, GATE_MEMBERS = 0.20, 10, 0.5, 3
CELLS = (len(POOL) + 1) * len(PHASERS) * 3 * len(GROUPS)  # 252
LEVEL = 0.05 / CELLS
FLOOR_BETA, MIN_TRUST = 1e-3, 0.5
TABLES = ROOT / "data" / "l6" / "tables"


def read_table(chrom: str) -> dict:
    """Variants and both phasers' `[allele, block]` (None if unphased), current platform."""
    import genome_common as common
    from phase_confidence import phasing_from_vcf

    sc = common._module("score_hg002.py")
    variants = sc.variant_table(common.hetsnv_path(chrom))
    positions = [p - 1 for p, _, _ in variants]
    table: dict = {"variants": [list(v) for v in variants]}
    for phaser in PHASERS:
        vcf = common.chrom_dir(chrom) / "external" / f"{phaser}_{DEPTH}.vcf"
        vcf = vcf if vcf.exists() else vcf.with_name(vcf.name + ".gz")
        phasing = phasing_from_vcf(str(vcf), common.VCF_SAMPLE, chrom, positions)
        table[phaser] = [None if a is None else [a, b]
                         for a, b in zip(phasing.alleles, phasing.blocks)]
    return table


def load_table(job: tuple) -> dict:
    platform, chrom = job
    path = TABLES / f"{platform}_{DEPTH}_{chrom}.json"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        env = dict(os.environ, PHASEBUDGET_PLATFORM=platform,
                   PYTHONPATH=f"{ROOT / 'src'}:{HERE}")
        part = path.with_name(path.name + ".part")
        subprocess.run([sys.executable, __file__, "--table", chrom, str(part)], env=env,
                       check=True)
        part.replace(path)
    return json.loads(path.read_text())


def entries(label: str, rows, own: dict, other: dict, phaser: str, stratum,
            unit: int = 100_000) -> list:
    """`(group, unit, score, wrong, agree, relates, matched)` per row `(chrom, a, score, wrong, b,
    _)`. `own` and `other` map chrom to a table; `other`'s values are found by (pos, ref, alt)."""
    import k1_assembly_agreement as k1
    import q6_assembly_phase as q6

    by_pos = {c: {} for c in own}
    for c, t in own.items():
        for v in t["variants"]:
            by_pos[c].setdefault(v[0], v)
    by_key = {c: {tuple(v): val for v, val in zip(t["variants"], t[phaser])}
              for c, t in other.items()}
    own_val = {c: {tuple(v): val for v, val in zip(t["variants"], t[phaser])}
               for c, t in own.items()}
    out = []
    for c, a, key, wrong, b, _ in rows:
        if c not in own:
            continue
        va, vb = tuple(by_pos[c][a]), tuple(by_pos[c][b])
        indel = any(len(x) != 1 for x in (va[1], va[2], vb[1], vb[2]))
        group = "indel" if indel else stratum(c, a, b)
        theirs = by_key.get(c, {})
        matched = va in theirs and vb in theirs
        mine = q6.relation(own_val[c][va], own_val[c][vb])
        check = q6.relation(theirs.get(va), theirs.get(vb))
        out.append((group, (label, c, a // unit), key, bool(wrong), k1.agrees(mine, check),
                    check is not None, matched))
    return out


def summary(mine: list) -> dict:
    """Counts for one (person, phaser, direction, group)."""
    n = len(mine)
    wrong = [r for r in mine if r[3]]
    judged = [r for r in wrong if r[5]]
    shared = sum(r[4] for r in judged)
    agree = [r for r in mine if r[4]]
    return {"pairs": n, "matched": sum(r[6] for r in mine), "relates": sum(r[5] for r in mine),
            "agree": len(agree), "wrong": len(wrong), "W": len(judged), "S": shared,
            "share": shared / len(judged) if judged else None,
            "wrong_agree": sum(r[3] for r in agree)}


def gate(stats: dict) -> tuple[bool, list]:
    """The design note's gate on (person, phaser) -> segdup summary, HiFi checked by ONT."""
    def passes(s):
        return (s["W"] >= GATE_MIN_WRONG and s["share"] < GATE_SHARE
                and s["agree"] >= GATE_AGREE * max(1, s["pairs"]))

    lines, ok = [], False
    for phaser in PHASERS:
        hg002 = passes(stats[("HG002", phaser)])
        members = sum(passes(stats[(s, phaser)]) for s in PEDIGREE)
        good = hg002 and members >= GATE_MEMBERS
        ok |= good
        lines.append(f"  {phaser}: HG002 {'passes' if hg002 else 'fails'}, Pedigree members "
                     f"passing {members} of {len(PEDIGREE)} -> {'passes' if good else 'fails'}")
    return ok, lines


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", nargs=2, metavar=("CHROM", "OUT"))
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    if args.table:
        Path(args.table[1]).write_text(json.dumps(read_table(args.table[0])))
        return

    import confidence_score as cs
    import stratified_score as ss
    from stratified_gate import strata

    from phasebudget.confidence import learn_then_test_by_group

    started = time.time()
    print(f"L6 (development, docs/design/l6-platform-agreement.md): {DEPTH}x, phasers {PHASERS}; "
          f"people {', '.join(f'{p} (HiFi {h}, ONT {o})' for p, (h, o) in PEOPLE.items())}; "
          f"{len(CHROMS)} chromosomes; pairs matched by (chrom, pos, ref, alt); agree = K1b rule "
          f"with the other platform's same phaser in place of the assembly; gate: segdup shared "
          f"wrong share < {GATE_SHARE} with W >= {GATE_MIN_WRONG} and agree share >= {GATE_AGREE} "
          f"on HG002 and >= {GATE_MEMBERS} of {len(PEDIGREE)} Pedigree members, one phaser; "
          f"step 2: learn_then_test_by_group on pooled agreeing HiFi pairs, unit (person, chrom, "
          f"{cs.UNIT // 1000} kb), betas {cs.BETAS}, delta {cs.DELTA}, {CELLS} cells, level "
          f"{LEVEL:.3e}", flush=True)

    stratum = strata(list(CHROMS))
    jobs = sorted({(p, c) for pair in PEOPLE.values() for p in pair for c in CHROMS})
    with ThreadPoolExecutor(args.workers) as ex:
        loaded = dict(zip(jobs, ex.map(load_table, jobs)))
    tables = {p: {c: loaded[(p, c)] for c in CHROMS} for pair in PEOPLE.values() for p in pair}

    hifi_entries: dict = {}
    segdup: dict = {}
    print("\nper person, phaser, direction, group: pairs, matched share, relates share, agree "
          "share, wrong, W (wrong the other platform relates), S (shared wrong), S/W, wrong rate "
          "among agreeing, wrong rate among all")
    for person, (hifi, ont) in PEOPLE.items():
        for phaser in PHASERS:
            for direction, (own, other) in (("HiFi by ONT", (hifi, ont)),
                                            ("ONT by HiFi", (ont, hifi))):
                rows = cs.load_rows(own, DEPTH, phaser, list(CHROMS))
                got = entries(person, rows, tables[own], tables[other], phaser, stratum)
                if direction == "HiFi by ONT":
                    hifi_entries[(person, phaser)] = got
                for g in GROUPS:
                    s = summary([r for r in got if r[0] == g])
                    if direction == "HiFi by ONT" and g == "segdup":
                        segdup[(person, phaser)] = s
                    n = max(1, s["pairs"])
                    share = "n/a" if s["share"] is None else f"{s['share']:.3f}"
                    print(f"  {person:<8} {phaser:<9} {direction:<11} {g:<9} pairs {s['pairs']:>8}"
                          f" matched {s['matched'] / n:.4f} relates {s['relates'] / n:.4f} agree "
                          f"{s['agree'] / n:.4f} wrong {s['wrong']:>6} W {s['W']:>6} S {s['S']:>5}"
                          f" S/W {share:>5} rate agree {s['wrong_agree'] / max(1, s['agree']):.2e}"
                          f" rate all {s['wrong'] / n:.2e}", flush=True)

    ok, lines = gate(segdup)
    print("\nGate (segdup, HiFi checked by ONT):")
    print("\n".join(lines))
    worst = max((s["share"] for s in segdup.values() if s["share"] is not None), default=None)
    print(f"  worst segdup S/W over all people and phasers: {worst}")
    print(f"Gate {'PASSES' if ok else 'FAILS'}")
    if not ok:
        print(f"runtime {time.time() - started:.0f} s")
        return

    print("\nStep 2: certificates on agreeing HiFi pairs")
    valid, useful = True, False
    for phaser in PHASERS:
        def cal(people):
            return [r[:4] for p in people for r in hifi_entries[(p, phaser)] if r[4]]

        full = {b: learn_then_test_by_group(cal(POOL), b, cs.DELTA) for b in cs.BETAS}
        seg = full[FLOOR_BETA].get("segdup") is not None
        print(f"  {phaser}: full pool segdup at {FLOOR_BETA:g} "
              f"{'certifies' if seg else 'abstains'}")
        tests = [("HG005", full)] + [
            (p, {b: learn_then_test_by_group(cal([q for q in POOL if q != p]), b, cs.DELTA)
                 for b in cs.BETAS}) for p in POOL]
        for person, thresholds in tests:
            got = hifi_entries[(person, phaser)]
            for beta in cs.BETAS:
                res = ss.cells(thresholds[beta], [r[:4] for r in got if r[4]], beta)
                for g in GROUPS:
                    n_all = sum(r[0] == g for r in got)
                    n, k, p, _, kept, _, wk, _ = res.get(g, (0, 0, 1.0, 0, 0, 0, 0, False))
                    bad = n > 0 and p < LEVEL
                    valid &= not bad
                    share = kept / max(1, n_all)
                    if person == "HG005" and g == "segdup" and beta == FLOOR_BETA and seg:
                        useful |= share >= MIN_TRUST
                    print(f"    {person:<8} {phaser:<9} {beta:<6g} {g:<9} units {n:>6} wrong {k:>3}"
                          f" p {p:.2e} trusted {kept:>8} share {share:.4f} wrong trusted {wk:>4}"
                          f"{'  VIOLATES' if bad else ''}")
    print(f"L6-2 validity: {'holds' if valid else 'fails'}")
    print(f"L6-3 usefulness: {'holds' if useful else 'fails'}")
    print(f"runtime {time.time() - started:.0f} s")


if __name__ == "__main__":
    main()
