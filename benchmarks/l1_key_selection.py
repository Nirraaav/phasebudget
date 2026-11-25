"""L1 (`docs/design/l1-key-selection.md`): does the consensus key select easy links?

Step A: HG002 HiFi rows (v5.0q wrong), split by whether an HPRC-style key (HPRC dipcall agreeing
with GIAB trio transmission) labels both variants. Step B: W1b's B1 cells rescored on all links,
on links a family key labels (child's hifiasm dipcall agreeing with transmission from the
parents' 1000 Genomes panel genotypes), and on a matched random subset. Step C: locus sharing
and a batch check, descriptive.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/l1_key_selection.py`
"""

from __future__ import annotations

import argparse
import json
import os
import random
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

ROOT = HERE.parent
KEYS = ROOT / "data" / "l1" / "keys"
CHILDREN = ("NA12879", "NA12881", "NA12882", "NA12885", "NA12886")
FATHER, MOTHER = "NA12877", "NA12878"
CAL = "hprc_HG002"
HARD = ("segdup", "lowmap", "tandem")
GROUPS = ("easy", "segdup", "lowmap", "tandem", "difficult", "indel")
B1_CELLS = 540
LEVEL = 0.05 / B1_CELLS
REPRODUCE = 163


def keyed(labels) -> set:
    """Positions where dipcall and transmission both label and give the same paternal allele."""
    return {p for p, d, t in labels if d is not None and d == t}


def pair_keyed(key: dict, a: int, b: int) -> bool:
    """Family key on one link: the assembly relates the pair (Q6 contig rule) and its relation
    equals the one transmission gives. The children's assemblies are not trio-binned, so only
    relations, not first alleles, can be compared."""
    from q6_assembly_phase import relation

    asm, tra = key["asm"], key["tra"]
    r = relation(asm.get(a), asm.get(b))
    ta, tb = tra.get(a), tra.get(b)
    return r is not None and ta is not None and tb is not None and r == int(ta == tb)


def restrict(rows, keys: dict) -> list:
    """Rows `(group, unit, score, wrong, chrom, a, b)` the key labels: both variants for a
    per-variant key (a set per chromosome), the link for a family key (a dict per chromosome)."""
    def ok(r):
        k = keys[r[4]]
        return pair_keyed(k, r[5], r[6]) if isinstance(k, dict) else r[5] in k and r[6] in k
    return [r for r in rows if ok(r)]


def matched_random(rows, kept, seed: int = 0) -> list:
    """Per group, keep each row with probability the key's kept share of that group."""
    out = []
    for g in {r[0] for r in rows}:
        mine = [r for r in rows if r[0] == g]
        share = sum(r[0] == g for r in kept) / len(mine)
        rng = random.Random(seed)
        out += [r for r in mine if rng.random() < share]
    return out


def family_transmission(chrom: str, variants) -> list:
    """Paternal allele of a heterozygous child of NA12877 and NA12878 from the panel genotypes."""
    from pysam import bcftools

    import genome_common as common

    cache = KEYS / f"family_{chrom}.json"
    if not cache.exists():
        panel = (common.KG_PANEL /
                 f"1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
        text = bcftools.view("-H", "-s", f"{FATHER},{MOTHER}", str(panel))
        found = []
        for line in text.splitlines():
            f = line.split("\t")
            x = common.transmitted(1, common._alt_count(f[9]), common._alt_count(f[10]))
            if x is not None:
                found.append([int(f[1]), f[3], f[4], x])
        part = cache.with_name(f"{cache.name}.{os.getpid()}.part")
        part.write_text(json.dumps(found))
        part.replace(cache)
    labels = {(p, r, a): x for p, r, a, x in json.loads(cache.read_text())}
    return [labels.get(tuple(v)) for v in variants]


def key_labels(chrom: str) -> list:
    """`[pos1, dipcall, transmission]` per input variant, in the current platform's environment."""
    import genome_common as common

    sc = common._module("score_hg002.py")
    variants = sc.variant_table(common.hetsnv_path(chrom))
    if common.SAMPLE == "HG002":
        d = common.dipcall_truth(chrom, variants)
        t = sc.parental_truth(common.truth_path(chrom), variants)
    else:  # assembly values come from K1b's cached tables (`k1b_score.read_table`)
        d = [None] * len(variants)
        t = family_transmission(chrom, variants)
    return [[v[0], a, b] for v, a, b in zip(variants, d, t)]


def load_labels(job) -> list:
    platform, chrom = job
    path = KEYS / f"{platform}_{chrom}.json"
    if not path.exists():
        KEYS.mkdir(parents=True, exist_ok=True)
        part = path.with_name(f"{path.name}.{os.getpid()}.part")
        env = dict(os.environ, PHASEBUDGET_PLATFORM=platform, PYTHONPATH=f"{ROOT / 'src'}:{HERE}")
        subprocess.run([sys.executable, __file__, "--labels", chrom, str(part)], env=env,
                       check=True)
        part.replace(path)
    return json.loads(path.read_text())


def rows_with_pos(platform, depth, phaser, stratum) -> list:
    """`(group, unit, score, wrong, chrom, a, b)` for every cached row of `platform`."""
    import confidence_score as cs
    import genome_common as common
    import stratified_score as ss

    chroms = list(common.CONFIRMATORY)
    ind = ss.indels(platform, chroms)
    return [("indel" if (c, a) in ind or (c, b) in ind else stratum(c, a, b), (c, a // cs.UNIT),
             key, wrong, c, a, b)
            for c, a, key, wrong, b, _ in cs.load_rows(platform, depth, phaser, chroms)]


def rate(rows) -> str:
    k = sum(r[3] for r in rows)
    return f"{k}/{len(rows)} = {k / max(1, len(rows)):.2e}"


def step_a(item) -> list:
    depth, phaser, keys = item
    import genome_common as common
    from stratified_gate import strata

    rows = rows_with_pos(CAL, depth, phaser, strata(list(common.CONFIRMATORY)))
    kept = restrict(rows, keys["key"])
    out = []
    for name, groups in [(g, (g,)) for g in GROUPS] + [("hard", HARD)]:
        mine = [r for r in rows if r[0] in groups]
        lab = [r for r in kept if r[0] in groups]
        rest = [r for r in mine if not (r[5] in keys["key"][r[4]] and r[6] in keys["key"][r[4]])]
        dip = restrict(mine, keys["dip"])
        tra = restrict(mine, keys["tra"])
        k_all, k_lab = sum(r[3] for r in mine), sum(r[3] for r in lab)
        r_all = k_all / max(1, len(mine))
        r_lab = k_lab / max(1, len(lab))
        ratio = r_all / r_lab if k_lab else float("inf") if k_all >= 3 else float("nan")
        out.append((depth, phaser, name, ratio,
                    f"  {depth}x {phaser:<9} {name:<9} links {len(mine):>8} key share "
                    f"{len(lab) / max(1, len(mine)):.4f} (dipcall {len(dip) / max(1, len(mine)):.4f}, "
                    f"transmission {len(tra) / max(1, len(mine)):.4f}); wrong all {rate(mine)}, key "
                    f"{rate(lab)}, rest {rate(rest)}; ratio {ratio:.2f}"))
    return out


def step_b(item) -> dict:
    depth, phaser, keys = item
    import confidence_score as cs
    import genome_common as common
    from stratified_gate import strata
    from stratified_score import cells, grouped
    from w1b_score import POOL, pooled

    from phasebudget.confidence import learn_then_test_by_group

    chroms = list(common.CONFIRMATORY)
    stratum = strata(chroms)
    pool = {"HG002": grouped(CAL, depth, phaser, chroms, stratum),
            **{s: grouped(f"hprc_{s}", depth, phaser, chroms, stratum) for s in POOL}}
    tests = {c: rows_with_pos(f"phifi_{c}", depth, phaser, stratum) for c in CHILDREN}
    kept = {c: restrict(r, keys[c]) for c, r in tests.items()}
    rand = {c: matched_random(r, kept[c]) for c, r in tests.items()}
    res = {"cells": [], "lines": []}
    for beta in cs.BETAS:
        th = learn_then_test_by_group(pooled(pool), beta, cs.DELTA)
        for c in CHILDREN:
            for mode, rows in (("all", tests[c]), ("key", kept[c]), ("random", rand[c])):
                for g, cell in cells(th, [r[:4] for r in rows], beta).items():
                    bad = bool(cell[0]) and cell[2] < LEVEL
                    res["cells"].append((mode, c, depth, phaser, beta, g, bad))
                    res["lines"].append(
                        f"  {mode:<6} {c} {depth}x {phaser:<9} beta {beta:<6g} {g:<9} units "
                        f"{cell[0]:>6} wrong {cell[1]:>4} p {cell[2]:.2e} pairs {cell[3]:>8} trusted "
                        f"{cell[4]:>8} wrong trusted {cell[6]:>5}{'  VIOLATES' if bad else ''}")
    if depth == "30":
        res["lines"] += step_c(phaser, pool, tests, kept)
    print(f"  step B {depth}x {phaser} done", flush=True)
    return res


def step_c(phaser, pool, tests, kept) -> list:
    out = [f"  step C {phaser} 30x: locus sharing of wrong units (chrom, 100 kb)"]
    wrong = {name: {} for name in list(pool) + list(tests)}
    for name, rows in list(pool.items()) + list(tests.items()):
        for r in rows:
            if r[3]:
                wrong[name].setdefault(r[0], set()).add(r[1])
    for g in GROUPS:
        def shared(names, others_of):
            hit = tot = 0
            for n in names:
                mine = wrong[n].get(g, set())
                other = set().union(*(wrong[o].get(g, set()) for o in others_of(n)))
                tot += len(mine)
                hit += len(mine & other)
            return f"{hit}/{tot} = {hit / max(1, tot):.3f}"
        kids, pl = list(tests), list(pool)
        out.append(f"    {g:<9} child wrong units also wrong in another child "
                   f"{shared(kids, lambda n: [k for k in kids if k != n])}; in a pool genome "
                   f"{shared(kids, lambda n: pl)}; pool genome in another pool genome "
                   f"{shared(pl, lambda n: [p for p in pl if p != n])}")
    out.append(f"  step C {phaser} 30x: batch, wrong rate on key-labelled links, children "
               "(family key, Pedigree truth) against the 20 HPRC pool genomes (consensus truth)")
    for g in GROUPS:
        kids = [r for rows in kept.values() for r in rows if r[0] == g]
        hprc = [r for n, rows in pool.items() if n != "HG002" for r in rows if r[0] == g]
        out.append(f"    {g:<9} children {rate(kids)}; HPRC pool {rate(hprc)}")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", nargs=2, metavar=("CHROM", "OUT"))
    args = parser.parse_args()
    if args.labels:
        Path(args.labels[1]).write_text(json.dumps(key_labels(args.labels[0])))
        return

    import confidence_score as cs
    import genome_common as common

    started = time.time()
    chroms = list(common.CONFIRMATORY)
    print(f"L1 key selection (docs/design/l1-key-selection.md, development, rules stated "
          f"before scoring). Rows: cached confidence_score rows; groups {GROUPS}; unit 100 kb; "
          f"betas {cs.BETAS}; delta {cs.DELTA}; depths {cs.DEPTHS}; phasers {cs.PHASERS}; "
          f"chromosomes {len(chroms)} confirmatory.\nStep A: {CAL}, wrong against v5.0q; key = "
          f"HPRC dipcall (data/hprc/HG002) == GIAB trio transmission (giab_wg/truth PATMAT); "
          f"rule: selection ratio over {HARD} > 3 for both phasers at 30x.\nStep B: W1b pooled "
          f"thresholds (HG002 v5.0q + 20 consensus-truth HPRC genomes), tests "
          f"{', '.join('phifi_' + c for c in CHILDREN)} against the Pedigree truth; family key "
          f"(amendment 1) = a link the child's hifiasm assembly relates (Q6 contig rule, K1b "
          f"tables data/k1/tables) with the relation transmission from {FATHER}, {MOTHER} "
          f"1000 Genomes panel genotypes gives; level 0.05/{B1_CELLS} = {LEVEL:.3e}; reproduction "
          f"must give {REPRODUCE}; rule: key <= {int(0.2 * REPRODUCE)} and random >= "
          f"{(REPRODUCE + 1) // 2}; random control seed 0 per cell.", flush=True)

    platforms = [CAL] + [f"phifi_{c}" for c in CHILDREN]
    labels = {}
    for batch in (platforms[:2], platforms[2:]):  # the first child fills the family cache
        jobs = [(p, c) for p in batch for c in chroms]
        with ThreadPoolExecutor(24) as ex:
            labels.update(zip(jobs, ex.map(load_labels, jobs)))
    keys = {"key": {c: keyed(labels[(CAL, c)]) for c in chroms},
            "dip": {c: {p for p, d, _ in labels[(CAL, c)] if d is not None} for c in chroms},
            "tra": {c: {p for p, _, t in labels[(CAL, c)] if t is not None} for c in chroms}}
    fam = {}
    for k in CHILDREN:
        fam[k] = {}
        for c in chroms:
            table = json.loads((ROOT / "data" / "k1" / "tables" /
                                f"phifi_{k}_30_{c}.json").read_text())
            fam[k][c] = {"asm": {v[0]: x for v, x in zip(table["variants"], table["assembly"])
                                 if x is not None},
                         "tra": {p: t for p, _, t in labels[(f"phifi_{k}", c)] if t is not None}}
        n = sum(len(labels[(f"phifi_{k}", c)]) for c in chroms)
        print(f"  family key {k}: of {n} input variants, assembly places "
              f"{sum(len(fam[k][c]['asm']) for c in chroms)}, transmission labels "
              f"{sum(len(fam[k][c]['tra']) for c in chroms)}", flush=True)

    print("\n== Step A")
    items = [(d, ph, keys) for d in cs.DEPTHS for ph in cs.PHASERS]
    with ProcessPoolExecutor(len(items)) as ex:
        a = [r for part in ex.map(step_a, items) for r in part]
    for r in a:
        print(r[4])
    hard30 = {r[1]: r[3] for r in a if r[0] == "30" and r[2] == "hard"}
    worst = min((r for r in a if r[0] == "30" and r[2] in HARD), key=lambda r: r[3])
    rule_a = all(v > 3 for v in hard30.values())
    print(f"Rule A: hard ratio 30x {', '.join(f'{k} {v:.2f}' for k, v in hard30.items())}; worst "
          f"single hard group {worst[1]} {worst[2]} {worst[3]:.2f} -> "
          f"{'supported' if rule_a else 'not supported'}", flush=True)

    print("\n== Step B and C")
    items = [(d, ph, fam) for d in cs.DEPTHS for ph in cs.PHASERS]
    with ProcessPoolExecutor(len(items)) as ex:
        b = list(ex.map(step_b, items))
    for part in b:
        print("\n".join(part["lines"]))
    cells_ = [c for part in b for c in part["cells"]]
    count = {m: sum(c[6] for c in cells_ if c[0] == m) for m in ("all", "key", "random")}
    print(f"\nviolating cells: all {count['all']}, key {count['key']}, random {count['random']}")
    for m in ("all", "key", "random"):
        by = {}
        for c in cells_:
            if c[0] == m and c[6]:
                by[c[5]] = by.get(c[5], 0) + 1
        print(f"  {m}: by group {by}")
    if count["all"] != REPRODUCE:
        print(f"Rule B: reproduction gives {count['all']}, not {REPRODUCE}: Step B stops")
    else:
        ok = count["key"] <= int(0.2 * REPRODUCE) and count["random"] >= (REPRODUCE + 1) // 2
        print(f"Rule B: {'supported' if ok else 'not supported'}")
    print(f"runtime {time.time() - started:.0f} s")


if __name__ == "__main__":
    main()
