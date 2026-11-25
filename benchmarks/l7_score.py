"""L7 (registered, `docs/design/l7-na12878-plan.md`): NA12878 HiFi, predictions P1 to P5.

P1 and P2: W1b's pooled thresholds on the test genome's links, all links and links the family
key labels (L1 amendment 1's key with the test genome's own parents), plus L1's random control.
P3: K1b's agreeing-pair certificates. P4: L2's end-variant calls, assembly arm. P5: wrong pairs
inside and outside the GIAB HG001 v4.2.1 benchmark BED at 30x.

`--no-assembly` scores P1 and P5 only. `--dry-run` puts an already-scored child, NA12879 with
parents NA12877 and NA12878, in NA12878's place: it checks the code and reads no NA12878 phase.
NA12879 is in the K1b and L2 pools, so the dry run's P3 and P4 numbers mean nothing.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/l7_score.py [--no-assembly] [--dry-run]`
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

ROOT = HERE.parent
KEYS = ROOT / "data" / "l7" / "keys"
GIAB_BED = ROOT / "data" / "giab_hg001" / "HG001_GRCh38_1_22_v4.2.1_benchmark.bed"
GROUPS = ("easy", "difficult", "segdup", "lowmap", "tandem", "indel")
P1_GROUPS = ("segdup", "lowmap", "indel")
P12_LEVEL = 0.05 / 108
P3_LEVEL, P3_BETA, P3_SHARE = 0.05 / 36, 1e-3, 0.5
P4_R, P4_SHARE, P4_BIN = 1e-2, 0.7, 100_000
P5_SHARE, P5_RATIO = 0.5, 3.0
TEST = {"real": ("NA12878", "phifi_NA12878", ("NA12891", "NA12892")),
        "dry": ("NA12879", "phifi_NA12879", ("NA12877", "NA12878"))}


# ---- pure rules (`tests/test_l7_score.py`) ----

def outside_stats(rows, inside) -> dict:
    """P5 on rows `(chrom, a, b, wrong)`: a pair is inside when both variants are in the BED."""
    n = {"inside": [0, 0], "outside": [0, 0]}
    for c, a, b, w in rows:
        side = "inside" if inside(c, a) and inside(c, b) else "outside"
        n[side][0] += 1
        n[side][1] += bool(w)
    rate = {s: k / max(1, m) for s, (m, k) in n.items()}
    wrong = n["inside"][1] + n["outside"][1]
    return {"pairs": n, "share_outside": n["outside"][1] / max(1, wrong),
            "ratio": rate["outside"] / rate["inside"] if rate["inside"] else float("inf")}


# ---- the family key with the test genome's own parents ----

def family_transmission(chrom: str, variants, parents) -> list:
    """Paternal allele of a heterozygous child of `parents` (father, mother) from the panel."""
    from pysam import bcftools

    import genome_common as common

    cache = KEYS / f"family_{parents[0]}_{parents[1]}_{chrom}.json"
    if not cache.exists():
        KEYS.mkdir(parents=True, exist_ok=True)
        panel = (common.KG_PANEL /
                 f"1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
        text = bcftools.view("-H", "-s", ",".join(parents), str(panel))
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


def family_key(platform: str, parents, chroms) -> dict:
    """Per chromosome `{"asm": pos -> assembly value, "tra": pos -> paternal allele}`."""
    import k1b_score as kb

    def one(c):
        table = kb.load_table((platform, c))
        tra = family_transmission(c, table["variants"], parents)
        return c, {"asm": {v[0]: x for v, x in zip(table["variants"], table["assembly"])
                           if x is not None},
                   "tra": {v[0]: t for v, t in zip(table["variants"], tra) if t is not None}}

    with ThreadPoolExecutor(8) as ex:
        return dict(ex.map(one, chroms))


# ---- P1 and P2 ----

def p12(item) -> list:
    depth, phaser, platform, key = item
    import confidence_score as cs
    import genome_common as common
    import l1_key_selection as l1
    from stratified_gate import strata
    from stratified_score import cells, grouped
    from w1b_score import POOL, pooled

    from phasebudget.confidence import learn_then_test_by_group

    chroms = list(common.CONFIRMATORY)
    stratum = strata(chroms)
    pool = {"HG002": grouped(l1.CAL, depth, phaser, chroms, stratum),
            **{s: grouped(f"hprc_{s}", depth, phaser, chroms, stratum) for s in POOL}}
    rows = l1.rows_with_pos(platform, depth, phaser, stratum)
    modes = {"all": rows}
    if key is not None:
        modes["key"] = l1.restrict(rows, key)
        modes["random"] = l1.matched_random(rows, modes["key"])
    out = []
    for beta in cs.BETAS:
        th = learn_then_test_by_group(pooled(pool), beta, cs.DELTA)
        for mode, sel in modes.items():
            for g, c in cells(th, [r[:4] for r in sel], beta).items():
                out.append((mode, depth, phaser, beta, g, c, bool(c[0]) and c[2] < P12_LEVEL))
    print(f"  P1/P2 {depth}x {phaser} done", flush=True)
    return out


# ---- P3 ----

def p3(label: str, platform: str) -> list:
    import confidence_score as cs
    import k1b_score as kb
    import stratified_score as ss
    from stratified_gate import strata

    from phasebudget.confidence import learn_then_test_by_group

    chroms = list(kb.CHROMS)
    stratum = strata(chroms)
    jobs = sorted({(p, c) for p in [*kb.POOL.values(), platform] for c in chroms})
    with ThreadPoolExecutor(12) as ex:
        tables = dict(zip(jobs, ex.map(kb.load_table, jobs)))
    out = []
    for phaser in kb.PHASERS:
        def rows_of(lab, p):
            return kb.pair_rows(lab, cs.load_rows(p, kb.DEPTH, phaser, chroms),
                                {c: tables[(p, c)] for c in chroms}, phaser, stratum)
        cal = [r for s, p in kb.POOL.items() for r in rows_of(s, p)]
        test = rows_of(label, platform)
        for beta in cs.BETAS:
            th = learn_then_test_by_group([r[:4] for r in cal if r[4]], beta, cs.DELTA)
            res = ss.cells(th, [r[:4] for r in test if r[4]], beta)
            for g in GROUPS:
                c = res.get(g, (0, 0, 1.0, 0, 0, 0, 0, False))
                share = c[4] / max(1, sum(r[0] == g for r in test))
                out.append((phaser, beta, g, th.get(g) is not None, c, share,
                            c[0] > 0 and c[2] < P3_LEVEL))
    return out


# ---- P4 ----

def p4(label: str, platform: str, workers: int) -> list:
    import l2_end_calls as l2

    chroms = list(l2.CHROMS)
    pool = [(s, p) for s, p in l2.POOL.items()]
    jobs = sorted({(p, c) for _, p in pool + [(label, platform)] for c in chroms})
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(l2.extract_job, jobs))
    work = [(lab, p, c) for lab, p in pool + [(label, platform)] for c in chroms]
    with ProcessPoolExecutor(workers) as ex:
        parts = list(ex.map(l2.tally, work))
    cal = l2.merge([p for p in parts if p["label"] != label])
    res = l2.merge([p for p in parts if p["label"] == label])
    level = 0.05 / (len(l2.PHASERS) * len(l2.TARGETS))
    out = []
    for ph in l2.PHASERS:
        for r in l2.TARGETS:
            groups = set()
            for g in l2.GROUPS:
                got = [v for (a, b, c, _), v in cal["picks"].items() if (a, b, c) == (ph, "assembly", g)]
                if l2.certified(len(got), sum(w for _, w in got), r):
                    groups.add(g)
            best: dict = {}
            for (a, b, g, (_, n)), pick in res["picks"].items():
                if (a, b) == (ph, "assembly") and g in groups and (n not in best or pick[0] < best[n][0]):
                    best[n] = pick
            n, k = len(best), sum(w for _, w in best.values())
            p = l2.p_exceed(k, n, r)
            long_share = sum(res["agree"].get((ph, "assembly", (P4_BIN, c)), 0)
                             for c in l2.CLASSES if (P4_BIN, c) in groups) / max(1, res["bins"][P4_BIN])
            out.append((ph, r, n, k, p, n > 0 and p < level, long_share, level))
    return out


# ---- P5 ----

def p5(platform: str) -> list:
    import confidence_score as cs
    import genome_common as common
    import q3_outside_benchmark as q3

    chroms = list(common.CONFIRMATORY)
    inside = q3.bed_lookup(GIAB_BED, chroms)
    return [(ph, outside_stats([(c, a, b, w) for c, a, _, w, b, _ in
                                cs.load_rows(platform, "30", ph, chroms)], inside))
            for ph in cs.PHASERS]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-assembly", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()

    import confidence_score as cs
    import genome_common as common

    started = time.time()
    label, platform, parents = TEST["dry" if args.dry_run else "real"]
    full = not args.no_assembly
    print(f"L7 ({'DRY RUN: ' + label + ' in place of NA12878, no NA12878 phase read' if args.dry_run else 'registered'}, "
          f"docs/design/l7-na12878-plan.md): test {label} ({platform}) against the Pedigree "
          f"truth; family key parents {parents} (1000 Genomes panel); depths {cs.DEPTHS}; "
          f"phasers {cs.PHASERS}; betas {cs.BETAS}; delta {cs.DELTA}; unit {cs.UNIT // 1000} kb; "
          f"P1 a cell in {P1_GROUPS} violates at {P12_LEVEL:.3e} on all links; P2 0 violating "
          f"key cells at {P12_LEVEL:.3e}; P3 K1b pool, 0 violating at {P3_LEVEL:.3e}, segdup at "
          f"{P3_BETA:g} trusted >= {P3_SHARE} for one phaser; P4 L2 pool, assembly arm, no "
          f"violating cell at 0.05/4, 10-100 kb share at r {P4_R:g} >= {P4_SHARE} for one "
          f"phaser; P5 30x, share of wrong pairs outside {GIAB_BED.name} >= {P5_SHARE} and "
          f"outside/inside rate >= {P5_RATIO} for both phasers; "
          f"{'all predictions' if full else 'P1 and P5 only (--no-assembly)'}", flush=True)
    verdicts = {}

    key = family_key(platform, parents, list(common.CONFIRMATORY)) if full else None
    if key is not None:
        print(f"  family key: assembly places {sum(len(k['asm']) for k in key.values())}, "
              f"transmission labels {sum(len(k['tra']) for k in key.values())} variants", flush=True)
    items = [(d, ph, platform, key) for d in cs.DEPTHS for ph in cs.PHASERS]
    with ProcessPoolExecutor(len(items)) as ex:
        cells_ = [r for part in ex.map(p12, items) for r in part]
    print("\n== P1, P2: mode depth phaser beta group | units, wrong units, p, pairs, trusted, "
          "wrong trusted")
    for mode, d, ph, beta, g, c, bad in cells_:
        print(f"  {mode:<6} {d}x {ph:<9} {beta:<6g} {g:<9} units {c[0]:>6} wrong {c[1]:>4} "
              f"p {c[2]:.2e} pairs {c[3]:>8} trusted {c[4]:>8} wrong trusted {c[6]:>5}"
              f"{'  VIOLATES' if bad else ''}")
    count = {m: sum(r[6] for r in cells_ if r[0] == m) for m in ("all", "key", "random")}
    print(f"  violating cells: {count}")
    verdicts["P1"] = any(r[6] for r in cells_ if r[0] == "all" and r[4] in P1_GROUPS)
    if full:
        verdicts["P2"] = count["key"] == 0

    if full:
        print("\n== P3: phaser beta group | certified, units, wrong units, p, trusted share")
        rows3 = p3(label, platform)
        for ph, beta, g, cert, c, share, bad in rows3:
            print(f"  {ph:<9} {beta:<6g} {g:<9} {'cert' if cert else 'none'} units {c[0]:>6} "
                  f"wrong {c[1]:>3} p {c[2]:.2e} share {share:.4f}{'  VIOLATES' if bad else ''}")
        verdicts["P3"] = (not any(r[6] for r in rows3) and any(
            r[3] and r[5] >= P3_SHARE for r in rows3 if r[1] == P3_BETA and r[2] == "segdup"))

        print("\n== P4: phaser r | genes picked, wrong picks, p, 10-100 kb share")
        rows4 = p4(label, platform, args.workers)
        for ph, r, n, k, p, bad, share, level in rows4:
            print(f"  {ph:<9} {r:<6g} picks {n:>5} wrong {k:>3} p {p:.2e} 10-100 kb share "
                  f"{share:.4f}{'  VIOLATES' if bad else ''} (level {level:.3e})")
        verdicts["P4"] = (not any(r[5] for r in rows4)
                          and any(r[6] >= P4_SHARE for r in rows4 if r[1] == P4_R))

    print("\n== P5: 30x, inside = both variants in the GIAB HG001 v4.2.1 BED")
    rows5 = p5(platform)
    for ph, s in rows5:
        (ni, ki), (no, ko) = s["pairs"]["inside"], s["pairs"]["outside"]
        print(f"  {ph:<9} inside {ki}/{ni} = {ki / max(1, ni):.2e}; outside {ko}/{no} = "
              f"{ko / max(1, no):.2e}; ratio {s['ratio']:.2f}; share of wrong outside "
              f"{s['share_outside']:.4f}")
    verdicts["P5"] = all(s["share_outside"] >= P5_SHARE and s["ratio"] >= P5_RATIO
                         for _, s in rows5)

    print()
    for name, ok in verdicts.items():
        print(f"{name}: {'supported' if ok else 'not supported'}{' (dry run, means nothing)' if args.dry_run else ''}")
    print(f"runtime {time.time() - started:.0f} s")


if __name__ == "__main__":
    main()
