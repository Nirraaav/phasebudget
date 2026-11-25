"""K1b (registered, `docs/design/k1b-agreement-plan.md`): assembly-agreement certificates.

Per (platform, chromosome), a subprocess with that platform's environment writes the input
variants, both phasers' 30x alleles and blocks, and the assembly's values (Q6's contig rule on
the genome's dipcall outputs) to `data/k1/tables/`. Scores and wrong flags are the cached rows
of `confidence_score.load_rows`. Pairs that do not agree are never trusted.

`--dry-run` reads no test genome: HG002's odd chromosomes, split into two pseudo-genomes, and
every Pedigree pool genome whose dipcall output exists form the pool, and HG002's even
chromosomes stand in for a test genome.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/k1b_score.py [--dry-run]`
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

import k1_assembly_agreement as k1  # noqa: E402
import q6_assembly_phase as q6  # noqa: E402

ROOT = HERE.parent
DEPTH = "30"
CHROMS = tuple(f"chr{i}" for i in range(1, 23) if i != 20)
PHASERS = ("whatshap", "longphase")
GROUPS = ("indel", "segdup", "lowmap", "tandem", "difficult", "easy")
POOL = {"HG002": "hprc_HG002", **{s: f"phifi_{s}" for s in
                                  ("NA12879", "NA12881", "NA12882", "NA12885", "NA12886")}}
TESTS = {"NA12877": "phifi_NA12877", "HG005": "k1_hg005"}
CELLS = len(TESTS) * len(PHASERS) * 3 * len(GROUPS)  # 72, three betas
LEVEL = 0.05 / CELLS
FLOOR_BETA, MIN_TRUST = 1e-3, 0.5
USEFUL = {"NA12877": ("segdup", "tandem", "indel"), "HG005": ("segdup", "tandem")}


def dip_of(sample: str) -> tuple[Path, str]:
    if sample == "HG002":
        return ROOT / "data" / "r6_q6" / "dip", "HG002"  # Q6's assembly
    return ROOT / "data" / "k1" / sample / "dip", sample


def read_table(chrom: str) -> dict:
    """One chromosome's variants, phasers and assembly values in the current platform's
    environment."""
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
    table["assembly"] = q6.dipcall_values(chrom, variants, *dip_of(common.SAMPLE))
    return table


def table_path(platform: str, chrom: str) -> Path:
    return ROOT / "data" / "k1" / "tables" / f"{platform}_{DEPTH}_{chrom}.json"


def load_table(job: tuple) -> dict:
    platform, chrom = job
    path = table_path(platform, chrom)
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        env = dict(os.environ, PHASEBUDGET_PLATFORM=platform,
                   PYTHONPATH=f"{ROOT / 'src'}:{HERE}")
        part = path.with_name(path.name + ".part")
        subprocess.run([sys.executable, __file__, "--table", chrom, str(part)], env=env,
                       check=True)
        part.replace(path)
    return json.loads(path.read_text())


def pair_rows(label: str, rows, tables: dict, phaser: str, stratum) -> list:
    """`(group, unit, score, wrong, agree, assembly relates, both wrong)` per cached row
    `(chrom, a, score, wrong, b, _)` of the chromosomes in `tables`; the unit is
    `(label, chrom, a // 100 kb)`."""
    import confidence_score as cs

    index = {c: {} for c in tables}
    for c, t in tables.items():
        for k, v in enumerate(t["variants"]):
            index[c].setdefault(v[0], k)
    out = []
    for c, a, key, wrong, b, _ in rows:
        if c not in tables:
            continue
        t, ia, ib = tables[c], index[c][a], index[c][b]
        va, vb = t["variants"][ia], t["variants"][ib]
        indel = any(len(x) != 1 for x in (va[1], va[2], vb[1], vb[2]))
        group = "indel" if indel else stratum(c, a, b)
        mine = q6.relation(t[phaser][ia], t[phaser][ib])
        asm = q6.relation(t["assembly"][ia], t["assembly"][ib])
        agree = k1.agrees(mine, asm)
        out.append((group, (label, c, a // cs.UNIT), key, bool(wrong), agree, asm is not None,
                    bool(wrong) and agree))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", nargs=2, metavar=("CHROM", "OUT"))
    parser.add_argument("--dry-run", action="store_true")
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
    odd = [c for c in CHROMS if int(c[3:]) % 2]
    even = [c for c in CHROMS if not int(c[3:]) % 2]
    if args.dry_run:
        # (label, platform, chromosomes)
        pool = [("HG002a", "hprc_HG002", odd[0::2]), ("HG002b", "hprc_HG002", odd[1::2])]
        pool += [(s, p, list(CHROMS)) for s, p in POOL.items()
                 if s != "HG002" and (dip_of(s)[0] / f"{s}.dip.vcf.gz").exists()]
        tests = [("HG002even", "hprc_HG002", even)]
    else:
        pool = [(s, p, list(CHROMS)) for s, p in POOL.items()]
        tests = [(s, p, list(CHROMS)) for s, p in TESTS.items()]
    print(f"K1b ({'DRY RUN, no test genome read' if args.dry_run else 'registered'}, "
          f"docs/design/k1b-agreement-plan.md): {DEPTH}x, whatshap and LongPhase; agreement = "
          f"the single-sample hifiasm 0.25.0 assembly relates the pair (Q6 contig rule on dipcall "
          f"0.3) with the phaser's relation; pool {', '.join(f'{s} ({p}, {len(c)} chromosomes)' for s, p, c in pool)}; "
          f"tests {', '.join(f'{s} ({p}, {len(c)} chromosomes)' for s, p, c in tests)}; "
          f"learn_then_test_by_group on pooled agreeing pairs, unit (genome, chromosome, "
          f"{cs.UNIT // 1000} kb), betas {cs.BETAS}, delta {cs.DELTA}; {CELLS} registered cells, "
          f"level {LEVEL:.3e}; K1b-2 segdup certifies at {FLOOR_BETA:g} for both phasers; K1b-3 "
          f"trusted share >= {MIN_TRUST} at {FLOOR_BETA:g} for one phaser in {USEFUL}", flush=True)

    stratum = strata(list(CHROMS))
    jobs = sorted({(p, c) for _, p, cs_ in pool + tests for c in cs_})
    with ThreadPoolExecutor(args.workers) as pool_exec:
        loaded = dict(zip(jobs, pool_exec.map(load_table, jobs)))

    def rows_of(label, platform, chroms, phaser):
        rows = cs.load_rows(platform, DEPTH, phaser, list(CHROMS))
        tables = {c: loaded[(platform, c)] for c in chroms}
        return pair_rows(label, rows, tables, phaser, stratum)

    verdict_1, verdict_2, verdict_3 = True, True, True
    useful: dict = {}
    for phaser in PHASERS:
        cal = [r for label, p, c in pool for r in rows_of(label, p, c, phaser)]
        print(f"\n== {phaser}: pool pairs {len(cal)}, agreeing {sum(r[4] for r in cal)}")
        thresholds = {}
        for beta in cs.BETAS:
            for mode in ("agree", "all"):
                sel = [r[:4] for r in cal if mode == "all" or r[4]]
                thresholds[(beta, mode)] = learn_then_test_by_group(sel, beta, cs.DELTA)
                units = {}
                for g, u, _, _ in sel:
                    units.setdefault(g, set()).add(u)
                print(f"  pool {mode:<5} beta {beta:<6g} " + ", ".join(
                    f"{g} {'none' if thresholds[(beta, mode)].get(g) is None else 'cert'}"
                    f" ({len(units.get(g, ()))} units)" for g in GROUPS))
        seg = thresholds[(FLOOR_BETA, "agree")].get("segdup")
        verdict_2 &= seg is not None
        print(f"  K1b-2 {phaser}: segdup at {FLOOR_BETA:g} {'certifies' if seg else 'abstains'}")

        for label, platform, chroms in tests:
            test = rows_of(label, platform, chroms, phaser)
            print(f"\n  test {label} {phaser}: pairs {len(test)}")
            print("  K1b-5 per group: pairs, assembly relates, agree, wrong among agreeing, "
                  "wrong among others, both wrong")
            for g in GROUPS:
                mine = [r for r in test if r[0] == g]
                n = max(1, len(mine))
                agree = [r for r in mine if r[4]]
                print(f"    {g:<9} pairs {len(mine):>8} relates {sum(r[5] for r in mine) / n:.4f} "
                      f"agree {len(agree) / n:.4f} wrong agreeing {sum(r[3] for r in agree):>5} "
                      f"others {sum(r[3] for r in mine if not r[4]):>5} both {sum(r[6] for r in mine):>4}")
            print("  cells: threshold, units with a trusted pair, wrong units, "
                  "P(X >= k | n, beta), trusted pairs, trusted share of the group's pairs, wrong "
                  "trusted")
            for beta in cs.BETAS:
                for mode in ("agree", "all"):
                    t = thresholds[(beta, mode)]
                    res = ss.cells(t, [r[:4] for r in test if mode == "all" or r[4]], beta)
                    for g in GROUPS:
                        n_all = sum(r[0] == g for r in test)
                        n, k, p, _, kept, _, wrong_kept, _ = res.get(
                            g, (0, 0, 1.0, 0, 0, 0, 0, False))
                        share = kept / max(1, n_all)
                        bad = mode == "agree" and n > 0 and p < LEVEL
                        if bad:
                            verdict_1 = False
                        if mode == "agree" and beta == FLOOR_BETA:
                            useful.setdefault((label, g), []).append(
                                (phaser, t.get(g) is not None, share))
                        th = t.get(g)
                        print(f"    {beta:<6g} {mode:<5} {g:<9} "
                              f"{'none' if th is None else f'({th[0]:.4f}, {th[1]})':<16} units {n:>6} "
                              f"wrong {k:>3} p {p:.2e} trusted {kept:>8} share {share:.4f} wrong "
                              f"trusted {wrong_kept:>4}{'  VIOLATES' if bad else ''}")

    print()
    if not args.dry_run:
        for label, groups in USEFUL.items():
            for g in groups:
                got = useful.get((label, g), [])
                ok = any(cert and share >= MIN_TRUST for _, cert, share in got)
                verdict_3 &= ok
                print(f"K1b-3 {label} {g}: " + ", ".join(
                    f"{ph} {'certifies' if c else 'none'} share {s:.4f}" for ph, c, s in got)
                    + f": {'holds' if ok else 'fails'}")
        print(f"K1b-1 validity: {'supported' if verdict_1 else 'not supported'}")
        print(f"K1b-2 floor: {'supported' if verdict_2 else 'not supported'}")
        print(f"K1b-3 usefulness: {'supported' if verdict_3 else 'not supported'}")
    else:
        print(f"dry run: violating agreement cells at {LEVEL:.3e}: "
              f"{'none' if verdict_1 else 'some'}; segdup at {FLOOR_BETA:g} "
              f"{'certifies' if verdict_2 else 'abstains'} for both phasers")
    print(f"runtime {time.time() - started:.0f} s")


if __name__ == "__main__":
    main()
