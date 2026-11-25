"""Q4 development gate: certified cis or trans calls in clinical genes.

`docs/design/q4-clinical-cis-trans.md` states the call, its certification at
a per-call rate `r`, the population, the cells and the gate. In short: a pair
of truth-phased variants in one clinical gene, at most 100 kb apart, gets the
phaser's call (the parity of its chain of links) when every link is trusted at
its group's threshold for `b(m, r)`, the largest grid beta with `m * beta <=
r`, where `m` counts the chain's distinct (group, 100 kb unit) cells.

Steps:
1. per (genome, depth, phaser, chromosome), a subprocess with the genome's
   environment writes every input variant as `[pos1, phaser allele, truth
   allele, link score or None, indel]` to `data/r6_q4/chains/`;
2. per cell, thresholds from `learn_then_test_by_group` on the cached
   calibration rows (two folds by chromosome parity, or `cal_v5` for the
   transfer to `hg001_v5`), then every pair is enumerated and tallied.

The data were seen before: this is development, not a test.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/q4_clinical_cis_trans.py`
"""

from __future__ import annotations

import bisect
import gzip
import json
import os
import statistics
import subprocess
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Iterator, Optional, Sequence

BETAS = (1e-4, 3e-4, 1e-3, 3e-3, 1e-2)
TARGETS = (1e-3, 1e-2)
WINDOW = 100_000
DEPTHS = ("30", "20", "10")
PHASERS = ("whatshap", "longphase")
GENOMES = ("cal_v5", "cal_r10", "hprc_HG002", "hg001_v5")
TRANSFER = ("cal_v5", "hg001_v5")  # (calibration, test), the third test's setting
CELLS = len(GENOMES) * len(DEPTHS) * len(PHASERS) * len(TARGETS) \
    + len(DEPTHS) * len(PHASERS) * len(TARGETS)
LEVEL_A = 0.05 / CELLS
GATE_B = ("hprc_HG002", "30", 1e-2, 0.80)
GROUPS = ("indel", "segdup", "lowmap", "tandem", "difficult", "easy")
NAMED = ("PMS2", "SMN1", "SMN2", "GBA1", "CYP21A2", "HBA1", "HBA2", "STRC", "TNXB", "NEB",
         "CYP11B1", "CYP11B2", "NCF1", "PKD1", "CYP2D6", "SDHA", "CHEK2", "OTOA", "IKBKG",
         "HYDIN")
MIN_GENE_PAIRS = 10
DISTANCES = (1_000, 10_000, WINDOW)  # post hoc bins, added after the smoke run


# ---- the call and its certificate (pure; `tests/test_q4_clinical_cis_trans.py`) ----

def b_max(m: int, r: float, betas: Sequence[float] = BETAS) -> Optional[int]:
    """Index of the largest beta with `m * beta <= r`, or None."""
    best = None
    for k, beta in enumerate(betas):
        if m * beta <= r * (1 + 1e-9):
            best = k
    return best


def link_level(key: tuple, thresholds: Sequence[Optional[tuple]]) -> int:
    """Smallest beta index whose threshold trusts `key`; `len(thresholds)` if none does."""
    for k, t in enumerate(thresholds):
        if t is not None and key >= t:
            return k
    return len(thresholds)


def check_nested(thresholds: Sequence[Optional[tuple]]) -> None:
    """Raise unless a threshold, once it exists, stays and only loosens as beta grows."""
    for a, b in zip(thresholds, thresholds[1:]):
        if a is not None and (b is None or b > a):
            raise ValueError(f"thresholds not nested: {thresholds}")


def pairs(variants: Sequence[tuple], reach: Sequence[int], window: int) -> Iterator[tuple]:
    """Every pair of the population, with what certification needs.

    `variants[k]` is `(pos1, allele, truth, link)`: `allele` the phaser's first
    haplotype allele (None when unphased), `truth` the truth's (None when it
    does not phase the variant), `link` `(level, cell, group)` of the link to
    the next phased variant, None when there is none. `reach[k]` is the
    largest end of a listed gene holding variant `k` (0 if none).

    Yields `(i, j, status, cells, level, groups, covered)` for truth-phased
    `i < j` with `pos_j <= min(pos_i + window, reach[i])`. `status` is
    "unphased", "break" or "chain"; for a chain, `cells` is the live set of
    its links' cells (do not keep it; `m` is its size), `level` the largest
    link level, `groups` a live `{group: largest level}`, and `covered`
    whether the truth phases every phased variant strictly inside.
    """
    n = len(variants)
    for i in range(n):
        pos_i, a_i, h_i, _ = variants[i]
        if h_i is None:
            continue
        limit = min(pos_i + window, reach[i])
        cells: set = set()
        groups: dict = {}
        level, broken, covered, prev = -1, a_i is None, True, i
        for j in range(i + 1, n):
            pos_j, a_j, h_j, _ = variants[j]
            if pos_j > limit:
                break
            if a_j is not None and not broken:
                link = variants[prev][3]
                if link is None:
                    broken = True
                else:
                    lv, cell, group = link
                    cells.add(cell)
                    if lv > level:
                        level = lv
                    if lv > groups.get(group, -1):
                        groups[group] = lv
                    prev = j
            if h_j is not None:
                if a_i is None or a_j is None:
                    yield i, j, "unphased", None, 0, None, True
                elif broken:
                    yield i, j, "break", None, 0, None, True
                else:
                    yield i, j, "chain", cells, level, groups, covered
            if a_j is not None and h_j is None:
                covered = False


# ---- per-chromosome inputs ----

def extract(depth: str, phaser: str, chrom: str) -> list:
    """`[pos1, allele, truth, link score or None, indel]` per input variant, in the current
    platform's environment (`PHASEBUDGET_PLATFORM`)."""
    import genome_common as common
    from phase_confidence import phasing_from_vcf

    from phasebudget.confidence import pair_support, score

    extractor = common._module("extract_fragments.py")
    sc = common._module("score_hg002.py")
    variants = sc.variant_table(common.hetsnv_path(chrom))
    indel = {p for p, r, a in variants if len(r) != 1 or len(a) != 1}
    matrix, positions, _, _ = extractor.load_fragments(
        str(common.fragments_path(chrom, depth, "realign")))
    if list(positions) != sorted(positions):
        raise ValueError(f"{chrom}: positions not sorted")
    truth = common.truth_at(chrom, variants, positions)
    vcf = common.chrom_dir(chrom) / "external" / f"{phaser}_{depth}.vcf"
    vcf = vcf if vcf.exists() else vcf.with_name(vcf.name + ".gz")
    phasing = phasing_from_vcf(str(vcf), common.VCF_SAMPLE, chrom, positions)
    link = {s.first: list(score(s)) for s in pair_support(matrix, phasing)}
    return [[p + 1, phasing.alleles[j], truth[j], link.get(j), p + 1 in indel]
            for j, p in enumerate(positions)]


def chain_path(root: Path, genome: str, depth: str, phaser: str, chrom: str) -> Path:
    return root / "data" / "r6_q4" / "chains" / f"{genome}_{depth}_{phaser}_{chrom}.json.gz"


def _extract_job(job) -> str:
    root, genome, depth, phaser, chrom = job
    path = chain_path(root, genome, depth, phaser, chrom)
    if not path.exists():
        here = Path(__file__).resolve().parent
        env = dict(os.environ, PHASEBUDGET_PLATFORM=genome,
                   PYTHONPATH=f"{here}:{root / 'src'}:{root / 'benchmarks'}")
        env.pop("PHASEBUDGET_TRUTH", None)
        subprocess.run([sys.executable, __file__, "--extract", depth, phaser, chrom, str(path)],
                       env=env, check=True)
    return str(path)


def clinical_genes(chroms) -> dict:
    """`{chrom: [(start, end, name)]}` of the listed genes (`chain_relations` helpers)."""
    import chain_relations as cr

    return cr.genes(chroms, cr.clinical_symbols())


def gene_index(variants: Sequence[tuple], genes: Sequence[tuple]) -> tuple[list, list]:
    """`(reach, holding)`: per variant the largest end of a listed gene holding it, and
    `[(end, name)]` of those genes."""
    positions = [v[0] for v in variants]
    holding: list = [[] for _ in variants]
    for start, end, name in genes:
        for k in range(bisect.bisect_left(positions, start), bisect.bisect_right(positions, end)):
            holding[k].append((end, name))
    return [max((e for e, _ in h), default=0) for h in holding], holding


# ---- calibration ----

def calibration_rows(genome: str, depth: str, phaser: str, chroms, stratum) -> list:
    """`[(group, unit, score, wrong)]` from the cached rows, grouped as `stratified_score`."""
    import confidence_score as cs
    import genome_common as common
    from stratified_score import indels

    want = set(chroms)
    odd, even = list(common.ODD), list(common.EVEN)
    for sets in ([list(common.CONFIRMATORY)], [odd, even]):
        paths = [cs.CACHE / f"{genome}_{depth}_{phaser}_{'-'.join(s)}.json" for s in sets]
        if all(p.exists() for p in paths):
            rows = [r for s in sets for r in cs.load_rows(genome, depth, phaser, s)]
            break
    else:
        raise SystemExit(f"no cached rows for {genome} {depth}x {phaser}")
    ind = indels(genome, sorted(want))
    return [("indel" if (c, a) in ind or (c, b) in ind else stratum(c, a, b), (c, a // cs.UNIT),
             key, wrong)
            for c, a, key, wrong, b, _ in rows if c in want]


def thresholds_for(rows: list) -> dict:
    """`{group: [threshold per beta]}`, nested as the design note requires."""
    import confidence_score as cs

    from phasebudget.confidence import learn_then_test_by_group

    by_beta = [learn_then_test_by_group(rows, beta, cs.DELTA) for beta in BETAS]
    out = {g: [t.get(g) for t in by_beta] for g in {r[0] for r in rows}}
    for t in out.values():
        check_nested(t)
    return out


# ---- one cell: thresholds, enumeration, tallies ----

def plan_for(setting: str, genome: str) -> list:
    """`[(calibration genome, calibration chromosomes, test chromosomes)]` of a setting."""
    import genome_common as common

    odd, even, every = list(common.ODD), list(common.EVEN), list(common.CONFIRMATORY)
    if setting == "crossfit":
        return [(genome, even, odd), (genome, odd, even)]
    if setting == "transfer":
        return [(TRANSFER[0], every, every)]
    raise ValueError(setting)


def run_cell(job) -> dict:
    """Tallies of one (setting, genome, depth, phaser) for every target.

    `job` is `(setting, genome, depth, phaser, root, plan, truth_bound)`, `plan` as
    `plan_for` returns it. With `truth_bound` (`{group: per-unit truth error bound}`) a
    certified call counts toward the validity test only when the bounds of its cells sum
    to at most r/3 (`docs/design/q4-clinical-cis-trans.md`); without, every certified call does.
    """
    setting, genome, depth, phaser, root, spec, truth_bound = job
    import confidence_score as cs
    import genome_common as common
    from stratified_gate import strata

    chroms = list(common.CONFIRMATORY)
    stratum = strata(chroms)
    plan = [(thresholds_for(calibration_rows(cal, depth, phaser, cal_chroms, stratum)), test)
            for cal, cal_chroms, test in spec]
    genes = clinical_genes(chroms)
    nt = len(TARGETS)
    total = {"pairs": 0, "unphased": 0, "break": 0, "chain": 0, "chain_wrong": 0,
             "testable": [0] * nt, "testable_wrong": [0] * nt,
             "m_sum": [0] * nt, "certified": [0] * nt, "wrong": [0] * nt,
             "uncovered": [0] * nt, "uncovered_wrong": [0] * nt, "no_beta": [0] * nt,
             "untrusted": [dict.fromkeys(GROUPS, 0) for _ in TARGETS],
             "only_indel": [0] * nt, "dist_pairs": [0] * len(DISTANCES),
             "dist_certified": [[0] * len(DISTANCES) for _ in TARGETS]}
    per_gene: dict = {}  # name: [pairs, unphased, break, per target (certified, untrusted...)]
    chosen: dict = {}  # (name, target index): (crc, wrong)
    thresholds_printed = []
    for thresholds, test in plan:
        thresholds_printed.append({g: [None if t is None else [t[0], t[1]] for t in v]
                                   for g, v in thresholds.items()})
        for chrom in test:
            path = chain_path(root, genome, depth, phaser, chrom)
            with gzip.open(path, "rt") as handle:
                raw = json.load(handle)
            phased = [k for k, v in enumerate(raw) if v[1] is not None]
            nxt = dict(zip(phased, phased[1:]))
            variants = []
            for k, (pos, allele, truth, key, indel) in enumerate(raw):
                link = None
                if key is not None:
                    other = raw[nxt[k]]
                    group = ("indel" if indel or other[4] else stratum(chrom, pos, other[0]))
                    t = thresholds.get(group, [None] * len(BETAS))
                    link = (link_level(tuple(key), t), (group, pos // cs.UNIT), group)
                variants.append((pos, allele, truth, link))
            reach, holding = gene_index(variants, genes[chrom])
            for i, j, status, cells, level, groups, covered in pairs(variants, reach, WINDOW):
                pos_j = variants[j][0]
                names = [g for e, g in holding[i] if e >= pos_j]
                total["pairs"] += 1
                d = bisect.bisect_left(DISTANCES, pos_j - variants[i][0])
                total["dist_pairs"][d] += 1
                rows = []
                for name in names:
                    row = per_gene.get(name)
                    if row is None:
                        row = per_gene[name] = [0, 0, 0] + [0] * (nt * (2 + len(GROUPS)))
                    row[0] += 1
                    rows.append(row)
                if status != "chain":
                    total[status] += 1
                    for row in rows:
                        row[1 if status == "unphased" else 2] += 1
                    continue
                a_i, h_i = variants[i][1], variants[i][2]
                wrong = (a_i ^ h_i) != (variants[j][1] ^ variants[j][2])
                total["chain"] += 1
                total["chain_wrong"] += wrong
                m = len(cells)
                for x, r in enumerate(TARGETS):
                    base = 3 + x * (2 + len(GROUPS))
                    b = b_max(m, r)
                    if b is not None and level <= b:
                        total["certified"][x] += 1
                        total["dist_certified"][x][d] += 1
                        total["m_sum"][x] += m
                        total["wrong"][x] += wrong
                        if not covered:
                            total["uncovered"][x] += 1
                            total["uncovered_wrong"][x] += wrong
                        testable = truth_bound is None or sum(
                            truth_bound.get(c[0], 1.0) for c in cells) <= r / 3
                        total["testable"][x] += testable
                        total["testable_wrong"][x] += testable and wrong
                        for name, row in zip(names, rows):
                            row[base] += 1
                            row[base + 1] += wrong
                            if not testable:
                                continue
                            crc = zlib.crc32(f"{name}:{variants[i][0]}:{pos_j}".encode())
                            old = chosen.get((name, x))
                            if old is None or crc < old[0]:
                                chosen[(name, x)] = (crc, wrong)
                        continue
                    if b is None:
                        total["no_beta"][x] += 1
                        continue
                    bad = [g for g, lv in groups.items() if lv > b]
                    total["only_indel"][x] += bad == ["indel"]
                    for g in bad:
                        total["untrusted"][x][g] += 1
                        for row in rows:
                            row[base + 2 + GROUPS.index(g)] += 1
    validity = []
    for x in range(nt):
        picks = [w for (name, y), (_, w) in chosen.items() if y == x]
        validity.append((len(picks), sum(picks)))
    return {"setting": setting, "genome": genome, "depth": depth, "phaser": phaser,
            "total": total, "per_gene": per_gene, "validity": validity,
            "thresholds": thresholds_printed}


# ---- report ----

def p_exceed(k: int, n: int, r: float) -> float:
    from phasebudget.confidence import binomial_cdf

    return 1.0 if k == 0 else 1.0 - binomial_cdf(k - 1, n, r)


def gene_share(row: list, x: int) -> float:
    return row[3 + x * (2 + len(GROUPS))] / row[0]


def report(results: list, args) -> None:
    from cmrg_withholding import duplicated
    from scipy.stats import mannwhitneyu

    dup = duplicated()
    print("\ncell table (pairs: truth-phased, one listed gene, <= 100 kb; unique pairs)")
    print("setting   genome      depth phaser     r      pairs      cert_share  certified  wrong"
          "  rate      mean_m  uncov_cert  uncov_wrong  unphased  break     no_beta  "
          "genes_n  genes_k  p_exceed  genes_all_cert  median_gene_share")
    violations, gate_b = [], {}
    for res in results:
        t = res["total"]
        for x, r in enumerate(TARGETS):
            n, k = res["validity"][x]
            p = p_exceed(k, n, r)
            if p < LEVEL_A:
                violations.append((res, r, n, k, p))
            genes = [row for row in res["per_gene"].values()]
            full = sum(gene_share(row, x) == 1.0 for row in genes)
            big = [gene_share(row, x) for row in genes if row[0] >= MIN_GENE_PAIRS]
            cert = t["certified"][x]
            share = cert / max(1, t["pairs"])
            if (res["setting"], res["genome"], res["depth"], r) == ("crossfit",) + GATE_B[:3]:
                gate_b[res["phaser"]] = share
            print(f"{res['setting']:<9} {res['genome']:<11} {res['depth']:>5} {res['phaser']:<10}"
                  f" {r:<6g} {t['pairs']:>9} {share:>11.4f} {cert:>10} {t['wrong'][x]:>6}"
                  f"  {t['wrong'][x] / max(1, cert):<8.2e} {t['m_sum'][x] / max(1, cert):>7.2f}"
                  f" {t['uncovered'][x]:>11} {t['uncovered_wrong'][x]:>12} {t['unphased']:>9}"
                  f" {t['break']:>8} {t['no_beta'][x]:>8} {n:>8} {k:>8}  {p:<8.2e}"
                  f" {full:>7}/{len(genes):<7} {statistics.median(big) if big else float('nan'):.4f}")
    print("\nabstention by untrusted group (pairs with an untrusted link in the group; a pair can "
          "count in several)")
    for res in results:
        for x, r in enumerate(TARGETS):
            u = res["total"]["untrusted"][x]
            print(f"  {res['setting']:<9} {res['genome']:<11} {res['depth']:>3}x {res['phaser']:<10}"
                  f" r {r:<6g} " + ", ".join(f"{g} {u[g]}" for g in GROUPS))
    print("\npost hoc, added after the smoke run on hprc_HG002 30x showed the indel group "
          "blocking most chains: pairs whose only untrusted links are indel links, and the "
          f"certified share by distance bin (upper ends {DISTANCES} bp)")
    for res in results:
        t = res["total"]
        for x, r in enumerate(TARGETS):
            print(f"  {res['setting']:<9} {res['genome']:<11} {res['depth']:>3}x {res['phaser']:<10}"
                  f" r {r:<6g} only indel {t['only_indel'][x] / max(1, t['pairs']):.4f}; by "
                  "distance " + ", ".join(
                      f"<= {e}: {c / max(1, n):.4f} of {n}" for e, n, c in
                      zip(DISTANCES, t["dist_pairs"], t["dist_certified"][x])))
    print("\nthresholds per fold (score, reads) per beta " + str(BETAS) + ", 30x only")
    for res in results:
        if res["depth"] != "30":
            continue
        for f, th in enumerate(res["thresholds"]):
            print(f"  {res['setting']} {res['genome']} {res['phaser']} fold {f}: "
                  + "; ".join(f"{g} " + " ".join("none" if t is None else f"{t[0]:.4f}/{t[1]}"
                                                 for t in th[g]) for g in sorted(th)))
    x = TARGETS.index(1e-2)
    print(f"\nabstaining genes at 30x, r = 1e-2: genes with >= {MIN_GENE_PAIRS} pairs and certified"
          " share below 0.5, lowest first (share, pairs; causes as shares of the gene's pairs: "
          "unphased u, break b, untrusted by group); D marks the GENCODE duplicated list")
    for res in results:
        if res["depth"] != "30":
            continue
        base = 3 + x * (2 + len(GROUPS))
        low = sorted((gene_share(row, x), name, row) for name, row in res["per_gene"].items()
                     if row[0] >= MIN_GENE_PAIRS and gene_share(row, x) < 0.5)
        print(f"  {res['setting']} {res['genome']} {res['phaser']}: {len(low)} genes")
        for share, name, row in low[:args.list]:
            causes = [f"u {row[1] / row[0]:.2f}", f"b {row[2] / row[0]:.2f}"] + [
                f"{g} {row[base + 2 + k] / row[0]:.2f}" for k, g in enumerate(GROUPS)
                if row[base + 2 + k]]
            print(f"    {name}{' D' if name in dup else ''} {share:.3f} ({row[0]}): "
                  + ", ".join(causes))
    print("\nnamed paralog genes at 30x, r = 1e-2: certified share (pairs); '-' no pairs")
    for res in results:
        if res["depth"] != "30":
            continue
        print(f"  {res['setting']} {res['genome']} {res['phaser']}: " + ", ".join(
            f"{g} " + (f"{gene_share(res['per_gene'][g], x):.2f} ({res['per_gene'][g][0]})"
                       if g in res["per_gene"] else "-") for g in NAMED))
    print(f"\nduplicated genes (GENCODE pseudogene list) against the rest, 30x, r = 1e-2, genes "
          f"with >= {MIN_GENE_PAIRS} pairs; one-sided Mann-Whitney, descriptive")
    for res in results:
        if res["depth"] != "30":
            continue
        shares = {g: gene_share(row, x) for g, row in res["per_gene"].items()
                  if row[0] >= MIN_GENE_PAIRS}
        a = [v for g, v in shares.items() if g in dup]
        b = [v for g, v in shares.items() if g not in dup]
        p = mannwhitneyu(a, b, alternative="less").pvalue if a and b else float("nan")
        print(f"  {res['setting']} {res['genome']} {res['phaser']}: duplicated {len(a)} genes "
              f"median {statistics.median(a) if a else float('nan'):.3f}, other {len(b)} median "
              f"{statistics.median(b) if b else float('nan'):.3f}, p {p:.2e}")
    print(f"\n(a) validity: {sum(len(r['validity']) for r in results)} cells of {CELLS} "
          f"registered, level {LEVEL_A:.2e}, violating {len(violations)}"
          + "".join(f"\n  {v[0]['setting']} {v[0]['genome']} {v[0]['depth']}x {v[0]['phaser']} "
                    f"r {v[1]:g}: {v[3]} wrong of {v[2]} genes, p {v[4]:.2e}" for v in violations))
    ok_b = any(s >= GATE_B[3] for s in gate_b.values())
    print(f"(b) usefulness: {GATE_B[0]} {GATE_B[1]}x r {GATE_B[2]:g} certified share "
          + ", ".join(f"{ph} {s:.4f}" for ph, s in gate_b.items())
          + f"; need >= {GATE_B[3]} for one phaser: {'holds' if ok_b else 'fails'}")
    passed = not violations and ok_b and len(results) == len(GENOMES) * 6 + 6
    print(f"Q4 gate: {'PASSES' if passed else 'DOES NOT PASS'}"
          + ("" if len(results) == len(GENOMES) * 6 + 6 else " (not every cell was run)"))


def main() -> None:
    import argparse

    import genome_common as common

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genomes", nargs="+", default=list(GENOMES))
    parser.add_argument("--depths", nargs="+", default=list(DEPTHS))
    parser.add_argument("--no-transfer", action="store_true")
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--list", type=int, default=40, help="abstaining genes printed per cell")
    args = parser.parse_args()
    root = common.ROOT
    chroms = list(common.CONFIRMATORY)
    start = time.time()
    print("Q4 development gate, docs/design/q4-clinical-cis-trans.md (data seen; not registered)")
    print(f"genomes {' '.join(args.genomes)}; depths {' '.join(args.depths)}; phasers "
          f"{' '.join(PHASERS)}; chromosomes {' '.join(chroms)}; calibration: two folds by "
          f"chromosome parity within each genome (crossfit), and {TRANSFER[1]} with thresholds "
          f"from {TRANSFER[0]} on all chromosomes (transfer)"
          + (" [transfer skipped]" if args.no_transfer else ""))
    print(f"beta grid {BETAS}; targets r {TARGETS}; certified when every link is trusted at "
          f"b(m, r), the largest beta with m * beta <= r, m = distinct (group, 100 kb unit) cells "
          f"of the chain's links; groups {', '.join(GROUPS)}; delta 0.05; Learn then Test grid "
          f"from 0.5 in steps of 0.0025")
    print(f"pairs: truth-phased input variants in one protein-coding gene of GENCODE v50 basic "
          f"whose ClinVar symbol has >= 10 P/LP alleles (chain_relations.clinical_symbols), at "
          f"most {WINDOW} bp apart; gene stats over genes with >= {MIN_GENE_PAIRS} pairs")
    print(f"gate (a): per cell one certified call per gene (smallest CRC32 of gene:pos_i:pos_j), "
          f"exact binomial P(X >= k | n, r) < 0.05/{CELLS} = {LEVEL_A:.2e} violates; gate (b): "
          f"{GATE_B[0]} {GATE_B[1]}x r {GATE_B[2]:g} crossfit certified share >= {GATE_B[3]} for "
          f"one phaser")
    jobs = []
    for genome in args.genomes:
        for depth in args.depths:
            for phaser in PHASERS:
                jobs += [(root, genome, depth, phaser, c) for c in chroms]
    with ProcessPoolExecutor(args.workers) as pool:
        for _ in pool.map(_extract_job, jobs):
            pass
    print(f"extraction done: {len(jobs)} chromosome files ({time.time() - start:.0f} s)",
          flush=True)
    cells = [("crossfit", g, d, ph, root, plan_for("crossfit", g), None)
             for g in args.genomes for d in args.depths for ph in PHASERS]
    if not args.no_transfer and TRANSFER[1] in args.genomes:
        cells += [("transfer", TRANSFER[1], d, ph, root, plan_for("transfer", TRANSFER[1]), None)
                  for d in args.depths for ph in PHASERS]
    with ProcessPoolExecutor(args.workers) as pool:
        results = list(pool.map(run_cell, cells))
    report(results, args)
    print(f"\nwall time {time.time() - start:.0f} s")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--extract":
        depth, phaser, chrom, out = sys.argv[2:]
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        data = extract(depth, phaser, chrom)
        tmp = out + f".{os.getpid()}.part"
        with gzip.open(tmp, "wt") as handle:
            json.dump(data, handle)
        os.replace(tmp, out)
    else:
        main()
