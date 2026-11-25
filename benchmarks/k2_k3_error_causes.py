"""Round 7, K2 and K3 (development): why phase relations are wrong, and HG001's benchmark gap.

`docs/design/k2-k3-error-causes.md` states the cause split, both gates and the platforms. Rows
per (platform, phaser, chromosome) come from `q3_outside_benchmark.build_rows`, run in a
subprocess with that platform's environment, and are cached in `data/k2/rows/`.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/k2_k3_error_causes.py`
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

import q3_outside_benchmark as q3  # noqa: E402

ROOT = HERE.parent
DEPTH = "30"
CHROMS = q3.CHROMS
PHASERS = q3.PHASERS
GROUPS = ("all",) + q3.GROUPS
GATE_PLATFORMS = ("w4_grch38", "cal_v5")
GATE_GROUPS = ("segdup", "lowmap", "tandem", "indel")
GATE_DATA, MIN_WRONG = 0.6, 20
K3_PLATFORM, K3_SHARE = "hg001_v5", 0.3
K3_BED = ROOT / "data" / "giab_hg001" / "HG001_GRCh38_1_22_v4.2.1_benchmark.bed"
PEDIGREE = ("NA12879", "NA12881", "NA12882", "NA12885", "NA12886")
# W1 genomes with no 1000 Genomes parents: no consensus truth, not scored in W1.
W1_UNSCORED = ("HG01123", "HG02109", "HG02486", "HG02559", "NA21309")


def platforms() -> list[str]:
    w1 = sorted(p.name[4:] for p in (ROOT / "data").glob("ext_hprc_*")
                if p.name[9:] != "HG002" and p.name[9:] not in W1_UNSCORED)
    return (list(GATE_PLATFORMS) + ["cal_r10", "hprc_HG002", K3_PLATFORM]
            + [f"ped_{s}" for s in PEDIGREE] + [f"phifi_{s}" for s in PEDIGREE] + w1)


def rows_path(platform: str, phaser: str, chrom: str) -> Path:
    return ROOT / "data" / "k2" / "rows" / f"{platform}_{phaser}_{DEPTH}_{chrom}.json"


def rows_job(job: tuple) -> Path:
    platform, phaser, chrom = job
    path = rows_path(platform, phaser, chrom)
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        env = dict(os.environ, PHASEBUDGET_PLATFORM=platform,
                   PYTHONPATH=f"{ROOT / 'src'}:{HERE}")
        env.pop("PHASEBUDGET_TRUTH", None)
        part = path.with_name(path.name + ".part")
        subprocess.run([sys.executable, __file__, "--rows", platform, phaser, chrom, str(part)],
                       env=env, check=True)
        part.replace(path)
    return path


def causes(wrongs) -> dict:
    """`{cause: count}` over wrong rows `[.., agree, disagree, .., costs]` (Q3's split)."""
    out = dict.fromkeys(q3.CAUSES, 0)
    for r in wrongs:
        out[q3.cause(q3.majority(r[2], r[3]), q3.optimum((r[7][0], r[7][1]), r[7][2]))] += 1
    return out


def line(label: str, counts: dict) -> str:
    n = sum(counts.values())
    return (f"  {label} wrong {n:>6}: " + ", ".join(
        f"{k} {v} ({v / max(1, n):.3f})" for k, v in counts.items()))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", nargs=4, metavar=("PLATFORM", "PHASER", "CHROM", "OUT"))
    parser.add_argument("--workers", type=int, default=24)
    parser.add_argument("--platforms", nargs="+")
    parser.add_argument("--k3-platform", default=K3_PLATFORM,
                        help="K3b (docs/design/k3b-hg005-plan.md) passes k3_hg005")
    parser.add_argument("--k3-bed", type=Path, default=K3_BED)
    args = parser.parse_args()
    if args.rows:
        platform, phaser, chrom, out = args.rows
        Path(out).write_text(json.dumps(q3.build_rows(chrom, phaser, platform)))
        return

    from stratified_gate import strata

    k3_platform, k3_bed = args.k3_platform, args.k3_bed
    started = time.time()
    names = args.platforms or platforms()
    print(f"K2 and K3 (development, docs/design/k2-k3-error-causes.md): {DEPTH}x runs of "
          f"{len(names)} platforms ({' '.join(names)}), whatshap and LongPhase as stored; "
          f"consecutive phased pairs both of whose variants the platform's primary truth phases; "
          f"cause window {q3.FLANK} variants per side, WhatsHap selection port at cap {q3.CAP}, "
          f"exact CoverageDP boundary costs; groups indel first, then GIAB v3.6 stratum; K2 gate: "
          f"data share >= {GATE_DATA} for both phasers in each of {', '.join(GATE_GROUPS)} with >= "
          f"{MIN_WRONG} wrong pairs, on {', '.join(GATE_PLATFORMS)}; K3 gate: {k3_platform}, share "
          f"of wrong pairs outside {k3_bed.name} >= {K3_SHARE} for both phasers; chromosomes "
          f"{' '.join(CHROMS)}", flush=True)
    jobs = [(p, ph, c) for p in names for ph in PHASERS for c in CHROMS]
    with ThreadPoolExecutor(args.workers) as pool:
        list(pool.map(rows_job, jobs))
    print(f"rows ready in {time.time() - started:.0f} s", flush=True)

    stratum = strata(list(CHROMS))
    inside = q3.bed_lookup(k3_bed, CHROMS)
    k2_ok, k3_ok = True, True
    print("\n1. Causes of wrong pairs per platform, phaser and group: none = no read calls both "
          "variants; data = spanning majority backs the wrong relation or ties, window optimum "
          "not right; search = window optimum strictly right; other = majority right, optimum "
          "wrong or tied")
    pooled: dict = {}
    for platform in names:
        for phaser in PHASERS:
            rows = []
            for chrom in CHROMS:
                for r in json.loads(rows_path(platform, phaser, chrom).read_text()):
                    group = "indel" if r[6] else stratum(chrom, r[0], r[1])
                    rows.append((chrom, group, r))
            for group in GROUPS:
                wrongs = [r for _, g, r in rows if r[5] and group in ("all", g)]
                counts = causes(wrongs)
                pairs = sum(group in ("all", g) for _, g, _ in rows)
                print(line(f"{platform:<16} {phaser:<9} {group:<9} pairs {pairs:>8}", counts))
                if platform.startswith("hprc_") and platform != "hprc_HG002":
                    acc = pooled.setdefault((phaser, group), dict.fromkeys(q3.CAUSES, 0))
                    for k, v in counts.items():
                        acc[k] += v
                n = sum(counts.values())
                if platform in GATE_PLATFORMS and group in GATE_GROUPS and n >= MIN_WRONG:
                    share = counts["data"] / n
                    if share < GATE_DATA:
                        k2_ok = False
                        print(f"    K2 gate fails here: data share {share:.3f} < {GATE_DATA}")
            if platform == k3_platform:
                print(f"\n2. K3: {platform} {phaser}, inside = both variants in {k3_bed.name}")
                for group in GROUPS:
                    sel = [(c, r) for c, g, r in rows if group in ("all", g)]
                    cut = {"inside": [(c, r) for c, r in sel if inside(c, r[0]) and inside(c, r[1])]}
                    cut["outside"] = [(c, r) for c, r in sel
                                      if not (inside(c, r[0]) and inside(c, r[1]))]
                    stats = {k: q3.unit_shares([(c, r[0], r[5]) for c, r in v])
                             for k, v in cut.items()}
                    ratio, share, _ = q3.gate(stats["inside"], stats["outside"])
                    for k, (n, w, u, wu) in stats.items():
                        print(f"  {group:<9} {k:<7} pairs {n:>8} wrong {w:>6} rate "
                              f"{w / max(1, n):.3e} units {u:>6} wrong {wu:>5}")
                        print("  " + line(f"{group:<9} {k:<7}",
                                          causes([r for _, r in cut[k] if r[5]])))
                    print(f"  {group:<9} outside/inside per pair {ratio:.2f}; share of wrong "
                          f"outside {share:.3f}")
                    if group == "all":
                        print(f"  K3 {phaser}: share {share:.3f} (>= {K3_SHARE}): "
                              f"{'holds' if share >= K3_SHARE else 'fails'}\n")
                        k3_ok &= share >= K3_SHARE
    if pooled:
        print("\n3. W1 genomes pooled (consensus truth)")
        for (phaser, group), counts in pooled.items():
            print(line(f"{phaser:<9} {group:<9}", counts))
    print(f"\nK2 gate: {'passed' if k2_ok else 'not passed'}")
    if k3_platform in names:
        print(f"K3 gate: {'finds the gap' if k3_ok else 'does not find the gap'}")
    print(f"runtime {time.time() - started:.0f} s")


if __name__ == "__main__":
    main()
