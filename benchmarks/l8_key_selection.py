"""L8 (`docs/design/l8-hprc-key-selection.md`): L1's key-selection prediction on HPRC genomes.

Each test genome's HiFi rows are scored against its own dipcall truth alone (HPRC release 2
assembly, `PHASEBUDGET_TRUTH=dipcall`). A dipcall-labelled link is `key` when 1000 Genomes
transmission labels both variants with dipcall's paternal allele (the W1 consensus key),
`disagree` when transmission labels a variant with the other allele (key noise, left out), and
`missing` otherwise. W1b's pooled thresholds (HG002 v5.0q and the 20 consensus-truth pool
genomes) are applied to all non-disagreeing links, to key links, and to a random subset of the
same share per group.

`--dev` tests the 20 pool genomes themselves (thresholds in sample, development only);
`--test` tests the 19 held-out genomes, and may run only after the design note is committed.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/l8_key_selection.py --dev|--test`
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

ROOT = HERE.parent
KEYS = ROOT / "data" / "l8" / "keys"
GROUPS = ("easy", "segdup", "lowmap", "tandem", "difficult", "indel")
HARD = ("segdup", "lowmap", "tandem")
CELLS_PER_GENOME = 3 * 2 * 3 * len(GROUPS)  # depths, phasers, betas, groups


def link_class(la, lb) -> str:
    """`key`, `disagree` or `missing` for a link whose two variants carry `(dipcall,
    transmission)` labels; dipcall labels both, since only such links have rows."""
    if any(t is not None and t != d for d, t in (la, lb)):
        return "disagree"
    if la[1] is not None and lb[1] is not None:
        return "key"
    return "missing"


def key_labels(chrom: str) -> list:
    """`[pos1, dipcall, transmission]` per input variant, in the current platform's environment."""
    import genome_common as common

    sc = common._module("score_hg002.py")
    variants = sc.variant_table(common.hetsnv_path(chrom))
    d = common.dipcall_truth(chrom, variants)
    t = common.kg_transmission(chrom, variants)
    return [[v[0], a, b] for v, a, b in zip(variants, d, t)]


def load_labels(job) -> dict:
    platform, chrom = job
    path = KEYS / f"{platform}_{chrom}.json"
    if not path.exists():
        KEYS.mkdir(parents=True, exist_ok=True)
        part = path.with_name(f"{path.name}.{os.getpid()}.part")
        env = dict(os.environ, PHASEBUDGET_PLATFORM=platform, PYTHONPATH=f"{ROOT / 'src'}:{HERE}")
        subprocess.run([sys.executable, __file__, "--labels", chrom, str(part)], env=env,
                       check=True)
        part.replace(path)
    return {p: (d, t) for p, d, t in json.loads(path.read_text())}


def dipcall_rows(platform, depth, phaser, stratum) -> list:
    """L1's `(group, unit, score, wrong, chrom, a, b)` rows, wrong against dipcall alone."""
    from l1_key_selection import rows_with_pos

    os.environ["PHASEBUDGET_TRUTH"] = "dipcall"
    try:
        return rows_with_pos(platform, depth, phaser, stratum)
    finally:
        del os.environ["PHASEBUDGET_TRUTH"]


def rate(k: int, n: int) -> str:
    return f"{k}/{n} = {k / max(1, n):.2e}"


def item(args) -> dict:
    depth, phaser, tests, level = args
    import confidence_score as cs
    import genome_common as common
    from l1_key_selection import matched_random
    from stratified_gate import strata
    from stratified_score import cells, grouped
    from w1_score import CAL
    from w1b_score import POOL, pooled

    from phasebudget.confidence import learn_then_test_by_group

    chroms = list(common.CONFIRMATORY)
    stratum = strata(chroms)
    pool = {"HG002": grouped(CAL, depth, phaser, chroms, stratum),
            **{s: grouped(f"hprc_{s}", depth, phaser, chroms, stratum) for s in POOL}}
    th = {beta: learn_then_test_by_group(pooled(pool), beta, cs.DELTA) for beta in cs.BETAS}
    del pool
    out = {"lines": [], "ratio": {}, "cells": []}
    for s in tests:
        labels = {c: load_labels((f"hprc_{s}", c)) for c in chroms}
        rows = dipcall_rows(f"hprc_{s}", depth, phaser, stratum)
        cls = [link_class(labels[r[4]].get(r[5], (None, None)),
                          labels[r[4]].get(r[6], (None, None))) for r in rows]
        kept = [r for r, c in zip(rows, cls) if c != "disagree"]
        key = [r for r, c in zip(rows, cls) if c == "key"]
        n_dis = sum(c == "disagree" for c in cls)
        k_dis = sum(r[3] for r, c in zip(rows, cls) if c == "disagree")
        out["lines"].append(f"  {s} {depth}x {phaser}: dipcall links {len(rows)}, disagree "
                            f"{n_dis} (wrong against dipcall {k_dis})")
        for name, groups in [(g, (g,)) for g in GROUPS] + [("hard", HARD)]:
            a = [r for r in kept if r[0] in groups]
            b = [r for r in key if r[0] in groups]
            ka, kb = sum(r[3] for r in a), sum(r[3] for r in b)
            km, nm = ka - kb, len(a) - len(b)
            ratio = ((ka / len(a)) / (kb / len(b)) if kb else
                     float("inf") if ka >= 3 else float("nan"))
            out["ratio"][(s, depth, phaser, name)] = ratio
            out["lines"].append(
                f"    {name:<9} links {len(a):>8} key share {len(b) / max(1, len(a)):.4f}; wrong "
                f"all {rate(ka, len(a))}, key {rate(kb, len(b))}, missing {rate(km, nm)}; "
                f"ratio all/key {ratio:.2f}")
        rand = matched_random(kept, key)
        for beta in cs.BETAS:
            for mode, rs in (("all", kept), ("key", key), ("random", rand)):
                for g, c in cells(th[beta], [r[:4] for r in rs], beta).items():
                    bad = bool(c[0]) and c[2] < level
                    out["cells"].append((mode, s, depth, phaser, beta, g, bad))
                    if bad or mode == "all":
                        out["lines"].append(
                            f"    {mode:<6} beta {beta:<6g} {g:<9} units {c[0]:>6} wrong {c[1]:>4} "
                            f"p {c[2]:.2e} trusted {c[4]:>8} wrong trusted {c[6]:>5}"
                            f"{'  VIOLATES' if bad else ''}")
    print(f"  {depth}x {phaser} done", flush=True)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", nargs=2, metavar=("CHROM", "OUT"))
    parser.add_argument("--dev", action="store_true")
    parser.add_argument("--test", action="store_true")
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    if args.labels:
        Path(args.labels[1]).write_text(json.dumps(key_labels(args.labels[0])))
        return
    if args.dev == args.test:
        raise SystemExit("give exactly one of --dev, --test")

    import confidence_score as cs
    import genome_common as common
    from w1b_score import HELD_OUT, POOL

    if args.test and not (ROOT / "docs" / "design" / "l8-hprc-key-selection.md").exists():
        raise SystemExit("commit docs/design/l8-hprc-key-selection.md before --test")
    started = time.time()
    chroms = list(common.CONFIRMATORY)
    tests = POOL if args.dev else HELD_OUT
    level = 0.05 / (len(tests) * CELLS_PER_GENOME)
    print(f"L8 key selection ({'DEV: the 20 W1b pool genomes, thresholds in sample' if args.dev else 'TEST: the 19 W1b held-out genomes'}; "
          f"docs/design/l8-hprc-key-selection.md). Rows: confidence_score rows against each "
          f"genome's HPRC release 2 dipcall truth alone; key = dipcall equal to 1000 Genomes "
          f"transmission at both variants; disagreeing links left out; thresholds W1b pooled "
          f"(HG002 v5.0q + 20 consensus-truth HPRC genomes); groups {GROUPS}; hard {HARD}; unit "
          f"100 kb; betas {cs.BETAS}; delta {cs.DELTA}; depths {cs.DEPTHS}; phasers "
          f"{cs.PHASERS}; {len(tests)} genomes; level 0.05/{len(tests) * CELLS_PER_GENOME} = "
          f"{level:.3e}; random control seed 0.", flush=True)

    jobs = [(f"hprc_{s}", c) for s in tests for c in chroms]
    with ThreadPoolExecutor(args.workers) as ex:
        list(ex.map(load_labels, jobs))
    os.environ["PHASEBUDGET_TRUTH"] = "dipcall"
    rjobs = [(f"hprc_{s}", d, ph) for s in tests for d in cs.DEPTHS for ph in cs.PHASERS]
    with ThreadPoolExecutor(args.workers) as ex:
        list(ex.map(lambda j: cs.load_rows(j[0], j[1], j[2], chroms), rjobs))
    del os.environ["PHASEBUDGET_TRUTH"]
    print(f"labels and rows cached ({time.time() - started:.0f} s)", flush=True)

    items = [(d, ph, tests, level) for d in cs.DEPTHS for ph in cs.PHASERS]
    with ProcessPoolExecutor(len(items)) as ex:
        res = list(ex.map(item, items))
    for r in res:
        print("\n".join(r["lines"]))
    ratio = {k: v for r in res for k, v in r["ratio"].items()}
    cells_ = [c for r in res for c in r["cells"]]
    print("\nSummary")
    for ph in cs.PHASERS:
        vals = [ratio[(s, "30", ph, "hard")] for s in tests]
        print(f"  30x {ph}: hard ratio all/key per genome "
              + ", ".join(f"{s} {v:.2f}" for s, v in zip(tests, vals))
              + f"; > 3 in {sum(v > 3 for v in vals)} of {len(tests)}; smallest {min(vals):.2f}")
    for mode in ("all", "key", "random"):
        bad = [c for c in cells_ if c[0] == mode and c[6]]
        by_g, by_s = {}, {}
        for c in bad:
            by_g[c[5]] = by_g.get(c[5], 0) + 1
            by_s[c[1]] = by_s.get(c[1], 0) + 1
        hard_genomes = {c[1] for c in bad if c[5] in (*HARD, "indel")}
        print(f"  {mode}: violating cells {len(bad)}; by group {by_g}; genomes with a violating "
              f"segdup, lowmap, tandem or indel cell {len(hard_genomes)} of {len(tests)}; by "
              f"genome {by_s}")
    print(f"runtime {time.time() - started:.0f} s")


if __name__ == "__main__":
    main()
