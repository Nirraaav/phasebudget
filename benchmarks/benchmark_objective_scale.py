"""Scale for the selection-objective ratio: WhatsHap, Fischer-Marschall and random.

The objective is the additive one `select_optimal` maximizes: the sum over
selected reads of `total_weight_score`, the total natural log-odds weight of a
phase-informative read's called alleles, under a per-variant coverage cap. The
ratio `sum objective / sum optimum` pooled over chromosomes is therefore the
fraction of the optimum's allele log-odds weight a selector keeps. The second
ratio, `calls`, counts called alleles instead: the number of called alleles in
the selected reads over the number in the optimum's selection. The optimum does
not maximize that count, so `calls` can exceed 1.

WhatsHap and fm_paper objectives and the optimum are read from the stored run
files where present (`lib_{depth}_{pipeline}_c{cap}.json.gz`); every selection
is also recomputed here, which the call counts need, and the recomputed
objective is checked against the stored one. Random selection is
`select_random` (phase-informative reads in shuffled order, kept while the cap
allows) under the same base weights. No phasing is run.
"""

from __future__ import annotations

import argparse
import datetime
import gzip
import importlib.util
import json
import os
import time
from multiprocessing import Pool
from pathlib import Path

from phasebudget.selection import (read_weights, select_fischer_marschall, select_optimal,
                                   select_random)

HERE = Path(__file__).resolve().parent
CORRECTION = datetime.datetime(2026, 9, 30, 19, 31)


def _load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


common = _load("genome_common")


def available_gb() -> float:
    with open("/proc/meminfo") as handle:
        for line in handle:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1024 ** 2
    raise RuntimeError("no MemAvailable in /proc/meminfo")


def run_file(chrom, depth, pipeline, cap):
    return common.chrom_dir(chrom) / "runs" / f"lib_{depth}_{pipeline}_c{cap}.json.gz"


def stored(chrom, depth, pipeline, cap):
    path = run_file(chrom, depth, pipeline, cap)
    if not path.exists():
        return None
    with gzip.open(path, "rt") as handle:
        record = json.load(handle)
    stale = datetime.datetime.fromtimestamp(path.stat().st_mtime) < CORRECTION
    return record["objective"], record["optimum"], stale


def task(args):
    chrom, depth, caps, seeds, min_gb, cache = args
    done = Path(cache) / f"{chrom}_{depth}.json" if cache else None
    if done and done.exists():
        rows = json.loads(done.read_text())
        if [r["cap"] for r in rows] == list(caps):
            return rows
    # The machine runs other long jobs: wait rather than start a chromosome short of memory.
    while available_gb() < min_gb:
        time.sleep(30)
    extractor = _load("extract_fragments")
    pipelines = _load("hg002_pipelines")
    matrix, positions, _, _ = extractor.load_fragments(
        str(common.fragments_path(chrom, depth, "realign")))
    matrix = extractor.phred_to_log_odds(matrix)
    weights = read_weights(matrix)
    calls = {f.read_id: f.n_observed for f in matrix}
    distinct = sorted({w for f in matrix for _, _, w in f.observed()})[:5]
    rows = []
    for cap in caps:
        started = time.time()
        optimum = select_optimal(matrix, cap, weights)
        row = {"chrom": chrom, "depth": depth, "cap": cap, "optimum": optimum.objective,
               "optimum_calls": sum(calls[r] for r in optimum.selected), "stale": 0,
               "mismatch": 0.0, "distinct": distinct}
        found = stored(chrom, depth, "optimal", cap)
        if found:
            row["mismatch"] = max(row["mismatch"], abs(found[1] - optimum.objective)
                                  / optimum.objective)
        selections = {
            "whatshap": pipelines.whatshap_selection(matrix, cap, positions),
            "fm_paper": select_fischer_marschall(matrix, cap, weights),
        }
        for seed in seeds:
            selections[f"random{seed}"] = select_random(matrix, cap, seed=seed, weights=weights)
        for name, selection in selections.items():
            objective = selection.objective
            source = "computed"
            found = stored(chrom, depth, name, cap) if not name.startswith("random") else None
            if found:
                row["mismatch"] = max(row["mismatch"], abs(found[0] - objective) / optimum.objective,
                                      abs(found[1] - optimum.objective) / optimum.objective)
                row["stale"] += found[2]
                objective, source = found[0], "run file"
            row[name] = objective
            row[name + "_calls"] = sum(calls[r] for r in selection.selected)
            row[name + "_source"] = source
        row["seconds"] = round(time.time() - started, 1)
        rows.append(row)
    if done:
        done.write_text(json.dumps(rows))
    return rows


def pooled(rows, name, key=""):
    return sum(r[name + key] for r in rows) / sum(r["optimum" + key] for r in rows)


def worst(rows, name, key=""):
    row = min(rows, key=lambda r: r[name + key] / r["optimum" + key])
    return row[name + key] / row["optimum" + key], row["chrom"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chromosomes", nargs="+", default=list(common.CONFIRMATORY))
    parser.add_argument("--depths", nargs="+", default=["full", "10"])
    parser.add_argument("--caps", nargs="+", type=int, default=[3, 5, 10, 15])
    parser.add_argument("--seeds", nargs="+", type=int, default=[0, 1, 2, 3, 4])
    parser.add_argument("--processes", type=int, default=2)
    parser.add_argument("--min-available-gb", type=float, default=4.0)
    parser.add_argument("--cache", help="directory for per-chromosome results, so a restart "
                        "resumes; the random seeds are part of the result, so clear it when "
                        "changing them")
    args = parser.parse_args()
    if available_gb() < args.min_available_gb:
        raise SystemExit(f"available memory {available_gb():.1f} GB, under "
                         f"{args.min_available_gb} GB")

    started = time.time()
    # Largest chromosomes first so the two workers finish together.
    if args.cache:
        os.makedirs(args.cache, exist_ok=True)
    jobs = sorted(((c, d, args.caps, args.seeds, args.min_available_gb, args.cache)
                   for c in args.chromosomes for d in args.depths),
                  key=lambda j: -common.fragments_path(j[0], j[1], "realign").stat().st_size)
    # One task per worker process, so a large chromosome's memory is returned before the next.
    with Pool(args.processes, maxtasksperchild=1) as pool:
        rows = [row for result in pool.imap_unordered(task, jobs) for row in result]
    elapsed = time.time() - started

    randoms = [f"random{s}" for s in args.seeds]
    print("# Objective scale: selector objective over the certified optimum (select_optimal)")
    print(f"# platform {common.PLATFORM}, fragments *_realign, weights logodds "
          f"(phred_to_log_odds), score total_weight_score")
    print(f"# chromosomes {' '.join(args.chromosomes)}")
    print(f"# depths {' '.join(args.depths)}; caps {' '.join(map(str, args.caps))}; "
          f"random seeds {' '.join(map(str, args.seeds))}; processes {args.processes}")
    print("# objective = sum over selected reads of the total natural log-odds weight of the")
    print("#   read's called alleles (phase-informative reads only); 'weight' = pooled")
    print("#   sum objective / sum optimum = fraction of the optimum's allele log-odds weight kept")
    print("# 'calls' = pooled called alleles in the selection / called alleles in the optimum's")
    print("#   selection; the optimum does not maximize this count, so it can exceed 1")
    print("# whatshap = WhatsHap 2.8 compiled readselection (hg002_pipelines.whatshap_selection);")
    print("#   fm_paper = select_fischer_marschall (published procedure, bridging on);")
    print("#   random = select_random, mean [min, max] over seeds of the pooled ratio")
    print("# whatshap/fm_paper weight ratios and optima come from run files where present, else")
    print("#   computed here; every selection is also recomputed for the call counts.")
    sources = {n: sum(r[n + "_source"] == "run file" for r in rows) for n in ("whatshap", "fm_paper")}
    print(f"# rows from run files: whatshap {sources['whatshap']}/{len(rows)}, "
          f"fm_paper {sources['fm_paper']}/{len(rows)}")
    print(f"# run files older than {CORRECTION:%Y-%m-%d %H:%M} (pre-correction HiFi library "
          f"run): {sum(r['stale'] for r in rows)} used; objectives do not depend on the "
          f"corrected reporting, so both are valid")
    print(f"# max |stored - recomputed| / optimum over objectives and optima: "
          f"{max(r['mismatch'] for r in rows):.2e}")
    distinct = sorted({w for r in rows for w in r["distinct"]})
    print(f"# distinct called-allele weights (first 5): {distinct[:5]}; with one value the"
          f" weight and calls ratios coincide")
    print(f"# wall time {elapsed / 60:.1f} min, including waits for {args.min_available_gb} GB "
          f"available memory before each chromosome; per-task seconds are in the last column")
    print()
    print(f"{'depth':>5} {'cap':>3} {'metric':>6} | {'whatshap':>8} {'worst':>13} | "
          f"{'fm_paper':>8} {'worst':>13} | {'random mean [min, max]':>24} {'worst':>13}")
    for depth in args.depths:
        for cap in args.caps:
            group = [r for r in rows if r["depth"] == depth and r["cap"] == cap]
            for key, label in (("", "weight"), ("_calls", "calls")):
                cells = []
                for name in ("whatshap", "fm_paper"):
                    value, chrom = worst(group, name, key)
                    cells.append(f"{pooled(group, name, key):8.4f} {value:7.4f} {chrom:>5}")
                values = [pooled(group, n, key) for n in randoms]
                per_chrom = [(sum(r[n + key] for n in randoms) / len(randoms) / r["optimum" + key],
                              r["chrom"]) for r in group]
                value, chrom = min(per_chrom)
                mean = sum(values) / len(values)
                print(f"{depth:>5} {cap:>3} {label:>6} | {cells[0]} | {cells[1]} | "
                      f"{mean:6.4f} [{min(values):.4f}, {max(values):.4f}] {value:7.4f} {chrom:>5}")
    print()
    print("# worst = lowest per-chromosome ratio and its chromosome (random: seed mean)")
    print()
    print("# per chromosome: weight ratios, then call ratios; random is the seed mean")
    print(f"{'chrom':>5} {'depth':>5} {'cap':>3} {'optimum':>12} {'wh':>7} {'fm':>7} {'rnd':>7} "
          f"{'wh_c':>7} {'fm_c':>7} {'rnd_c':>7} {'sec':>6}")
    order = {c: i for i, c in enumerate(args.chromosomes)}
    for r in sorted(rows, key=lambda r: (args.depths.index(r["depth"]), r["cap"],
                                         order[r["chrom"]])):
        rnd = sum(r[n] for n in randoms) / len(randoms) / r["optimum"]
        rnd_c = sum(r[n + "_calls"] for n in randoms) / len(randoms) / r["optimum_calls"]
        print(f"{r['chrom']:>5} {r['depth']:>5} {r['cap']:>3} {r['optimum']:12.1f} "
              f"{r['whatshap'] / r['optimum']:7.4f} {r['fm_paper'] / r['optimum']:7.4f} "
              f"{rnd:7.4f} {r['whatshap_calls'] / r['optimum_calls']:7.4f} "
              f"{r['fm_paper_calls'] / r['optimum_calls']:7.4f} {rnd_c:7.4f} {r['seconds']:6.1f}")


if __name__ == "__main__":
    main()
