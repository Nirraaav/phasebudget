"""Score the registered genome-wide HG002 run and decide G1 to G10.

`docs/design/genome-run-plan.md` fixes the analysis and its amendment 1 the
details the plan left open. This script reads the run files `genome_run.py`
left in `data/genome`, the three truths of `genome_common.load_truths`, and
prints:

  0. the data: depth of each slice, truth sizes, the recalibration table;
  1. per depth, cap and pipeline, pooled over the confirmatory chromosomes:
     switch error against the filtered truth (the primary), flips, phased
     fraction, block N50 and QAN50, the objective over the certified optimum
     (pooled and the worst chromosome), switch error against the unfiltered
     transmission truth and the hifiasm phase transfer, runtime and peak
     memory;
  2. the paired sign tests, with the difference in error rate per 10,000
     common pairs and its 1 Mb window-bootstrap interval, against the
     filtered truth and, beside it, the unfiltered one;
  3. G1 to G10 with their registered decision rules.

Per-chromosome rows go to a TSV file. Chromosome 20 is the pilot and is not
read here. `--truth v5` puts GIAB v5.0q in the filtered truth's place, as
amendment 2 fixes; that analysis is post hoc.
"""

from __future__ import annotations

import argparse
import gzip
import json
import multiprocessing
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402
import genome_run  # noqa: E402

from phasebudget.bootstrap import statistic_interval, window_counts  # noqa: E402
from phasebudget.metrics import (block_lengths, comparable_pairs, flip_errors, n50,  # noqa: E402
                                 phased_fraction, qan50_lengths, switch_positions)
from phasebudget.paired import common_pairs, sign_test, window_discordance  # noqa: E402
from permutation import window_sign_flip_test  # noqa: E402

score = common._module("score_hg002.py")

WINDOW = 1_000_000
RESAMPLES = 2000
SEED = 0
TRUTHS = ("filtered", "transmission", "hifiasm", "v5")
# The truth every test and decision uses: "filtered" as registered, "v5" under amendment 2.
PRIMARY = "filtered"
CHROMS = common.CONFIRMATORY
DEPTHS = common.DEPTHS
CAPS = genome_run.CAPS
TOOLS = ("whatshap phase", "HapCUT2")
# (first, second, caps): every paired comparison the hypotheses make, plus
# `whatshap` against `optimal` at caps 10 and 15 for G10's description.
COMPARISONS = (
    ("whatshap", "optimal", (3, 5, 10, 15)),
    ("merge_optimal_1e3", "optimal", (3, 5)),
    ("two_level", "optimal", (3,)),
    ("connected", "optimal", (3,)),
    ("loop_oriented", "loop", (3,)),
    ("loop", "optimal", (3, 5)),
)
if common.PLATFORM == "ont":
    # `docs/design/ont-run-plan.md`, O1 to O7.
    TOOLS = ("whatshap phase", "HapCUT2", "LongPhase")
    COMPARISONS = (
        ("whatshap", "optimal", (3, 5, 10, 15)),
        ("merge_optimal_1e3", "optimal", (3, 5)),
        ("two_level", "optimal", (3,)),
        ("error_optimal_1e2", "optimal", (3,)),
    )
if common.PLATFORM == "hg005":
    # `docs/design/hg005-run-plan.md`, H5-1 to H5-3.
    TOOLS = ("whatshap phase",)
    COMPARISONS = (
        ("whatshap", "optimal", (3, 5, 10, 15)),
        ("two_level", "optimal", (3,)),
        ("connected", "optimal", (3,)),
    )
KINDS = ("lib", "bq") if common.PLATFORM == "hifi" else ("lib",)

TRUTH: dict = {}
POSITIONS: dict = {}


def load_all_truths() -> None:
    for chrom in CHROMS:
        variants = score.variant_table(common.hetsnv_path(chrom))
        positions = [p - 1 for p, _, _ in variants]
        POSITIONS[chrom] = positions
        TRUTH[chrom] = dict(zip(TRUTHS, (*common.load_truths(chrom, positions),
                                         common.v5_truth(chrom, variants))))


def pipelines_of(kind: str, cap) -> tuple:
    if kind == "lib":
        return genome_run.LIBRARY + (genome_run.LOOPS if cap in genome_run.LOOP_CAPS else ())
    if kind == "bq":
        return ("optimal", "whatshap")
    return TOOLS


def _existing(path):
    """`path`, or its gzipped copy: finished tool outputs are compressed to save disk."""
    gz = path.with_name(path.name + ".gz")
    return path if path.exists() or not gz.exists() else gz


def load_phasing(kind, chrom, depth, pipeline, cap):
    """`(phasing, record)`; the record holds objective, time and memory."""
    d = common.chrom_dir(chrom)
    if kind in ("lib", "bq"):
        with gzip.open(d / "runs" / f"{kind}_{depth}_{pipeline}_c{cap}.json.gz", "rt") as handle:
            record = json.load(handle)
        return score.run_phasing(record), record
    variants = score.variant_table(common.hetsnv_path(chrom))
    if pipeline == "whatshap phase":
        phasing = score.whatshap_phasing(_existing(d / "external" / f"whatshap_{depth}.vcf"), variants)
        log = genome_run.log_of(f"{chrom}:whatshap:{depth}")
    elif pipeline == "LongPhase":
        phasing = score.whatshap_phasing(_existing(d / "external" / f"longphase_{depth}.vcf"), variants)
        log = genome_run.log_of(f"{chrom}:longphase:{depth}")
    else:
        phasing = score.hapcut2_phasing(d / "external" / f"hapcut2_{depth}.blocks", variants)
        log = genome_run.log_of(f"{chrom}:hapcut2:{depth}")
    seconds, peak = score.time_log(log)
    return phasing, {"objective": None, "optimum": None, "seconds": seconds, "peak_mb": peak}


def score_group(group):
    """Everything one (kind, depth, cap) contributes, summed over chromosomes."""
    kind, depth, cap = group
    names = pipelines_of(kind, cap)
    comparisons = [(a, b) for a, b, caps in COMPARISONS if kind == "lib" and cap in caps]
    pooled = {n: {"switches": {t: 0 for t in TRUTHS}, "pairs": {t: 0 for t in TRUTHS},
                  "flips": 0, "phased": 0, "variants": 0, "blocks": [], "qblocks": [],
                  "objective": 0.0, "optimum": 0.0, "worst": None, "seconds": 0.0,
                  "peak_mb": 0.0, "windows": {PRIMARY: [], "transmission": []}}
              for n in names}
    paired = {(a, b, t): {"pairs": 0, "first": 0, "second": 0, "both": 0, "windows": []}
              for a, b in comparisons for t in tuple(dict.fromkeys((PRIMARY, "transmission")))}
    rows = []
    signals = []
    for chrom in CHROMS:
        truth, positions = TRUTH[chrom], POSITIONS[chrom]
        phasings = {}
        for name in names:
            phasing, record = load_phasing(kind, chrom, depth, name, cap)
            if len(phasing) != len(positions):
                raise ValueError(f"{chrom} {kind} {depth} {name} {cap}: wrong variant count")
            phasings[name] = phasing
            p = pooled[name]
            row = {"chrom": chrom, "kind": kind, "depth": depth, "cap": cap, "pipeline": name}
            for t in TRUTHS:
                s = len(switch_positions(truth[t], phasing))
                c = comparable_pairs(phasing, truth[t])
                p["switches"][t] += s
                p["pairs"][t] += c
                row[f"switches_{t}"], row[f"pairs_{t}"] = s, c
            for t in tuple(dict.fromkeys((PRIMARY, "transmission"))):
                p["windows"][t].extend(window_counts(truth[t], phasing, positions, WINDOW))
            flips = len(flip_errors(truth[PRIMARY], phasing))
            phased = round(phased_fraction(phasing) * len(phasing))
            lengths = block_lengths(phasing, positions)
            qlengths = qan50_lengths(truth[PRIMARY], phasing, positions)
            p["flips"] += flips
            p["phased"] += phased
            p["variants"] += len(phasing)
            p["blocks"].extend(lengths)
            p["qblocks"].extend(qlengths)
            p["seconds"] += record["seconds"]
            p["peak_mb"] = max(p["peak_mb"], record["peak_mb"])
            ratio = None
            if record.get("objective") is not None:
                p["objective"] += record["objective"]
                p["optimum"] += record["optimum"]
                ratio = record["objective"] / record["optimum"]
                p["worst"] = ratio if p["worst"] is None else min(p["worst"], ratio)
            row.update(flips=flips, phased=phased / len(phasing), n50_kb=n50(lengths) / 1000,
                       qan50_kb=n50(qlengths) / 1000, objective_ratio=ratio,
                       seconds=record["seconds"], peak_mb=record["peak_mb"])
            rows.append(row)
            if kind == "bq" and depth == "full" and cap == 5 and name == "optimal" \
                    and chrom in common.EVEN:
                n_windows = len(window_counts(truth[PRIMARY], phasing, positions, WINDOW))
                for b, called, predicted, post_margin, mec_margin in record["signals"]:
                    f = truth[PRIMARY]
                    if f[b] is None or f[b + 1] is None:
                        continue
                    k = (positions[b + 1] - positions[0]) // WINDOW
                    signals.append((chrom, k, n_windows, called != (f[b] ^ f[b + 1]), predicted,
                                    post_margin, mec_margin))
        for a, b in comparisons:
            for t in tuple(dict.fromkeys((PRIMARY, "transmission"))):
                found = common_pairs(truth[t], phasings[a], phasings[b])
                q = paired[(a, b, t)]
                q["pairs"] += len(found)
                q["first"] += sum(1 for _, x, y in found if x and not y)
                q["second"] += sum(1 for _, x, y in found if y and not x)
                q["both"] += sum(1 for _, x, y in found if x and y)
                # The same windows as `window_counts`, so the list lines up across pipelines.
                q["windows"].extend(window_discordance(found, positions, WINDOW))
    for p in pooled.values():
        p["n50_kb"] = n50(p.pop("blocks")) / 1000
        p["qan50_kb"] = n50(p.pop("qblocks")) / 1000
    return group, pooled, paired, rows, signals


def rate(counts, picks):
    pairs = sum(counts[k][1] for k in picks)
    return None if pairs == 0 else sum(counts[k][0] for k in picks) / pairs


def difference_interval(windows):
    """First's error rate minus second's per 10,000 common pairs, window bootstrap."""
    def statistic(picks):
        pairs = sum(windows[k][0] for k in picks)
        if pairs == 0:
            return None
        return 1e4 * (sum(windows[k][1] for k in picks) - sum(windows[k][2] for k in picks)) / pairs
    return statistic_interval(statistic, len(windows), resamples=RESAMPLES, seed=SEED)


def fmt(interval, digits=2):
    return (f"{interval.estimate:+.{digits}f} [{interval.low:+.{digits}f}, "
            f"{interval.high:+.{digits}f}]")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--tsv", default="benchmarks/results/runs/genome_per_chromosome.tsv")
    parser.add_argument("--truth", choices=("filtered", "v5", "transmission"),
                        default={"ont": "v5", "hg005": "transmission"}.get(common.PLATFORM,
                                                                          "filtered"),
                        help="the truth for every test and decision; on HiFi v5 is amendment 2, "
                             "post hoc; on ONT v5 is registered")
    args = parser.parse_args()
    global PRIMARY
    PRIMARY = args.truth
    started = time.time()
    load_all_truths()

    print("genome-wide HG002 scoring (docs/design/genome-run-plan.md, amendment 1)")
    registered = PRIMARY == {"ont": "v5", "hg005": "transmission"}.get(common.PLATFORM,
                                                                        "filtered")
    print(f"platform {common.PLATFORM}; primary truth: {PRIMARY}" +
          (" (registered)" if registered else " (POST HOC, beside the registered analysis)"))
    print(f"confirmatory chromosomes {len(CHROMS)}: {', '.join(CHROMS)}; chromosome 20 is the "
          "pilot and is not scored here")
    print(f"depths {', '.join(DEPTHS)}; caps {', '.join(map(str, CAPS))}; loop pipelines at "
          f"caps {', '.join(map(str, genome_run.LOOP_CAPS))}")
    print(f"bootstrap: {WINDOW // 1_000_000} Mb windows cut per chromosome and pooled, "
          f"{RESAMPLES} resamples, seed {SEED}; sign tests exact, two-sided")
    status = json.loads((common.DATA / "status.json").read_text())
    print(f"run status {status['time']}: {status['done']} jobs done, "
          f"{len(status['failed'])} failed")

    # 0. Data.
    print()
    print("0. data")
    variants = 0
    counts = {t: 0 for t in TRUTHS}
    removed = 0
    for chrom in CHROMS:
        depth = json.loads((common.chrom_dir(chrom) / "depth.json").read_text())
        skips = [d for d in DEPTHS if genome_run.skipped(chrom, d)]
        truth = TRUTH[chrom]
        n = len(POSITIONS[chrom])
        variants += n
        for t in TRUTHS:
            counts[t] += sum(x is not None for x in truth[t])
        cut = sum(a is not None and b is None
                  for a, b in zip(truth["transmission"], truth["filtered"]))
        removed += cut
        v5, tr = truth["v5"], truth["transmission"]
        both = [j for j in range(n) if v5[j] is not None and tr[j] is not None]
        same = sum(1 for p, q in zip(both, both[1:]) if (v5[p] ^ v5[q]) == (tr[p] ^ tr[q]))
        print(f"{chrom:>5} depth {depth['depth']:5.1f}, {n:>7} het SNVs, filtered truth "
              f"{sum(x is not None for x in truth['filtered']):>7} (removed {cut:>4}), "
              f"hifiasm {sum(x is not None for x in truth['hifiasm']):>7}, v5 "
              f"{sum(x is not None for x in v5):>7}; v5 and transmission relations agree at "
              f"{same}/{max(len(both) - 1, 0)} consecutive shared sites; skipped depths: "
              f"{', '.join(skips) or 'none'}")
    print(f"total {variants} het SNVs; truth sites: " +
          ", ".join(f"{t} {counts[t]}" for t in TRUTHS) + f"; filter removed {removed}")
    table = json.loads((common.DATA / "recalibration.json").read_text()) \
        if common.PLATFORM == "hifi" else {"chromosomes": [], "bins": []}
    if table["bins"]:
        print(f"recalibration fitted on {', '.join(table['chromosomes'])}: "
              f"{len(table['bins'])} bins")
    for b in table["bins"]:
        print(f"  quality {b['low']:>3}-{b['high']:<3} alleles {b['alleles']:>10} wrong "
              f"{b['wrong']:>6} error {b['error']:.6f} weight {b['weight']:.3f}")

    groups = [(k, d, c) for k in KINDS for d in DEPTHS for c in CAPS] + \
        [("tool", d, None) for d in DEPTHS]
    with multiprocessing.get_context("fork").Pool(args.jobs) as pool:
        results = pool.map(score_group, groups, chunksize=1)
    pooled, paired, rows, signals = {}, {}, [], []
    for group, p, q, r, s in results:
        for name, value in p.items():
            pooled[(*group, name)] = value
        for (a, b, t), value in q.items():
            paired[(group[1], group[2], a, b, t)] = value
        rows.extend(r)
        signals.extend(s)

    # Per-chromosome TSV.
    tsv = Path(args.tsv)
    tsv.parent.mkdir(parents=True, exist_ok=True)
    columns = list(rows[0])
    with open(tsv, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join("" if row[c] is None else
                                   f"{row[c]:.6g}" if isinstance(row[c], float) else str(row[c])
                                   for c in columns) + "\n")
    print(f"per-chromosome rows: {tsv} ({len(rows)} rows)")

    # 1. Pooled metrics.
    OTHER = "v5" if PRIMARY == "filtered" else "filtered"
    header = (f"{'depth':>5} {'cap':>4} {'pipeline':>20} {'weights':>8} {'switch':>8} "
              f"{'sw/pairs':>15} {'flips':>5} {'phased':>6} {'N50 kb':>7} {'QAN50 kb':>8} "
              f"{'obj/opt':>7} {'worst':>7} {'unfilt':>8} {'hifiasm':>8} {OTHER:>8} {'hours':>6} "
              f"{'peak MB':>7}")
    print()
    print(f"1. pooled over the confirmatory chromosomes; switch against the {PRIMARY} truth, "
          f"unfilt, hifiasm and {OTHER} against the others; obj/opt is the pooled ratio and worst "
          "the lowest chromosome, n/a where reads are merged or for external tools")
    print(header)
    print("-" * len(header))
    for depth in DEPTHS:
        for cap in CAPS + (None,):
            for kind in KINDS + ("tool",):
                if (kind == "tool") != (cap is None):
                    continue
                for name in pipelines_of(kind, cap):
                    p = pooled[(kind, depth, cap, name)]
                    sw = {t: p["switches"][t] / p["pairs"][t] if p["pairs"][t] else 0.0
                          for t in TRUTHS}
                    ratio = f"{p['objective'] / p['optimum']:.4f}" if p["optimum"] else "n/a"
                    worst = f"{p['worst']:.4f}" if p["worst"] is not None else "n/a"
                    weights = {"lib": "logodds", "bq": "recal", "tool": "tool"}[kind]
                    print(f"{depth:>5} {cap if cap else '-':>4} {name:>20} {weights:>8} "
                          f"{sw[PRIMARY]:>8.6f} "
                          f"{p['switches'][PRIMARY]:>6}/{p['pairs'][PRIMARY]:<8} "
                          f"{p['flips']:>5} {p['phased'] / p['variants']:>6.3f} "
                          f"{p['n50_kb']:>7.0f} {p['qan50_kb']:>8.0f} {ratio:>7} {worst:>7} "
                          f"{sw['transmission']:>8.6f} {sw['hifiasm']:>8.6f} {sw[OTHER]:>8.6f} "
                          f"{p['seconds'] / 3600:>6.2f} {p['peak_mb']:>7.0f}")
            print()

    # 2. Paired tests.
    print("2. paired sign tests on common pairs; first/second = pairs where only that "
          "pipeline is wrong; diff = first's rate minus second's per 10,000 common pairs, "
          f"with its window-bootstrap interval ({PRIMARY} truth); flip p = window sign-flip "
          "test (benchmarks/permutation.py, added after the mock review, no decision uses it); "
          "unfiltered truth beside")
    print(f"{'first':>18} {'second':>9} {'depth':>5} {'cap':>3} {'pairs':>8} {'first':>6} "
          f"{'second':>6} {'both':>5} {'p':>9} {'diff per 10k [95%]':>26} {'flip p':>9} | "
          f"{'unfilt first/second':>19} {'p':>9}")
    tests = {}
    for a, b, caps in COMPARISONS:
        for cap in caps:
            for depth in DEPTHS:
                q = paired[(depth, cap, a, b, PRIMARY)]
                u = paired[(depth, cap, a, b, "transmission")]
                pv = sign_test(q["first"], q["first"] + q["second"])
                pu = sign_test(u["first"], u["first"] + u["second"])
                interval = difference_interval(q["windows"])
                flip = window_sign_flip_test(q["windows"]).p_value
                tests[(a, b, depth, cap)] = (q, pv, interval)
                print(f"{a:>18} {b:>9} {depth:>5} {cap:>3} {q['pairs']:>8} {q['first']:>6} "
                      f"{q['second']:>6} {q['both']:>5} {pv:>9.2e} {fmt(interval):>26} "
                      f"{flip:>9.2e} | {u['first']:>9}/{u['second']:<9} {pu:>9.2e}")
        print()

    # 3. Hypotheses.
    def favours_second(key, level):
        q, pv, _ = tests[key]
        return pv < level and q["first"] > q["second"]

    def favours_first(key, level):
        q, pv, _ = tests[key]
        return pv < level and q["first"] < q["second"]

    print("3. hypotheses")
    if common.PLATFORM == "hg005":
        decide_hg005(pooled, tests, favours_first, favours_second)
        print(f"\nscored in {(time.time() - started) / 60:.1f} min")
        return
    if common.PLATFORM == "ont":
        decide_ont(pooled, tests, favours_first, favours_second)
        print(f"\nscored in {(time.time() - started) / 60:.1f} min")
        return
    p = pooled[("lib", "full", 15, "whatshap")]
    g1 = p["objective"] / p["optimum"]
    print(f"G1 whatshap / certified optimum, full depth, cap 15: pooled {g1:.4f}, worst "
          f"chromosome {p['worst']:.4f} -> {'supported' if g1 < 0.95 else 'not supported'} "
          "(rule: below 0.95)")

    level = 0.05 / 8
    g2 = {(d, c): tests[("whatshap", "optimal", d, c)] for d in DEPTHS for c in (3, 5)}
    low = all(favours_second(("whatshap", "optimal", d, c), level)
              for d in ("15", "10") for c in (3, 5))
    none = not any(g2[k][1] < level for k in g2)
    print(f"G2 whatshap vs optimal, caps 3 and 5, level {level:.5f}: significant for optimal at "
          + ", ".join(f"{d}x c{c}" for (d, c) in g2
                      if favours_second(("whatshap", "optimal", d, c), level))
          + " -> " + ("supported" if low else "refuted" if none else
                      "not supported, not refuted")
          + " (rule: all four tests at 10- and 15-fold for optimal; refuted if none significant)")

    g3 = all(favours_first(("merge_optimal_1e3", "optimal", d, c), level)
             for d in ("15", "10") for c in (3, 5))
    print("G3 merge_optimal_1e3 vs optimal, caps 3 and 5, level 0.00625: significant for merge "
          "at " + ", ".join(f"{d}x c{c}" for d in DEPTHS for c in (3, 5)
                           if favours_first(("merge_optimal_1e3", "optimal", d, c), level))
          + "; significant for optimal at "
          + ", ".join(f"{d}x c{c}" for d in DEPTHS for c in (3, 5)
                      if favours_second(("merge_optimal_1e3", "optimal", d, c), level))
          + f" -> {'supported' if g3 else 'not supported'}")

    level4 = 0.05 / 4
    for label, name in (("G4", "two_level"), ("G5", "connected")):
        gains = {d: pooled[("lib", d, 3, name)]["n50_kb"] / pooled[("lib", d, 3, "optimal")]["n50_kb"]
                 for d in DEPTHS}
        against = [d for d in DEPTHS if favours_second((name, "optimal", d, 3), level4)]
        ok = all(g >= 1.25 for g in gains.values()) and not against
        print(f"{label} {name} vs optimal, cap 3: block N50 ratio "
              + ", ".join(f"{d}x {g:.2f}" for d, g in gains.items())
              + f"; sign tests for optimal at level {level4}: {', '.join(against) or 'none'} -> "
              + ("supported" if ok else "not supported")
              + " (rule: N50 up at least 25 per cent at every depth, no test for optimal)")

    fixed = [d for d in DEPTHS if favours_first(("loop_oriented", "loop", d, 3), level4)]
    g6 = all(d in fixed for d in ("15", "10"))
    print(f"G6 loop_oriented vs loop, cap 3, level {level4}: significant for the repair at "
          f"{', '.join(fixed) or 'no depth'} -> {'supported' if g6 else 'not supported'} "
          "(rule: at 10- and 15-fold)")
    loop_better = [(d, c) for d in DEPTHS for c in (3, 5)
                   if favours_first(("loop", "optimal", d, c), level)]
    loop_worse = [(d, c) for d in DEPTHS for c in (3, 5)
                  if favours_second(("loop", "optimal", d, c), level)]
    print("G6 report, loop vs optimal two-sided at level 0.00625: loop better at "
          + (", ".join(f"{d}x c{c}" for d, c in loop_better) or "none")
          + "; loop worse at " + (", ".join(f"{d}x c{c}" for d, c in loop_worse) or "none"))

    if not loop_better:
        print("G7 without its premise: loop is not significantly better than optimal at any "
              "cap and depth")
    for d, c in loop_better:
        o, l, s = (pooled[("lib", d, c, n)]["windows"][PRIMARY]
                   for n in ("optimal", "loop", "shuffled"))

        def kept(picks, o=o, l=l, s=s):
            ro, rl, rs = rate(o, picks), rate(l, picks), rate(s, picks)
            if None in (ro, rl, rs) or ro == rl:
                return None
            return (ro - rs) / (ro - rl)

        interval = statistic_interval(kept, len(o), resamples=RESAMPLES, seed=SEED)
        print(f"G7 {d}x cap {c}: share of the loop's reduction kept by shuffled "
              f"{fmt(interval, 3)} -> {'supported' if interval.estimate >= 0.8 else 'not supported'}"
              " (rule: at least 0.8)")

    # G8 and G9.
    offsets, total = {}, 0
    for chrom in common.EVEN:
        n_windows = next((r[2] for r in signals if r[0] == chrom), 0)
        offsets[chrom] = total
        total += n_windows
    per_window = [[0.0, 0.0] for _ in range(total)]
    for chrom, k, _, wrong, predicted, _, _ in signals:
        per_window[offsets[chrom] + k][0] += wrong
        per_window[offsets[chrom] + k][1] += predicted

    def calibration(picks):
        predicted = sum(per_window[k][1] for k in picks)
        return None if predicted == 0 else sum(per_window[k][0] for k in picks) / predicted

    observed = sum(r[3] for r in signals)
    predicted = sum(r[4] for r in signals)
    interval = statistic_interval(calibration, total, resamples=RESAMPLES, seed=SEED)
    bins = [[0, 0.0, 0.0] for _ in range(10)]
    for r in signals:
        k = min(int(r[4] * 10), 9)
        bins[k][0] += 1
        bins[k][1] += r[3]
        bins[k][2] += r[4]
    ece = sum(abs(w - p) for _, w, p in bins) / len(signals)
    calibrated = interval.low <= 1 <= interval.high and ece < 0.01
    print(f"G8 recalibrated optimal, full depth, cap 5, even chromosomes: {len(signals)} "
          f"boundaries, observed wrong {observed}, predicted {predicted:.1f}; observed / "
          f"predicted {fmt(interval)}; ECE {ece:.5f} -> "
          f"{'calibrated' if calibrated else 'not calibrated'} (rule: interval contains 1 and "
          "ECE below 0.01)")
    for k, (count, wrong, pred) in enumerate(bins):
        if count:
            print(f"  predicted in [{k / 10:.1f}, {(k + 1) / 10:.1f}): {count:>7} boundaries, "
                  f"observed rate {wrong / count:.5f}, mean predicted {pred / count:.5f}")
    bad = [r[5] for r in signals if r[3]]
    good = [r[5] for r in signals if not r[3]]
    value = score.auc(bad, good)
    mec = score.auc([r[6] for r in signals if r[3]], [r[6] for r in signals if not r[3]])
    print(f"G9 AUC of the posterior margin for a wrong call: {value:.3f} over {len(bad)} wrong "
          f"and {len(good)} right calls -> {'supported' if value > 0.6 else 'not supported'} "
          f"(rule: above 0.6); MEC margin {mec:.3f}")

    print("G10 whatshap minus optimal, error per 10,000 common pairs, by depth (descriptive):")
    for c in CAPS:
        print(f"  cap {c:>2}: " + "; ".join(
            f"{d}x {fmt(tests[('whatshap', 'optimal', d, c)][2])}" for d in DEPTHS))
    print(f"\nscored in {(time.time() - started) / 60:.1f} min")


def decide_hg005(pooled, tests, favours_first, favours_second) -> None:
    """H5-1 and H5-2 of `docs/design/hg005-run-plan.md`; H5-3 is scored by predict_regime."""
    p = pooled[("lib", "full", 15, "whatshap")]
    r = p["objective"] / p["optimum"]
    print(f"H5-1 whatshap / certified optimum, full depth, cap 15: pooled {r:.4f}, worst "
          f"chromosome {p['worst']:.4f} -> {'supported' if r < 0.95 else 'not supported'} "
          "(rule: below 0.95)")
    level = 0.05 / 4
    cap3 = [d for d in DEPTHS if favours_second(("whatshap", "optimal", d, 3), level)]
    cap5 = [d for d in DEPTHS if tests[("whatshap", "optimal", d, 5)][1] < level]
    supported = len(cap3) == len(DEPTHS) and not cap5
    refuted = not any(tests[("whatshap", "optimal", d, 3)][1] < level for d in DEPTHS)
    print(f"H5-2 whatshap vs optimal, level {level}: cap 3 significant for optimal at "
          f"{', '.join(cap3) or 'none'}; cap 5 significant at {', '.join(cap5) or 'none'} -> "
          + ("supported" if supported else "refuted" if refuted else "not supported, not refuted")
          + " (rule: both cap-3 tests for optimal, neither cap-5 test significant)")
    for name in ("two_level", "connected"):
        print(f"descriptive {name} vs optimal, cap 3, block N50 ratio: " + ", ".join(
            f"{d}x {pooled[('lib', d, 3, name)]['n50_kb'] / pooled[('lib', d, 3, 'optimal')]['n50_kb']:.2f}"
            for d in DEPTHS))
    print("regime (descriptive), whatshap minus optimal per 10,000 common pairs:")
    for c in CAPS:
        print(f"  cap {c:>2}: " + "; ".join(
            f"{d}x {fmt(tests[('whatshap', 'optimal', d, c)][2])}" for d in DEPTHS))


def decide_ont(pooled, tests, favours_first, favours_second) -> None:
    """O1 to O7 of `docs/design/ont-run-plan.md`."""
    p = pooled[("lib", "full", 15, "whatshap")]
    o1 = p["objective"] / p["optimum"]
    print(f"O1 whatshap / certified optimum, full depth, cap 15: pooled {o1:.4f}, worst "
          f"chromosome {p['worst']:.4f} -> {'supported' if o1 < 0.95 else 'not supported'} "
          "(rule: below 0.95)")

    level = 0.05 / 12
    cells = [(d, c) for c in (3, 5, 10) for d in DEPTHS]
    for_optimal = [k for k in cells if favours_second(("whatshap", "optimal", *k), level)]
    significant = [k for k in cells if tests[("whatshap", "optimal", *k)][1] < level]
    supported = all((d, 5) in for_optimal for d in DEPTHS)
    refuted = not any(c in (5, 10) for _, c in significant)
    print(f"O2 whatshap vs optimal, caps 3, 5, 10, level {level:.5f}: significant for optimal at "
          + (", ".join(f"{d}x c{c}" for d, c in for_optimal) or "none")
          + "; any direction at " + (", ".join(f"{d}x c{c}" for d, c in significant) or "none")
          + " -> " + ("supported" if supported else "refuted" if refuted else
                      "not supported, not refuted")
          + " (rule: cap 5 for optimal at every depth; refuted if nothing at caps 5 or 10)")

    level = 0.05 / 8
    merge = [(d, c) for d in DEPTHS for c in (3, 5)
             if favours_first(("merge_optimal_1e3", "optimal", d, c), level)]
    against = [(d, c) for d in DEPTHS for c in (3, 5)
               if favours_second(("merge_optimal_1e3", "optimal", d, c), level)]
    print("O3 merge_optimal_1e3 vs optimal, caps 3 and 5, level 0.00625: for the merge at "
          + (", ".join(f"{d}x c{c}" for d, c in merge) or "none") + "; for optimal at "
          + (", ".join(f"{d}x c{c}" for d, c in against) or "none") + " -> "
          + ("supported" if all((d, 5) in merge for d in DEPTHS) else "not supported")
          + " (rule: cap 5 for the merge at every depth)")

    ok = True
    parts = []
    for d in DEPTHS:
        t, o = pooled[("lib", d, 3, "two_level")], pooled[("lib", d, 3, "optimal")]
        gain = t["n50_kb"] / o["n50_kb"]
        extra_pairs = t["pairs"][PRIMARY] - o["pairs"][PRIMARY]
        extra = t["switches"][PRIMARY] - o["switches"][PRIMARY]
        own = o["switches"][PRIMARY] / o["pairs"][PRIMARY]
        joins = extra / extra_pairs if extra_pairs > 0 else float("nan")
        good = gain >= 1.25 and extra_pairs > 0 and joins <= 2 * own
        ok = ok and good
        parts.append(f"{d}x N50 x{gain:.2f}, join rate {extra}/{extra_pairs} = {joins:.6f} "
                     f"against 2 x {own:.6f}")
    print("O4 two_level vs optimal, cap 3: " + "; ".join(parts) + " -> "
          + ("supported" if ok else "not supported")
          + " (rule: N50 up 25 per cent and join rate at most twice optimal's, every depth)")

    level = 0.05 / 4
    worse = [d for d in DEPTHS if favours_second(("error_optimal_1e2", "optimal", d, 3), level)]
    better = [d for d in DEPTHS if favours_first(("error_optimal_1e2", "optimal", d, 3), level)]
    print(f"O5 error_optimal_1e2 vs optimal, cap 3, level {level}: for optimal at "
          f"{', '.join(worse) or 'none'}; for error_optimal at {', '.join(better) or 'none'} -> "
          + ("replicated (worse, as on HiFi)" if len(worse) >= 2 else
             "reversed" if len(better) >= 2 else "neither")
          + " (rule: two or more depths)")

    print("O6 capped at cap 15 against uncapped tools (descriptive): table 1, rows cap 15 and '-'")
    print("O7 whatshap minus optimal, error per 10,000 common pairs, by depth (descriptive; the "
          "HiFi rows are in the v5 re-score of the HiFi run):")
    for c in CAPS:
        print(f"  cap {c:>2}: " + "; ".join(
            f"{d}x {fmt(tests[('whatshap', 'optimal', d, c)][2])}" for d in DEPTHS))


if __name__ == "__main__":
    main()
