"""Score the registered HG002 pipelines against the truth and decide H1 to H6.

Reads the run files `hg002_pipelines.py` wrote, the outputs of `whatshap
phase` and HapCUT2, the variant table, and the two truths of
`docs/design/real-data-plan.md`, and prints:

  1. per depth, pipeline and cap: switch error (switches and comparable
     pairs), flips, phased fraction, block N50 and QAN50 in kb, the additive
     objective over the certified optimum, full-matrix wMEC, runtime and peak
     memory;
  2. the hypotheses with their registered decision rules, intervals from the
     1 Mb window bootstrap (2,000 resamples, seed 0);
  3. the secondary table against the hifiasm phase transfer.

The primary truth is transmission from the parents; variants it leaves
unphased are skipped by every metric (`docs/design/metrics.md`).
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import importlib.util
import json
import math
import re
from pathlib import Path
from typing import Optional

from phasebudget._types import Phasing
from phasebudget.bootstrap import ratio_interval, statistic_interval, window_counts
from phasebudget.metrics import comparable_pairs, summarize, switch_positions

WINDOW = 1_000_000
RESAMPLES = 2000
SEED = 0
LIBRARY = ("whatshap", "optimal", "merge_optimal", "merge_optimal_1e3", "loop", "merge_loop",
           "merge_loop_1e3", "shuffled",
           # added after the chromosome-20 registration, run post hoc on 27 September:
           "whatshap_two_level", "two_level", "merge_two_level_1e3", "connected",
           "loop_oriented", "merge_loop_oriented_1e3")
CAPS = (3, 5, 10, 15)


def _extractor():
    path = Path(__file__).with_name("extract_fragments.py")
    spec = importlib.util.spec_from_file_location("extract_fragments", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def variant_table(vcf):
    """`[(pos1, ref, alt)]` in file order, 1-based."""
    rows = []
    with _open(vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.split("\t", 5)
            rows.append((int(f[1]), f[3], f[4]))
    return rows


def parental_truth(vcf, variants, *, phase_set="PATMAT"):
    """Paternal allele per variant from a `pat|mat` VCF, None where not phased so."""
    index = {v: j for j, v in enumerate(variants)}
    truth: list[Optional[int]] = [None] * len(variants)
    with _open(vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            key = (int(f[1]), f[3], f[4])
            if key not in index:
                continue
            fields = dict(zip(f[8].split(":"), f[9].split(":")))
            gt = fields.get("GT", "")
            if gt not in ("0|1", "1|0") or fields.get("PS") != phase_set:
                continue
            truth[index[key]] = int(gt[0])
    return truth


def whatshap_phasing(vcf, variants):
    """Alleles and blocks from a `whatshap phase` VCF; unphased calls are None."""
    index = {v: j for j, v in enumerate(variants)}
    alleles: list[Optional[int]] = [None] * len(variants)
    blocks: list[Optional[int]] = [None] * len(variants)
    with _open(vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            key = (int(f[1]), f[3], f[4])
            if key not in index:
                continue
            fields = dict(zip(f[8].split(":"), f[9].split(":")))
            gt = fields.get("GT", "")
            if gt in ("0|1", "1|0") and "PS" in fields and fields["PS"] != ".":
                alleles[index[key]] = int(gt[0])
                blocks[index[key]] = int(fields["PS"])
    return Phasing(tuple(alleles), tuple(blocks))


def hapcut2_phasing(path, variants):
    """Alleles and blocks from a HapCUT2 block file; `-` alleles are unphased."""
    index = {v: j for j, v in enumerate(variants)}
    alleles: list[Optional[int]] = [None] * len(variants)
    blocks: list[Optional[int]] = [None] * len(variants)
    block = -1
    with open(path) as handle:
        for line in handle:
            if line.startswith("BLOCK"):
                block += 1
                continue
            if line.startswith("*") or not line.strip():
                continue
            f = line.split("\t")
            if f[1] == "-" or f[1] not in ("0", "1"):
                continue
            key = (int(f[4]), f[5], f[6])
            if key in index:
                alleles[index[key]] = int(f[1])
                blocks[index[key]] = block
    return Phasing(tuple(alleles), tuple(blocks))


def run_phasing(record):
    return Phasing(tuple(None if a < 0 else a for a in record["alleles"]),
                   tuple(None if b < 0 else b for b in record["blocks"]))


def time_log(path):
    """`(seconds, peak MB)` summed and maximized over the GNU time blocks in a log."""
    seconds, peak = 0.0, 0.0
    text = Path(path).read_text() if Path(path).exists() else ""
    for match in re.finditer(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)", text):
        parts = [float(p) for p in match.group(1).split(":")]
        seconds += sum(p * 60 ** i for i, p in enumerate(reversed(parts)))
    for match in re.finditer(r"Maximum resident set size \(kbytes\): (\d+)", text):
        peak = max(peak, int(match.group(1)) / 1024)
    return seconds, peak


def auc(bad, good):
    """P(a value at a wrong call < a value at a right one), ties one half, by sorting."""
    if not bad or not good:
        return float("nan")
    ordered = sorted(good)
    wins = 0.0
    for value in bad:
        below = bisect.bisect_left(ordered, value)
        equal = bisect.bisect_right(ordered, value) - below
        wins += (len(ordered) - below - equal) + 0.5 * equal
    return wins / (len(bad) * len(good))


def describe(truth, phasing, positions):
    report = summarize(truth, phasing, positions)
    return {
        "switch": report["switch_error_rate"],
        "switches": len(switch_positions(truth, phasing)),
        "pairs": comparable_pairs(phasing, truth),
        "flips": int(report["flip_errors"]),
        "phased": report["phased_fraction"],
        "n50": report["block_n50"] / 1000,
        "qan50": report["qan50"] / 1000,
    }


def fmt_interval(interval, digits=3):
    return (f"{interval.estimate:.{digits}f} [{interval.low:.{digits}f}, "
            f"{interval.high:.{digits}f}]"
            + (f", undefined in {interval.undefined}" if interval.undefined else ""))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="data/hg002")
    parser.add_argument("--giab", default="data/giab")
    parser.add_argument("--exclude", default="",
                        help="post hoc: variant index ranges a-b (inclusive, comma separated) "
                             "whose truth is set to unphased; the registered analysis uses none")
    args = parser.parse_args()
    data, giab = Path(args.data), Path(args.giab)
    extractor = _extractor()

    variants = variant_table(data / "HG002_chr20_hetsnv_bed.vcf")
    matrices, positions = {}, None
    for depth, name in (("full", "HG002_chr20_fragments.json.gz"),
                        ("half", "HG002_chr20_half_fragments.json.gz")):
        matrix, pos, _, _ = extractor.load_fragments(str(data / name))
        matrices[depth] = extractor.phred_to_log_odds(matrix)
        positions = pos if positions is None else positions
        assert pos == positions
    assert [p - 1 for p, _, _ in variants] == positions, "variant table and fragments disagree"
    truth = parental_truth(giab / "HG002_chr20_trio_truth.vcf.gz", variants)
    secondary = parental_truth(giab / "HG002_chr20_hifiasm_phasetransfer.vcf.gz", variants)
    excluded = 0
    for part in filter(None, args.exclude.split(",")):
        low, high = (int(x) for x in part.split("-"))
        for j in range(low, high + 1):
            excluded += truth[j] is not None
            truth[j] = None
    if args.exclude:
        print(f"POST HOC: truth set to unphased at variants {args.exclude} ({excluded} phased "
              "sites); not the registered analysis")
    both = [(a, b) for a, b in zip(truth, secondary) if a is not None and b is not None]
    print(f"variants {len(variants)}; primary truth phases {sum(t is not None for t in truth)}; "
          f"hifiasm PATMAT phases {sum(t is not None for t in secondary)}, of which "
          f"{sum(t is not None and p is None for t, p in zip(secondary, truth))} are unphased "
          f"by the primary; discordant where both phase: "
          f"{sum(a != b for a, b in both)} of {len(both)}")
    print(f"bootstrap: {WINDOW // 1_000_000} Mb windows, {RESAMPLES} resamples, seed {SEED}")

    runs = {}
    for path in sorted((data / "runs").glob("*.json.gz")):
        with gzip.open(path, "rt") as handle:
            record = json.load(handle)
        depth = path.name.split("_", 1)[0]
        runs[(depth, record["pipeline"], record["cap"], record["weights"])] = record
    commits = sorted({r["commit"] for r in runs.values()})
    print(f"run files {len(runs)}, pipeline commits {', '.join(commits)}")

    phasings = {}
    for key, record in runs.items():
        phasings[key] = run_phasing(record)
    external = {}
    for depth in ("full", "half"):
        for cap in CAPS:
            vcf = data / "external" / f"whatshap_{depth}_c{cap}.vcf"
            if vcf.exists():
                external[(depth, "whatshap phase", cap)] = (
                    whatshap_phasing(vcf, variants),
                    time_log(data / "external" / f"whatshap_{depth}_c{cap}.log"))
        blocks = data / "external" / f"hapcut2_{depth}.blocks"
        if blocks.exists():
            seconds, peak = time_log(data / "external" / f"hapcut2_{depth}.extract.log")
            s2, p2 = time_log(data / "external" / f"hapcut2_{depth}.log")
            external[(depth, "HapCUT2", None)] = (hapcut2_phasing(blocks, variants),
                                                   (seconds + s2, max(peak, p2)))

    # 1. The main table.
    header = (f"{'depth':>5} {'cap':>4} {'pipeline':>18} {'weights':>7} {'switch':>8} "
              f"{'sw/pairs':>12} {'flips':>5} {'phased':>6} {'N50 kb':>7} {'QAN50 kb':>8} "
              f"{'obj/opt':>7} {'MEC full':>9} {'rounds':>6} {'sec':>7} {'peak MB':>7}")
    print()
    print("1. every pipeline against the transmission truth; MEC full is under log-odds "
          "weights, n/a for external tools; obj/opt only where the reads are not merged")
    print(header)
    print("-" * len(header))
    counts = {}
    for depth in ("full", "half"):
        for cap in CAPS:
            for pipeline in LIBRARY:
                for weights in ("logodds", "phred"):
                    key = (depth, pipeline, cap, weights)
                    if key not in runs:
                        continue
                    record = runs[key]
                    m = describe(truth, phasings[key], positions)
                    counts[key] = window_counts(truth, phasings[key], positions, WINDOW)
                    mec = matrices[depth].cost(tuple(0 if a < 0 else a for a in record["alleles"]))
                    ratio = ("n/a" if record["objective"] is None
                             else f"{record['objective'] / record['optimum']:.4f}")
                    print(f"{depth:>5} {cap:>4} {pipeline:>18} {weights:>7} {m['switch']:>8.5f} "
                          f"{m['switches']:>5}/{m['pairs']:<6} {m['flips']:>5} "
                          f"{m['phased']:>6.3f} {m['n50']:>7.0f} {m['qan50']:>8.0f} "
                          f"{ratio:>7} {mec:>9.0f} {record['rounds']:>6} "
                          f"{record['seconds']:>7.1f} {record['peak_mb']:>7.0f}")
            for (d, name, c), (phasing, (seconds, peak)) in sorted(
                    external.items(), key=lambda x: (x[0][0], x[0][1], x[0][2] or 0)):
                if d != depth or (c is not None and c != cap) or (c is None and cap != 15):
                    continue
                m = describe(truth, phasing, positions)
                counts[(d, name, c, "tool")] = window_counts(truth, phasing, positions, WINDOW)
                label = name if c is not None else f"{name} (no cap)"
                note = "" if (c == 15 or c is None) else " *"
                print(f"{depth:>5} {cap if c else '-':>4} {label + note:>18} {'tool':>7} "
                      f"{m['switch']:>8.5f} {m['switches']:>5}/{m['pairs']:<6} {m['flips']:>5} "
                      f"{m['phased']:>6.3f} {m['n50']:>7.0f} {m['qan50']:>8.0f} {'n/a':>7} "
                      f"{'n/a':>9} {'':>6} {seconds:>7.1f} {peak:>7.0f}")
            print()
    print("* not registered: whatshap phase with --internal-downsampling other than its default of 15 "
          "(WhatsHap 2.8 ignores --max-coverage)")

    # 2. The hypotheses.
    def ci(numerator, denominator):
        """The ratio's interval, or None when the ratio is undefined on the full data."""
        try:
            return ratio_interval(counts[numerator], counts[denominator], resamples=RESAMPLES,
                                  seed=SEED)
        except ValueError:
            return None

    def show(interval, digits=3):
        if interval is None:
            return "undefined: the denominator pipeline has no switch error"
        return fmt_interval(interval, digits)

    def above(interval):
        return interval is not None and interval.low > 1

    def contains(interval):
        return interval is not None and interval.low <= 1 <= interval.high

    print()
    print("2. hypotheses, full depth unless stated, log-odds weights (amendment 3); the "
          "registered unscaled weights follow where they differ")
    for depth in ("full", "half"):
        print(f"-- depth {depth}")
        key = (depth, "whatshap", 15, "logodds")
        if key in runs:
            r = runs[key]["objective"] / runs[key]["optimum"]
            print(f"H1 whatshap / optimum objective at cap 15: {r:.4f} -> "
                  f"{'supported' if r < 0.95 else 'not supported'} (rule: below 0.95)")
        h2 = {}
        for cap in (3, 5, 15):
            a, b = (depth, "whatshap", cap, "logodds"), (depth, "optimal", cap, "logodds")
            if a in counts and b in counts:
                h2[cap] = ci(a, b)
                print(f"H2 cap {cap:>2} switch whatshap / optimal: {show(h2[cap])}")
        if len(h2) == 3:
            small = all(above(h2[c]) for c in (3, 5))
            large = contains(h2[15])
            refuted = all(contains(h2[c]) for c in (3, 5))
            print("H2 decision: " + ("supported" if small and large else
                                     "refuted" if refuted else "not supported, not refuted")
                  + f" (caps 3 and 5 above 1: {small}; cap 15 contains 1: {large})")
        for weights in ("logodds", "phred"):
            h3 = {}
            for cap in (3, 5):
                a, b = (depth, "optimal", cap, "logodds"), (depth, "merge_loop", cap, weights)
                if a in counts and b in counts:
                    h3[cap] = ci(a, b)
                    print(f"H3 [{weights}] cap {cap} switch optimal / merge_loop: "
                          f"{show(h3[cap])}")
            if len(h3) == 2:
                supported = all(above(i) for i in h3.values())
                refuted = all(contains(i) for i in h3.values())
                print(f"H3 [{weights}] decision: " + ("supported" if supported else
                      "refuted" if refuted else "not supported, not refuted"))
            for cap in (3, 5):
                keys = [(depth, "optimal", cap, "logodds"), (depth, "loop", cap, weights),
                        (depth, "shuffled", cap, weights)]
                if not all(k in counts for k in keys):
                    continue
                o, l, s = (counts[k] for k in keys)

                def kept(picks, o=o, l=l, s=s):
                    rate = lambda c: (sum(c[k][0] for k in picks) / sum(c[k][1] for k in picks)
                                      if sum(c[k][1] for k in picks) else None)
                    ro, rl, rs = rate(o), rate(l), rate(s)
                    if None in (ro, rl, rs) or ro == rl:
                        return None
                    return (ro - rs) / (ro - rl)

                try:
                    interval = statistic_interval(kept, len(o), resamples=RESAMPLES, seed=SEED)
                except ValueError:
                    print(f"H4 [{weights}] cap {cap}: undefined, the loop and static selection "
                          "have the same switch error")
                    continue
                verdict = ("kept, as on simulation" if interval.estimate >= 0.8 else
                           "OQ-13 reopened" if interval.estimate < 0.5 else "in between")
                everything = range(len(o))
                rate = lambda c: sum(c[k][0] for k in everything) / sum(c[k][1] for k in everything)
                if sum(c[0] for c in l) > sum(c[0] for c in o):
                    verdict += ("; but the loop raised switch error, from "
                                f"{rate(o):.5f} to {rate(l):.5f}, so there was no reduction "
                                "to keep and the rule's premise fails")
                print(f"H4 [{weights}] cap {cap} (optimal - shuffled) / (optimal - loop): "
                      f"{fmt_interval(interval)} -> {verdict} (rule: at least 0.8; below 0.5 "
                      "at both caps reopens OQ-13)")
        for weights in ("logodds", "phred"):
            key = (depth, "optimal", 5, weights)
            if key not in runs or "signals" not in runs[key]:
                continue
            rows = []
            for b, called, predicted, post_margin, mec_margin in runs[key]["signals"]:
                if truth[b] is None or truth[b + 1] is None:
                    continue
                wrong = called != (truth[b] ^ truth[b + 1])
                rows.append((b, wrong, predicted, post_margin, mec_margin))
            n_windows = (positions[-1] - positions[0]) // WINDOW + 1
            per_window = [[0.0, 0.0] for _ in range(n_windows)]
            for b, wrong, predicted, _, _ in rows:
                k = (positions[b + 1] - positions[0]) // WINDOW
                per_window[k][0] += wrong
                per_window[k][1] += predicted
            observed = sum(r[1] for r in rows)
            predicted = sum(r[2] for r in rows)

            def calibration(picks):
                p = sum(per_window[k][1] for k in picks)
                return None if p == 0 else sum(per_window[k][0] for k in picks) / p

            try:
                interval = statistic_interval(calibration, n_windows, resamples=RESAMPLES,
                                              seed=SEED)
            except ValueError:
                interval = None
            print(f"H5 [{weights}] cap 5 optimal, {len(rows)} truth-phased spanned boundaries: "
                  f"observed wrong {observed}, predicted {predicted:.2f}; observed / predicted "
                  f"{show(interval, 2)} -> "
                  + ("supported (under-prediction)" if above(interval) else "not supported"))
            bad = [r[3] for r in rows if r[1]]
            good = [r[3] for r in rows if not r[1]]
            value = auc(bad, good)
            mec = auc([r[4] for r in rows if r[1]], [r[4] for r in rows if not r[1]])
            print(f"H6 [{weights}] cap 5 AUC of the posterior margin for a wrong call: "
                  f"{value:.3f} -> {'supported' if value > 0.6 else 'not supported'} "
                  f"(rule: above 0.6); MEC margin {mec:.3f}")

    # 4 is printed after 3 but computed here, while the phasings are at hand.
    everything = [phasings[k] for k in runs if k[3] == "logodds"] + \
        [phasing for phasing, _ in external.values()]
    recurring = recurring_switches(truth, secondary, positions, everything,
                                   max(2, len(everything) // 2))

    # 3. Secondary truth.
    print()
    print("3. against the hifiasm phase transfer (PATMAT sites; not a hypothesis test)")
    for key in sorted(k for k in runs if k[3] == "logodds"):
        m = describe(secondary, phasings[key], positions)
        print(f"{key[0]:>5} {key[2]:>4} {key[1]:>18} {m['switch']:>8.5f} "
              f"{m['switches']:>5}/{m['pairs']:<6}")
    for (d, name, c), (phasing, _) in sorted(external.items(),
                                             key=lambda x: (x[0][0], x[0][1], x[0][2] or 0)):
        m = describe(secondary, phasing, positions)
        print(f"{d:>5} {c if c else '-':>4} {name:>18} {m['switch']:>8.5f} "
              f"{m['switches']:>5}/{m['pairs']:<6}")

    print()
    print(f"4. switch sites shared by at least half of the {len(everything)} phasings, and "
          "whether hifiasm gives the same phase relation as the transmission truth between the "
          "site and the previous site both phase")
    for q, position, count, relation in recurring:
        verdict = {True: "same relation", False: "OPPOSITE relation", None: "not comparable"}
        print(f"variant {q:>6} at {position:>9}: {count:>3} phasings; {verdict[relation]}")



def recurring_switches(truth, secondary, positions, phasings, minimum):
    """Switch sites shared by at least `minimum` phasings, with both truths there.

    A site where nearly every method switches is more likely an error of the
    truth than of every method at once; the hifiasm call at the site is the
    independent check.
    """
    hits: dict[int, int] = {}
    for phasing in phasings:
        for q in switch_positions(truth, phasing):
            hits[q] = hits.get(q, 0) + 1
    rows = []
    for q, count in sorted(hits.items(), key=lambda x: (-x[1], x[0])):
        if count < minimum:
            break
        previous = next((p for p in range(q - 1, -1, -1)
                         if truth[p] is not None and secondary[p] is not None), None)
        if previous is None or secondary[q] is None:
            relation = None
        else:
            relation = (truth[previous] ^ truth[q]) == (secondary[previous] ^ secondary[q])
        rows.append((q, positions[q] + 1, count, relation))
    return rows


if __name__ == "__main__":
    main()
