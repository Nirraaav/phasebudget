"""Run one registered pipeline of the real-data plan on an extracted fragment file.

`docs/design/real-data-plan.md` fixes the pipelines, their parameters and the
caps; this script runs one pipeline at one cap on one fragment file and
writes what `score_hg002.py` needs: the phase call with its blocks, the
selection, the objective against the certified optimum, the full-matrix wMEC,
wall time and peak memory. Nothing here reads a truth.

Weights: `logodds` converts WhatsHap's phred qualities to the library's
natural log-odds (amendment 3); `phred` uses them unscaled, as first
registered. With `--signals` the static pipelines also write, for every
boundary the selected reads span, the called relation, the posterior's
probability that it is wrong, the posterior margin and the MEC margin, for
H5 and H6.

Run once per configuration, in its own process, so peak memory is the
pipeline's own.
"""

from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
import math
import resource
import subprocess
import time
from pathlib import Path

from phasebudget.adaptive import adaptive_selection
from phasebudget.dp import CoverageDP
from phasebudget.hapchat import merge_reads
from phasebudget.metrics import phasing_from_solution
from phasebudget.posterior import boundary_log_sums
from phasebudget.selection import (read_weights, select_connected, select_error_optimal,
                                   select_fischer_marschall,
                                   select_optimal, select_whatshap, selection_objective,
                                   switch_model_costs)
from phasebudget.twolevel import orient_blocks

PIPELINES = ("whatshap", "optimal", "merge_optimal", "merge_optimal_1e3", "loop",
             "merge_loop", "merge_loop_1e3", "shuffled",
             # added after the chromosome-20 registration, 27 September:
             "whatshap_two_level", "two_level", "merge_two_level_1e3", "connected",
             "loop_oriented", "merge_loop_oriented_1e3",
             # D1 pilot (docs/design/error-optimal-selection.md), 27 September:
             "error_optimal_1e2", "error_optimal_1e4", "error_optimal_1e6",
             # Selector fix, 30 September: Fischer and Marschall's published procedure.
             "fm_paper")
# Pooled allele error of the recalibration fit on the odd chromosomes.
ALLELE_ERROR = 0.00117
LOOP = dict(iterations=4, boost=2.0, error_probability=0.05, acceptance="stitched", window=25,
            shuffle_seed=0)


def _extractor():
    path = Path(__file__).with_name("extract_fragments.py")
    spec = importlib.util.spec_from_file_location("extract_fragments", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def whatshap_selection(matrix, cap, positions):
    """WhatsHap 2.8's own compiled `readselection`, when WhatsHap is installed.

    `select_whatshap` reproduces it read for read on every instance checked,
    including all eight HG002 depth and cap pairs, but is a Python port and
    took 20 minutes at cap 15 on chromosome 20. Weights are scaled by 10 and
    rounded, since WhatsHap's qualities are integers; with WhatsHap's own
    constant qualities the order of reads is unchanged.
    """
    from phasebudget.selection import SelectionResult, is_phase_informative

    try:
        import whatshap  # noqa: F401
    except ImportError:
        return select_whatshap(matrix, cap, read_weights(matrix), positions=positions)
    bridge_path = Path(__file__).with_name("whatshap_bridge.py")
    spec = importlib.util.spec_from_file_location("whatshap_bridge", bridge_path)
    bridge = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bridge)
    informative = matrix.subset([f.read_id for f in matrix if is_phase_informative(f)])
    selected = bridge.whatshap_selection(bridge.integer_weights(informative, scale=10.0), cap,
                                         positions)
    weights = read_weights(matrix)
    return SelectionResult(tuple(selected), selection_objective(weights, selected),
                           "whatshap", cap)


def logistic(x: float) -> float:
    return math.exp(-x) / (1.0 + math.exp(-x)) if x >= 0 else 1.0 / (1.0 + math.exp(x))


def signals(subset, solution):
    """Per spanned boundary: `(boundary, called, predicted wrong, posterior margin, MEC margin)`."""
    sums = boundary_log_sums(subset)
    mec = CoverageDP(subset).margins()
    spanned = subset.spanned_boundaries()
    blocks = subset.blocks()
    rows = []
    for b in range(subset.n_variants - 1):
        # A boundary next to a variant no selected read calls carries no evidence.
        if not spanned[b] or blocks[b] is None or blocks[b + 1] is None:
            continue
        called = solution.haplotype[b] ^ solution.haplotype[b + 1]
        log_odds = sums[b][called] - sums[b][1 - called]
        rows.append((b, called, logistic(log_odds), abs(sums[b][0] - sums[b][1]),
                     mec[b].margin))
    return rows


def run(matrix, pipeline, cap, positions, want_signals):
    base = read_weights(matrix)
    optimum = select_optimal(matrix, cap, base)
    record = {"optimum": optimum.objective}
    merged_threshold = 1e3 if pipeline.endswith("_1e3") else None
    merging = pipeline.startswith("merge")
    work = matrix
    if pipeline.startswith("merge"):
        merged = merge_reads(matrix, **({"threshold": merged_threshold}
                                        if merged_threshold else {}))
        work = merged.matrix
        record["merged_reads"] = len(work)
    oriented = None
    if pipeline in ("whatshap", "optimal", "merge_optimal", "merge_optimal_1e3",
                    "whatshap_two_level", "two_level", "merge_two_level_1e3", "connected") \
            or pipeline.startswith("error_optimal") or pipeline == "fm_paper":
        if pipeline.startswith("whatshap"):
            selection = whatshap_selection(work, cap, positions)
        elif pipeline == "connected":
            weights = read_weights(work)
            positive = sorted(w for w in weights.values() if w > 0)
            # Penalty of one median read weight per missing spanning read; on
            # simulation every multiple from 0.5 to 10 gave the same selection.
            selection = select_connected(work, cap, weights, k=1,
                                         penalty=positive[len(positive) // 2] if positive else 0.0)
        elif pipeline == "fm_paper":
            selection = select_fischer_marschall(work, cap, read_weights(work))
        elif pipeline.startswith("error_optimal"):
            scale = float(pipeline.rsplit("_", 1)[1])
            selection = select_error_optimal(work, cap, read_weights(work),
                                             costs=switch_model_costs(cap, ALLELE_ERROR, scale))
        else:
            selection = select_optimal(work, cap, read_weights(work))
        subset = work.subset(selection.selected)
        solution = CoverageDP(subset).solve()
        selected = selection.selected
        support = selected
        record["rounds"] = 1
        if "two_level" in pipeline:
            # Orientation always uses the original reads, merged or not.
            oriented = orient_blocks(matrix, subset, solution)
            record["dropped_reads"] = oriented.dropped_reads
            record["orientation_window"] = oriented.window
        if want_signals:
            record["signals"] = signals(subset, solution)
    else:
        signal = "shuffled" if pipeline == "shuffled" else "posterior"
        loop = dict(LOOP)
        if "oriented" in pipeline:
            loop["acceptance"] = "stitched_oriented"
        result = adaptive_selection(work, cap, signal=signal, **loop)
        selected = result.selection.selected
        support = result.support or selected
        subset = work.subset(support)
        solution = result.solution
        record["rounds"] = result.iterations
        record["stitched"] = bool(result.support) and result.support != selected
    phasing = oriented.phasing if oriented is not None else phasing_from_solution(subset, solution)
    if oriented is not None:
        solution = type(solution)(cost=oriented.cost, haplotype=oriented.haplotype,
                                  assignment=solution.assignment)
    record.update({
        "selected": len(selected),
        "support": len(support),
        "objective": (selection_objective(base, selected) if not merging else None),
        "mec_full": matrix.cost(solution.haplotype),
        "alleles": [-1 if a is None else a for a in phasing.alleles],
        "blocks": [-1 if b is None else b for b in phasing.blocks],
    })
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fragments", required=True)
    parser.add_argument("--pipeline", required=True, choices=PIPELINES)
    parser.add_argument("--cap", type=int, required=True)
    parser.add_argument("--weights", choices=["logodds", "phred", "recalibrated"],
                        default="logodds")
    parser.add_argument("--table", help="recalibration table for --weights recalibrated")
    parser.add_argument("--signals", action="store_true")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    extractor = _extractor()
    matrix, positions, _, provenance = extractor.load_fragments(args.fragments)
    if args.weights == "logodds":
        matrix = extractor.phred_to_log_odds(matrix)
    elif args.weights == "recalibrated":
        spec = importlib.util.spec_from_file_location(
            "recalibrate", Path(__file__).with_name("recalibrate.py"))
        recalibrate = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(recalibrate)
        table = json.loads(Path(args.table).read_text())["bins"]
        matrix = recalibrate.apply_table(matrix, table)
    started = time.time()
    record = run(matrix, args.pipeline, args.cap, positions, args.signals)
    record["seconds"] = round(time.time() - started, 2)
    record["peak_mb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)
    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    record.update({"pipeline": args.pipeline, "cap": args.cap, "weights": args.weights,
                   "fragments": Path(args.fragments).name, "commit": commit,
                   "extraction": provenance})
    with gzip.open(args.out, "wt") as handle:
        json.dump(record, handle)
    print(f"{args.pipeline} cap {args.cap} {args.weights}: {record['seconds']} s, "
          f"{record['peak_mb']} MB, rounds {record['rounds']}, selected {record['selected']}")


if __name__ == "__main__":
    main()
