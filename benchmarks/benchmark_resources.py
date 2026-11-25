"""Runtime and peak memory of selection and of the dynamic program.

The proposal's milestone 2 asks for time and memory against the coverage cap
and the chromosome length. Each configuration runs in its own Python process,
so one measurement cannot inherit another's heap. Peak memory is the kernel's
high-water mark of resident memory (`VmHWM`), reset through
`/proc/self/clear_refs` just before the measured step and reported as the rise
above the resident size at the reset. That counts the native pass's
allocations, which `tracemalloc` cannot see. It needs Linux.

Instances come from the genomic simulator with the `ont` profile of
`benchmark_selection_at_scale.py`. The dynamic program runs on the reads that
`select_optimal` keeps at the cap, so its active coverage is the cap. Four
steps are timed:

  * `select_optimal` and `select_whatshap` on the full read set;
  * the native pass, with checkpoints every `ceil(sqrt(n))` columns (the
    default) and with none (one segment, so every column is stored), each for
    the optimum alone and with the backward pass that margins need;
  * the Python pass, for the optimum alone, where it finishes in reasonable
    time.

`prep s` is the Python time to lay out the columns for the native pass and is
included in `total s`; `pass s` is the native call alone.
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import pickle
import platform
import subprocess
import sys
import tempfile
import time
from pathlib import Path

PROFILE = {"read_length": 10_000, "error_rate": 0.05}


def _status(field: str) -> int:
    """A `/proc/self/status` field in kB."""
    with open("/proc/self/status") as handle:
        for line in handle:
            if line.startswith(field + ":"):
                return int(line.split()[1])
    raise RuntimeError(f"{field} not in /proc/self/status")


def _reset_peak() -> int:
    with open("/proc/self/clear_refs", "w") as handle:
        handle.write("5")
    return _status("VmRSS")


def prepare(n_variants: int, depth: float, cap: int, seed: int, path: str) -> None:
    from phasebudget.selection import read_weights, select_optimal
    from phasebudget.simulate import simulate_genomic

    instance = simulate_genomic(n_variants, coverage=depth, heterozygosity=0.001,
                                quality_spread=2.0, seed=seed, **PROFILE)
    matrix = instance.matrix
    selection = select_optimal(matrix, cap, read_weights(matrix))
    with open(path, "wb") as handle:
        pickle.dump({"matrix": matrix, "positions": instance.positions,
                     "selected": selection.selected}, handle)


def measure(step: str, path: str, cap: int) -> dict:
    from phasebudget import _native
    from phasebudget.dp import CoverageDP
    from phasebudget.selection import read_weights, select_optimal, select_whatshap

    with open(path, "rb") as handle:
        data = pickle.load(handle)
    matrix = data["matrix"]
    subset = matrix.subset(data["selected"])
    weights = read_weights(matrix)
    record = {"reads": len(matrix), "kept": len(subset),
              "active": subset.max_coverage(), "prep_s": 0.0}

    baseline = _reset_peak()
    started = time.perf_counter()
    if step == "select_optimal":
        select_optimal(matrix, cap, weights)
    elif step == "select_whatshap":
        select_whatshap(matrix, cap, weights, positions=data["positions"])
    elif step == "python_solve":
        CoverageDP(subset, backend="python").solve()
    else:
        margins = step.endswith("margins")
        interval = subset.n_variants if "full" in step else 0
        dp = CoverageDP(subset, backend="native")
        offsets, active, c0, c1 = dp.native_inputs()
        record["prep_s"] = time.perf_counter() - started
        _native.run(subset.n_variants, offsets, active, c0, c1, len(subset),
                    want_margins=margins, interval=interval)
    record["total_s"] = time.perf_counter() - started
    record["peak_mb"] = (_status("VmHWM") - baseline) / 1024
    return record


STEPS = ("select_optimal", "select_whatshap", "native_sqrt_solve", "native_sqrt_margins",
         "native_full_solve", "native_full_margins", "python_solve")


def feasible(step: str, n_variants: int, cap: int, memory_gb: float, python_cap: int,
             python_variants: int) -> bool:
    if step == "python_solve":
        return cap <= python_cap and n_variants <= python_variants
    if "full" in step:
        # Every column stored: 2^cap doubles each, plus the traceback's bytes.
        return (2 ** cap) * n_variants * 9 / 1e9 <= memory_gb
    return True


def run_child(arguments: list[str]) -> dict:
    output = subprocess.run([sys.executable, __file__, *arguments], capture_output=True,
                            text=True, env={**os.environ})
    if output.returncode != 0:
        raise RuntimeError(output.stderr)
    return json.loads(output.stdout.strip().splitlines()[-1])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", type=int, nargs="+", default=[500, 2000, 8000])
    parser.add_argument("--caps", type=int, nargs="+", default=[5, 10, 15, 20])
    parser.add_argument("--depth", type=float, default=30.0)
    parser.add_argument("--seed", type=int, default=8000)
    parser.add_argument("--memory-gb", type=float, default=4.0,
                        help="skip the unchecked native pass when its table would exceed this")
    parser.add_argument("--python-cap", type=int, default=10)
    parser.add_argument("--python-variants", type=int, default=2000)
    parser.add_argument("--worker", nargs=3, metavar=("STEP", "PATH", "CAP"),
                        help=argparse.SUPPRESS)
    parser.add_argument("--prepare", nargs=5, metavar=("N", "DEPTH", "CAP", "SEED", "PATH"),
                        help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.prepare:
        n, depth, cap, seed, path = args.prepare
        prepare(int(n), float(depth), int(cap), int(seed), path)
        print("{}")
        return
    if args.worker:
        step, path, cap = args.worker
        print(json.dumps(measure(step, path, int(cap))))
        return

    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    from phasebudget import _native

    print(f"commit {commit}; {platform.system()} {platform.release()} {platform.machine()}; "
          f"Python {platform.python_version()}; native library "
          f"{'built' if _native.available() else 'absent'}; {os.cpu_count()} logical CPUs")
    print(f"variants {' '.join(map(str, args.variants))}, caps {' '.join(map(str, args.caps))}, "
          f"depth {args.depth:g}, heterozygosity 0.001, quality spread 2, ont profile "
          f"({PROFILE['read_length']} bp, error {PROFILE['error_rate']}), seed {args.seed}; "
          f"one instance per length; simulated data")
    print(f"unchecked native pass skipped above {args.memory_gb:g} GB of table; Python pass "
          f"only at cap <= {args.python_cap} and <= {args.python_variants} variants")
    print("peak MB is the rise of VmHWM over VmRSS at the start of the step, in its own process")
    header = (f"{'variants':>8} {'reads':>6} {'cap':>4} {'active':>6} {'step':>20} "
              f"{'prep s':>8} {'total s':>9} {'peak MB':>9}")
    print(header)
    print("-" * len(header))
    sys.stdout.flush()
    with tempfile.TemporaryDirectory() as scratch:
        for n_variants, cap in itertools.product(args.variants, args.caps):
            path = str(Path(scratch) / f"{n_variants}_{cap}.pickle")
            run_child(["--prepare", str(n_variants), str(args.depth), str(cap),
                       str(args.seed), path])
            for step in STEPS:
                if not feasible(step, n_variants, cap, args.memory_gb, args.python_cap,
                                args.python_variants):
                    print(f"{n_variants:>8} {'':>6} {cap:>4} {'':>6} {step:>20} "
                          f"{'skipped':>8}")
                    continue
                record = run_child(["--worker", step, path, str(cap)])
                print(f"{n_variants:>8} {record['reads']:>6} {cap:>4} {record['active']:>6} "
                      f"{step:>20} {record['prep_s']:>8.2f} {record['total_s']:>9.2f} "
                      f"{record['peak_mb']:>9.1f}")
                sys.stdout.flush()
            print()


if __name__ == "__main__":
    main()
