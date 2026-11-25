"""The registered tool-level test (`docs/design/tool-run-plan.md`).

Per chromosome, smallest first: stream the HiFi slice, run stock and
exact-selection `whatshap phase` at caps 3, 5 and 15, HiPhase and LongPhase,
four at a time under memory ceilings, then delete the slice. A run is done
when its output exists and its marker is written, so the script resumes.
`--score` scores the outputs against GIAB v5.0q and decides W1 and W2.
With `PHASEBUDGET_PLATFORM=ont` it runs `docs/design/tool-run-ont-plan.md`:
stock and exact WhatsHap only, with `--ignore-read-groups`, deciding WO1
to WO3.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/tool_run.py`
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402

ONT = common.PLATFORM == "ont"
TOOLS_DIR = common.ROOT / "data" / ("tools_ont" if ONT else "tools")
# The plan requires this much free disk before a download (ONT chromosome 1 is about 37 GB).
FREE_GB = 45 if ONT else 20
CAPS = (3, 5, 15)
BIN = Path.home() / ".local" / "bin"
PY = str(common.ROOT / ".venv" / "bin" / "python")


def runs(chrom, bam):
    d = TOOLS_DIR / chrom
    ref = str(common.REFERENCE / f"{chrom}.fa")
    vcf = str(common.chrom_dir(chrom) / "hetsnv.vcf.gz")
    plain = str(common.chrom_dir(chrom) / "hetsnv.vcf")
    out = {}
    groups = ["--ignore-read-groups"] if ONT else []
    for cap in CAPS:
        out[f"stock_c{cap}"] = [str(common.ROOT / ".venv" / "bin" / "whatshap"), "phase",
                                "--reference", ref, "--internal-downsampling", str(cap), *groups,
                                "-o", str(d / f"stock_c{cap}.vcf"), vcf, bam]
        out[f"exact_c{cap}"] = [PY, str(HERE / "whatshap_exact.py"), "phase", "--reference", ref,
                                "--internal-downsampling", str(cap), *groups,
                                "-o", str(d / f"exact_c{cap}.vcf"), vcf, bam]
    if ONT:
        return out
    out["hiphase"] = [str(BIN / "hiphase"), "--bam", bam, "--vcf", vcf, "--output-vcf",
                      str(d / "hiphase.vcf.gz"), "--reference", ref, "--ignore-read-groups",
                      "--threads", "1"]
    out["longphase"] = [str(BIN / "longphase"), "phase", "-s", plain, "-b", bam, "-r", ref,
                        "--pb", "-t", "1", "-o", str(d / "longphase")]
    return out


def run_chromosome(chrom, jobs):
    d = TOOLS_DIR / chrom
    d.mkdir(parents=True, exist_ok=True)
    bam = d / "full.bam"
    todo = {name: cmd for name, cmd in runs(chrom, str(bam)).items()
            if not (d / f"{name}.ok").exists()}
    if not todo:
        return
    if not bam.exists():
        import time

        # The plan requires free disk before a download; wait for it rather than fail.
        while shutil.disk_usage(d).free < FREE_GB * 1e9:
            time.sleep(300)
        import pysam

        here = os.getcwd()
        os.chdir(d)  # htslib caches the remote index in the working directory
        try:
            part = d / "full.part.bam"
            pysam.view("-b", "--output-fmt-option", "level=1", "-o", str(part),
                       common.BAM_URL, chrom, catch_stdout=False)
            part.rename(bam)
            pysam.index(str(bam))
        finally:
            os.chdir(here)
    env = dict(os.environ, PYTHONPATH=f"{common.ROOT / 'src'}")

    def one(item):
        name, cmd = item
        with open(d / f"{name}.log", "w") as log:
            code = subprocess.run(["systemd-run", "--user", "--scope", "-q", "-p", "MemoryMax=6G",
                                   "-p", "MemorySwapMax=0", "nice", "-n", "19",
                                   "/usr/bin/time", "-v", *cmd],
                                  stdout=log, stderr=subprocess.STDOUT, env=env).returncode
        if code == 0:
            (d / f"{name}.ok").write_text("done\n")
        return name, code

    with ThreadPoolExecutor(jobs) as pool:
        results = list(pool.map(one, todo.items()))
    failed = [n for n, code in results if code]
    print(f"{chrom}: {len(results) - len(failed)} runs done, failed: {failed or 'none'}",
          flush=True)
    if not failed:
        for path in list(d.glob("*.bam")) + list(d.glob("*.bai")):
            path.unlink()


def score(chromosomes):
    from phasebudget.bootstrap import statistic_interval
    from phasebudget.metrics import block_lengths, comparable_pairs, n50, switch_positions
    from phasebudget.paired import common_pairs, sign_test, window_discordance

    sc = common._module("score_hg002.py")
    files = {f"{k}_c{c}": f"{k}_c{c}.vcf" for k in ("stock", "exact") for c in CAPS}
    if not ONT:
        files.update({"hiphase": "hiphase.vcf.gz", "longphase": "longphase.vcf"})
    pooled = {n: [0, 0, [], 0.0] for n in files}
    paired = {c: [0, 0, 0, []] for c in CAPS}
    for chrom in chromosomes:
        variants = sc.variant_table(common.hetsnv_path(chrom))
        truth = common.v5_truth(chrom, variants)
        positions = [p - 1 for p, _, _ in variants]
        def existing(path):
            # Finished outputs are gzipped to save disk.
            gz = path.with_name(path.name + ".gz")
            return path if path.exists() or not gz.exists() else gz

        ph = {n: sc.whatshap_phasing(existing(TOOLS_DIR / chrom / f), variants)
              for n, f in files.items()}
        for n, p in ph.items():
            pooled[n][0] += len(switch_positions(truth, p))
            pooled[n][1] += comparable_pairs(p, truth)
            pooled[n][2].extend(block_lengths(p, positions))
            pooled[n][3] += sc.time_log(TOOLS_DIR / chrom / f"{n}.log")[0]
        for c in CAPS:
            rows = common_pairs(truth, ph[f"stock_c{c}"], ph[f"exact_c{c}"])
            q = paired[c]
            q[0] += len(rows)
            q[1] += sum(1 for _, x, y in rows if x and not y)
            q[2] += sum(1 for _, x, y in rows if y and not x)
            q[3].extend(window_discordance(rows, positions, 1_000_000))
    plan = "tool-run-ont-plan.md" if ONT else "tool-run-plan.md"
    print(f"tool-level test, HG002 {'ONT' if ONT else 'HiFi'} full depth, truth GIAB v5.0q, "
          f"{len(chromosomes)} chromosomes (docs/design/{plan})")
    print(f"{'run':>12} {'switches':>9} {'pairs':>9} {'per 10k':>8} {'N50 kb':>7} {'hours':>6}")
    for n, (s, p, blocks, secs) in pooled.items():
        print(f"{n:>12} {s:>9} {p:>9} {1e4 * s / max(p, 1):>8.3f} {n50(blocks) / 1e3:>7.0f} "
              f"{secs / 3600:>6.2f}")
    for c in CAPS:
        pairs, a, b, windows = paired[c]

        def stat(picks, windows=windows):
            n = sum(windows[k][0] for k in picks)
            return None if n == 0 else 1e4 * (sum(windows[k][1] for k in picks)
                                              - sum(windows[k][2] for k in picks)) / n

        interval = statistic_interval(stat, len(windows), resamples=2000, seed=0)
        pv = sign_test(a, a + b)
        if ONT:
            label = {3: "WO1", 5: "WO2", 15: "WO3"}[c]
            directional = c in (3, 5)
        else:
            label, directional = ("W1" if c == 3 else "W2"), c == 3
        verdict = ("supported" if (pv < 0.05 / 3 and b < a) else "not supported") if directional \
            else ("significant" if pv < 0.05 / 3 else "not significant")
        print(f"{label} cap {c}: only stock wrong {a}, only exact wrong {b}, of {pairs} common "
              f"pairs, p {pv:.3g}; stock minus exact per 10k {interval.estimate:+.2f} "
              f"[{interval.low:+.2f}, {interval.high:+.2f}] -> {verdict} (level 0.0167)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chromosomes", nargs="+", default=list(common.CONFIRMATORY))
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--score", action="store_true")
    args = parser.parse_args()
    if common.PLATFORM not in ("hifi", "ont"):
        raise SystemExit("the registered tool-level tests are HiFi and ONT only")
    if args.score:
        score(args.chromosomes)
        return
    order = sorted(args.chromosomes, key=lambda c: (TOOLS_DIR / c).exists() and 0 or 1)
    import genome_run

    for chrom in sorted(order, key=lambda c: genome_run.LENGTH_MB[c]):
        run_chromosome(chrom, args.jobs)


if __name__ == "__main__":
    main()
