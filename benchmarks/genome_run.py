"""Run the registered genome-wide HG002 analysis, chromosome by chromosome, safely.

`docs/design/genome-run-plan.md` fixes what runs. This script schedules it
on one machine without endangering it:

* every job runs under `systemd-run --user --scope` with a memory ceiling and
  no swap, at `nice 19` and idle I/O priority, so a misbehaving job is killed
  alone and the desktop stays responsive;
* a job starts only if its estimated memory plus a reserve is available,
  at most `--jobs` run at once, and at most two downloads;
* a download starts only if at least `--min-free-gb` of disk is free and at
  most `--max-bam-chromosomes` chromosomes hold reads on disk;
* a chromosome's reads are deleted as soon as every job that reads them has
  finished; fragments, phasings and logs are kept;
* a job is done only when it exited 0, which writes a marker, so the run
  resumes after an interruption by skipping marked jobs.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/genome_run.py`
with `--chromosomes` to restrict the run; `--task` is used internally.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402

PYTHON = str(Path(sys.executable))
LIBRARY = ("whatshap", "optimal", "merge_optimal", "merge_optimal_1e3", "two_level",
           "whatshap_two_level", "merge_two_level_1e3", "connected")
LOOPS = ("loop", "loop_oriented", "merge_loop_1e3", "shuffled")
TOOLS = ("whatshap", "hapcut2")
if common.PLATFORM == "ont":
    LIBRARY = ("whatshap", "optimal", "merge_optimal_1e3", "two_level", "connected",
               "error_optimal_1e2")
    LOOPS = ()
    TOOLS = ("whatshap", "hapcut2", "longphase")
elif common.PLATFORM == "hg005":
    LIBRARY = ("whatshap", "optimal", "two_level", "merge_optimal_1e3", "connected")
    LOOPS = ()
    TOOLS = ("whatshap",)
if common.SEEDS:
    # Extra subsample draws are only for the greedy-against-optimal comparison.
    LIBRARY, LOOPS, TOOLS = ("whatshap", "optimal"), (), ()
CAPS = (3, 5, 10, 15)
LOOP_CAPS = (3, 5)
TARGETS = {d: common.depth_seed(d)[0] for d in common.DEPTHS if d != "full"}
# Doubling the estimates for ONT's deeper slices stalled the run on 29 September:
# every job asked for 7 GB plus the reserve, while measured peaks were at most 1.6 GB.
MEMORY_SCALE = 1.0
# ONT reads, full slice plus the three subsamples, took 0.14 GB per Mb on
# chromosome 20; a download waits until the disk can hold its chromosome.
READ_GB_PER_MB = {"ont": 0.15, "hg005": 0.03}.get(common.PLATFORM)
IGNORE_READ_GROUPS = common.PLATFORM in ("ont", "hg005")
LENGTH_MB = {f"chr{i}": n for i, n in enumerate(
    [0, 249, 242, 198, 190, 182, 171, 159, 145, 138, 134, 135, 133, 114, 107, 102, 90, 83, 80,
     59, 64, 47, 51], 0) if i}


# ---- tasks run inside a job ------------------------------------------------

def task_download(chrom: str) -> None:
    import pysam

    d = common.chrom_dir(chrom)
    d.mkdir(parents=True, exist_ok=True)
    # htslib caches the remote index in the working directory; keep it per chromosome.
    os.chdir(d)
    bam = d / "full.bam"
    part = d / "full.part.bam"
    for attempt in range(1, 7):
        try:
            pysam.view("-b", "--output-fmt-option", "level=1", "-o", str(part), common.BAM_URL,
                       chrom, catch_stdout=False)
            break
        except pysam.utils.SamtoolsError as error:
            # Remote reads drop under load; a cached index from a failed read is not trusted.
            if attempt == 6:
                raise
            print(f"{chrom}: download failed ({error}); retry {attempt} in 120 s", flush=True)
            for cached in list(d.glob("*.bai")) + list(d.glob("*.csi")):
                cached.unlink()
            time.sleep(120)
    part.rename(bam)
    pysam.index(str(bam))
    handle = pysam.AlignmentFile(str(bam))
    mapped = sum(s.mapped for s in handle.get_index_statistics() if s.contig == chrom)
    length = handle.get_reference_length(chrom)
    lengths = []
    for read in handle.fetch(chrom, length // 3, 2 * length // 3):
        if read.is_secondary or read.is_supplementary or read.mapping_quality < 20:
            continue
        lengths.append(read.reference_length)
        if len(lengths) == 2000:
            break
    depth = mapped * (sum(lengths) / len(lengths)) / length
    record = {"chrom": chrom, "records": mapped,
              "mean_aligned_length": sum(lengths) / len(lengths), "length": length,
              "depth": depth}
    if common.PLATFORM == "ont":
        # ONT amendment 1: the index counts secondary alignments, which this BAM has many of.
        aligned = sum(read.reference_length or 0 for read in handle.fetch(chrom)
                      if not read.is_secondary and not read.is_unmapped)
        record.update(index_depth=depth, depth=aligned / length)
        depth = aligned / length
    (d / "depth.json").write_text(json.dumps(record))
    print(f"{chrom}: {mapped} records, depth {depth:.1f}")


def task_subsample(chrom: str, depth: str) -> None:
    import pysam

    d = common.chrom_dir(chrom)
    measured = json.loads((d / "depth.json").read_text())["depth"]
    fraction = TARGETS[depth] / measured
    if fraction >= 1.0:
        (d / f"{depth}.skip").write_text(f"slice depth {measured:.1f} is below {depth}\n")
        print(f"{chrom}: skip {depth}, slice depth {measured:.1f}")
        return
    part = d / f"{depth}.part.bam"
    pysam.view("-b", "--output-fmt-option", "level=1", "--subsample", f"{fraction:.6f}",
               "--subsample-seed", str(common.depth_seed(depth)[1]), "-o", str(part),
               str(d / "full.bam"), catch_stdout=False)
    part.rename(d / f"{depth}.bam")
    pysam.index(str(d / f"{depth}.bam"))
    print(f"{chrom}: {depth}-fold subsample, fraction {fraction:.4f}")


def task_hetsnv(chrom: str) -> None:
    extractor = common._module("extract_fragments.py")
    d = common.chrom_dir(chrom)
    d.mkdir(parents=True, exist_ok=True)
    split = common.GIAB / "split"
    sample = common.SAMPLE
    kept = extractor.het_snvs(str(split / f"{sample}_{chrom}_v4.2.1.vcf.gz"), sample, chrom,
                              str(split / f"{sample}_{chrom}_v4.2.1.bed"),
                              str(common.hetsnv_path(chrom)))
    import gzip

    with gzip.open(common.hetsnv_path(chrom), "rt") as source, \
            open(d / "hetsnv.vcf", "w") as target:
        shutil.copyfileobj(source, target)
    print(f"{chrom}: {kept} heterozygous SNVs")


def skipped(chrom: str, depth: str) -> bool:
    return depth != "full" and (common.chrom_dir(chrom) / f"{depth}.skip").exists()


def bam_of(chrom: str, depth: str) -> Path:
    return common.chrom_dir(chrom) / f"{depth}.bam"


# ---- scheduling -------------------------------------------------------------

@dataclass
class Job:
    name: str
    command: list
    deps: set = field(default_factory=set)
    kind: str = "compute"
    memory_gb: float = 1.0
    ceiling_gb: float = 4.0
    priority: tuple = ()
    chrom: str = ""
    depth: str = ""


def marker(name: str) -> Path:
    return common.DATA / "done" / (name.replace(":", "__") + ".ok")


def log_of(name: str) -> Path:
    return common.DATA / "logs" / (name.replace(":", "__") + ".log")


def self_task(*args) -> list:
    return [PYTHON, str(Path(__file__).resolve()), "--task", *args]


def build_jobs(chromosomes, include_bq=True) -> list[Job]:
    jobs: list[Job] = []
    order = sorted(chromosomes, key=lambda c: LENGTH_MB[c])
    for rank, chrom in enumerate(order):
        scale = LENGTH_MB[chrom] / 64.0
        d = common.chrom_dir(chrom)
        runs = d / "runs"
        ext = d / "external"
        ref = str(common.REFERENCE / f"{chrom}.fa")
        jobs.append(Job(f"{chrom}:hetsnv", self_task("hetsnv", chrom), kind="compute",
                        memory_gb=0.5, ceiling_gb=2, priority=(rank, 0), chrom=chrom))
        jobs.append(Job(f"{chrom}:download", self_task("download", chrom), kind="download",
                        memory_gb=0.5, ceiling_gb=2, priority=(rank, 0), chrom=chrom))
        readers = []
        for depth in common.DEPTHS:
            source = f"{chrom}:download"
            if depth != "full":
                source = f"{chrom}:subsample:{depth}"
                jobs.append(Job(source, self_task("subsample", chrom, depth),
                                deps={f"{chrom}:download"}, memory_gb=0.5, ceiling_gb=2,
                                priority=(rank, 1), chrom=chrom, depth=depth))
            for mode in common.MODES:
                name = f"{chrom}:extract:{depth}:{mode}"
                readers.append(name)
                jobs.append(Job(name, self_task("extract", chrom, depth, mode),
                                deps={source, f"{chrom}:hetsnv"},
                                memory_gb=(1.0 + 1.5 * scale) * MEMORY_SCALE, ceiling_gb=12,
                                priority=(rank, 2), chrom=chrom, depth=depth))
            for tool in TOOLS:
                name = f"{chrom}:{tool}:{depth}"
                readers.append(name)
                jobs.append(Job(name, self_task(tool, chrom, depth),
                                deps={source, f"{chrom}:hetsnv"},
                                memory_gb=(1.0 + 1.5 * scale) * MEMORY_SCALE, ceiling_gb=12,
                                priority=(rank, 2), chrom=chrom, depth=depth))
        jobs.append(Job(f"{chrom}:cleanup", self_task("cleanup", chrom), deps=set(readers),
                        memory_gb=0.2, ceiling_gb=1, priority=(rank, 3), chrom=chrom))
        for depth in common.DEPTHS:
            extract = f"{chrom}:extract:{depth}:realign"
            for pipeline in LIBRARY + LOOPS:
                caps = LOOP_CAPS if pipeline in LOOPS else CAPS
                for cap in caps:
                    name = f"{chrom}:lib:{depth}:{pipeline}:{cap}"
                    jobs.append(Job(name, self_task("lib", chrom, depth, pipeline, str(cap)),
                                    deps={extract},
                                    memory_gb=(0.5 + 0.6 * scale) * MEMORY_SCALE,
                                    ceiling_gb=6, priority=(rank, 4), chrom=chrom, depth=depth))
    for job in jobs:
        # Jobs that touch reads come first everywhere, so reads leave the disk early.
        job.priority = (0 if job.priority[1] <= 3 else 1, *job.priority)
    if not include_bq:
        return jobs
    odd = [c for c in chromosomes if c in common.ODD]
    jobs.append(Job("recalibrate", self_task("recalibrate"),
                    deps={f"{c}:extract:full:basequality" for c in odd},
                    memory_gb=3.0, ceiling_gb=8, priority=(2, 99, 5)))
    for rank, chrom in enumerate(order):
        scale = LENGTH_MB[chrom] / 64.0
        for depth in common.DEPTHS:
            for pipeline in ("optimal", "whatshap"):
                for cap in CAPS:
                    name = f"{chrom}:bq:{depth}:{pipeline}:{cap}"
                    jobs.append(Job(name, self_task("bq", chrom, depth, pipeline, str(cap)),
                                    deps={"recalibrate", f"{chrom}:extract:{depth}:basequality"},
                                    memory_gb=0.5 + 0.6 * scale, ceiling_gb=6,
                                    priority=(3, rank, 6), chrom=chrom, depth=depth))
    return jobs


def mem_available_gb() -> float:
    with open("/proc/meminfo") as handle:
        for line in handle:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1024 / 1024
    return 0.0


def free_disk_gb() -> float:
    return shutil.disk_usage(common.DATA).free / 1e9


def bam_chromosomes() -> int:
    return sum(1 for d in common.DATA.glob("chr*") if any(d.glob("*.bam")))


def read_growth_gb(starting=()) -> float:
    """Disk the chromosomes holding reads, or in `starting`, are still expected to take.

    `starting` names chromosomes whose download is running but may not have
    written a file yet; without them two downloads launched together each
    saw the other's disk as free (30 September).
    """
    total = 0.0
    for d in common.DATA.glob("chr*"):
        bams = list(d.glob("*.bam"))
        if bams or d.name in starting:
            used = sum(b.stat().st_size for b in bams) / 1e9
            total += max(0.0, READ_GB_PER_MB * LENGTH_MB[d.name] - used)
    return total


def run_scheduler(jobs, max_jobs, min_free_gb, max_bam_chromosomes, reserve_gb) -> int:
    (common.DATA / "done").mkdir(parents=True, exist_ok=True)
    (common.DATA / "logs").mkdir(parents=True, exist_ok=True)
    # No download is running yet, so any partial read file is left from a killed run. Left in
    # place it counts toward max_bam_chromosomes and stalls every download (27 November).
    for stale in common.DATA.glob("chr*/full.part.bam"):
        stale.unlink()
    by_name = {j.name: j for j in jobs}
    done = {j.name for j in jobs if marker(j.name).exists()}
    failed: set = set()
    running: dict = {}
    pending = sorted((j for j in jobs if j.name not in done), key=lambda j: j.priority)
    env = dict(os.environ, PYTHONPATH=f"{common.ROOT / 'src'}:{HERE}")
    last_status = 0.0
    while pending or running:
        for name, (process, handle, started) in list(running.items()):
            if process.poll() is None:
                continue
            handle.close()
            del running[name]
            if process.returncode == 0:
                marker(name).write_text(f"{time.time() - started:.1f} s\n")
                done.add(name)
            else:
                failed.add(name)
        blocked = {j.name for j in pending if j.deps & failed}
        failed |= blocked
        pending = [j for j in pending if j.name not in blocked]
        launched = True
        while launched:
            launched = False
            if len(running) >= max_jobs:
                break
            downloads = sum(1 for n in running if by_name[n].kind == "download")
            max_downloads = int(os.environ.get("PHASEBUDGET_MAX_DOWNLOADS", "2"))
            for job in pending:
                if not job.deps <= done:
                    continue
                if job.kind == "download":
                    if downloads >= max_downloads or free_disk_gb() < min_free_gb:
                        continue
                    if bam_chromosomes() >= max_bam_chromosomes:
                        continue
                    starting = {by_name[n].chrom for n in running if by_name[n].kind == "download"}
                    if READ_GB_PER_MB and (free_disk_gb() - read_growth_gb(starting)
                                           < READ_GB_PER_MB * LENGTH_MB[job.chrom] + 5.0):
                        continue
                if READ_GB_PER_MB and ":subsample:" in job.name:
                    # A subsample writes up to its share of the full slice; keep 5 GB spare.
                    full = common.chrom_dir(job.chrom) / "full.bam"
                    measured = json.loads((common.chrom_dir(job.chrom) / "depth.json").read_text())
                    size = full.stat().st_size / 1e9 * TARGETS[job.depth] / measured["depth"] \
                        if full.exists() else 0.0
                    if free_disk_gb() < size + 5.0:
                        continue
                # With nothing running, an estimate above the free memory would stall the run
                # for good; the job's own memory ceiling still protects the machine.
                if running and mem_available_gb() < job.memory_gb + reserve_gb:
                    continue
                handle = open(log_of(job.name), "w")
                command = ["systemd-run", "--user", "--scope", "-q",
                           "-p", f"MemoryMax={job.ceiling_gb}G", "-p", "MemorySwapMax=0",
                           "nice", "-n", "19", "ionice", "-c", "3", *job.command]
                process = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT,
                                           env=env, cwd=common.ROOT)
                running[job.name] = (process, handle, time.time())
                pending.remove(job)
                launched = True
                break
        if time.time() - last_status > 30:
            status = {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "done": len(done),
                      "running": sorted(running), "pending": len(pending),
                      "failed": sorted(failed), "free_disk_gb": round(free_disk_gb(), 1),
                      "mem_available_gb": round(mem_available_gb(), 1)}
            (common.DATA / "status.json").write_text(json.dumps(status, indent=1))
            last_status = time.time()
        time.sleep(2)
    (common.DATA / "status.json").write_text(json.dumps(
        {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "done": len(done), "failed": sorted(failed),
         "finished": True}, indent=1))
    return len(failed)


# ---- job commands -------------------------------------------------------------

def dispatch(args) -> None:
    task, rest = args.task[0], args.task[1:]
    extractor_script = str(HERE / "extract_fragments.py")
    pipelines_script = str(HERE / "hg002_pipelines.py")
    if task == "download":
        task_download(rest[0])
    elif task == "hetsnv":
        task_hetsnv(rest[0])
    elif task == "subsample":
        task_subsample(rest[0], rest[1])
    elif task == "extract":
        chrom, depth, mode = rest
        if skipped(chrom, depth):
            print("skipped depth")
            return
        subprocess.run([PYTHON, extractor_script, "--bam", str(bam_of(chrom, depth)),
                        "--vcf", str(common.hetsnv_path(chrom)), "--reference",
                        str(common.REFERENCE / f"{chrom}.fa"), "--chromosome", chrom,
                        "--alleles", mode, "--out",
                        str(common.fragments_path(chrom, depth, mode)),
                        "--sample", common.SAMPLE,
                        *(["--ignore-read-groups"] if IGNORE_READ_GROUPS else [])],
                       check=True)
    elif task in ("whatshap", "hapcut2", "longphase"):
        chrom, depth = rest
        if skipped(chrom, depth):
            print("skipped depth")
            return
        d = common.chrom_dir(chrom)
        ext = d / "external"
        ext.mkdir(exist_ok=True)
        ref = str(common.REFERENCE / f"{chrom}.fa")
        bam = str(bam_of(chrom, depth))
        if task == "longphase":
            subprocess.run(["/usr/bin/time", "-v", str(Path.home() / ".local" / "bin" / "longphase"),
                            "phase", "-s", str(d / "hetsnv.vcf"), "-b", bam, "-r", ref, "--ont",
                            "-t", "1", "-o", str(ext / f"longphase_{depth}")], check=True)
        elif task == "whatshap":
            whatshap = str(Path(PYTHON).with_name("whatshap"))
            subprocess.run(["/usr/bin/time", "-v", whatshap, "phase", "--reference", ref,
                            *(["--ignore-read-groups"] if IGNORE_READ_GROUPS else []),
                            "-o", str(ext / f"whatshap_{depth}.vcf"),
                            str(common.hetsnv_path(chrom)), bam], check=True)
        else:
            home = Path.home() / ".local" / "bin"
            frag = str(ext / f"hapcut2_{depth}.frag")
            platform = "--ont" if common.PLATFORM == "ont" else "--pacbio"
            subprocess.run(["/usr/bin/time", "-v", str(home / "extractHAIRS"), "--bam", bam,
                            "--VCF", str(d / "hetsnv.vcf"), platform, "1", "--ref", ref,
                            "--out", frag], check=True)
            subprocess.run(["/usr/bin/time", "-v", str(home / "HAPCUT2"), "--fragments", frag,
                            "--VCF", str(d / "hetsnv.vcf"), "--output",
                            str(ext / f"hapcut2_{depth}.blocks")], check=True)
            Path(frag).unlink()
    elif task == "cleanup":
        d = common.chrom_dir(rest[0])
        for path in list(d.glob("*.bam")) + list(d.glob("*.bai")):
            path.unlink()
        print(f"{rest[0]}: reads deleted")
    elif task in ("lib", "bq"):
        chrom, depth, pipeline, cap = rest
        mode = "realign" if task == "lib" else "basequality"
        fragments = common.fragments_path(chrom, depth, mode)
        if not fragments.exists():
            print("no fragments (skipped depth)")
            return
        runs = common.chrom_dir(chrom) / "runs"
        runs.mkdir(exist_ok=True)
        weights = "logodds" if task == "lib" else "recalibrated"
        command = [PYTHON, pipelines_script, "--fragments", str(fragments), "--pipeline",
                   pipeline, "--cap", cap, "--weights", weights, "--out",
                   str(runs / f"{task}_{depth}_{pipeline}_c{cap}.json.gz")]
        if task == "bq":
            command += ["--table", str(common.DATA / "recalibration.json")]
        if pipeline == "optimal" and cap == "5":
            command.append("--signals")
        subprocess.run(command, check=True)
    elif task == "recalibrate":
        subprocess.run([PYTHON, str(HERE / "recalibrate.py")], check=True)
    else:
        raise SystemExit(f"unknown task {task}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chromosomes", nargs="+", default=list(common.CONFIRMATORY))
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--min-free-gb", type=float, default=20.0)
    parser.add_argument("--max-bam-chromosomes", type=int, default=3)
    parser.add_argument("--reserve-gb", type=float, default=3.0)
    parser.add_argument("--skip-bq", action="store_true",
                        help="leave out the recalibration and the base-quality pipelines")
    parser.add_argument("--task", nargs="+")
    args = parser.parse_args()
    if args.task:
        dispatch(args)
        return
    common.DATA.mkdir(parents=True, exist_ok=True)
    jobs = build_jobs(args.chromosomes, include_bq=not args.skip_bq and common.PLATFORM == "hifi"
                      and not common.SEEDS)
    print(f"{common.PLATFORM}: {len(jobs)} jobs over {', '.join(args.chromosomes)}", flush=True)
    failures = run_scheduler(jobs, args.jobs, args.min_free_gb, args.max_bam_chromosomes,
                             args.reserve_gb)
    print(f"finished with {failures} failed jobs", flush=True)


if __name__ == "__main__":
    main()
