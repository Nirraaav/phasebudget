"""Data steps of `docs/design/confidence-run-plan.md`, one chromosome at a time.

With `PHASEBUDGET_PLATFORM=hg001_ont` or `ont_hac`, for each confirmatory
chromosome, smallest first:
1. het SNVs: Platinum NA12878 inside its confident regions (HG001), or the
   ONT run's `hetsnv.vcf.gz` (HG002 hac);
2. estimate the depth of the remote BAM as the median non-secondary aligned
   depth of ten 500 kb windows spaced over the chromosome (amendment 1);
3. draw the 30-fold slice from the remote BAM (seed 0); measure it as the
   median depth over 1 Mb windows; draw the 20- and 10-fold slices from it
   with fractions 20 and 10 over that depth and seeds 1 and 2, which are
   independent of the first draw; and measure each slice;
4. at each depth: extract `realign` fragments as `whatshap phase` reads
   them, and run `whatshap phase` and LongPhase with the calibration runs'
   commands;
5. compress the outputs and delete the BAMs.

Each step writes a marker, so the script resumes where it stopped. With
`--workers k`, k chromosomes are processed at once, each in its own process.
Before a chromosome's draws, its process reserves the chromosome's peak disk
use in `data/reservations/` and only starts when free disk minus every other
reservation covers it plus a margin, so runs of several platforms and
chromosomes download side by side, each on its own connection, without
overfilling the disk. A reservation is released after cleanup; a run stopped
for good leaves its files behind, and they must be deleted by hand. At most
`PHASE_SLOTS` chromosomes, over all runs, phase at once, which bounds memory
(`PHASEBUDGET_PHASE_SLOTS` and `PHASEBUDGET_SLOT_FREE_GB` override the defaults).

With `PHASEBUDGET_PLATFORM=hprc_<sample>` (round 4, W1,
`docs/design/w1-ancestry-plan.md`) the reads are HPRC HiFi aligned with
winnowmap, LongPhase runs with `--pb`, and the variants are HG002's v5.0q
calls or, for a test genome, dipcall's calls on its assemblies inside
dipcall's BED, both treated as below.

With `PHASEBUDGET_PLATFORM=phifi_<sample>` (W1b, `docs/design/w1b-pooled-plan.md`) the
reads are the Platinum Pedigree HiFi alignments, LongPhase runs with `--pb`, and the
variants are the Pedigree truthset's, as for `ped_<sample>`.

With `PHASEBUDGET_PLATFORM=cal_r10` or `ped_<sample>` (track B of
the extension plan (development repository)), and with `cal_v5` or `hg001_v5`
(`docs/design/third-test-plan.md`), the variants are heterozygous SNVs and
indels of at most 50 bases, with their phase removed: GIAB v5.0q inside its
BED for HG002, or the Platinum Pedigree truthset inside its high-quality
regions (NA12878 in `hg001_v5` too). Fragments keep indels, and LongPhase runs with `--indels`;
`whatshap phase` phases indels by default.

Usage: `PHASEBUDGET_PLATFORM=hg001_ont PYTHONPATH=src:benchmarks .venv/bin/python
benchmarks/confidence_run.py`
"""

from __future__ import annotations

import argparse
import fcntl
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402

LENGTH_MB = {f"chr{i}": n for i, n in enumerate(
    [0, 249, 242, 198, 190, 182, 171, 159, 145, 138, 134, 135, 133, 114, 107, 102, 90, 83, 80,
     59, 64, 47, 51], 0) if i}
# Peak disk per Mb: the 36-fold draw and the 30-fold trim side by side. Measured on HG002 sup
# chromosome 21 (2.29 GB drawn at level 1 for 46.7 Mb, 0.049 GB per Mb) as 0.049 + 0.041.
PEAK_GB_PER_MB = 0.09
FREE_MARGIN_GB = 3
# Chromosome 1 at 30x peaks near 5 GB over its three phasing jobs, so a slot beyond the second
# is taken only while the machine has SLOT_FREE_GB available. Neither changes any output.
PHASE_SLOTS = int(os.environ.get("PHASEBUDGET_PHASE_SLOTS", "4"))
SLOT_FREE_GB = float(os.environ.get("PHASEBUDGET_SLOT_FREE_GB", "6"))
SEEDS = {"20": "1", "10": "2"}  # amendment 1: a seed reused across draws nests the subsamples
WHATSHAP = str(Path(sys.executable).with_name("whatshap"))
LONGPHASE = str(Path.home() / ".local" / "bin" / "longphase")


def done(chrom: str, step: str) -> Path:
    return common.DATA / "done" / f"{chrom}__{step}.ok"


def step(chrom: str, name: str, action) -> None:
    marker = done(chrom, name)
    if marker.exists():
        return
    started = time.time()
    action()
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(f"{time.time() - started:.1f} s\n")
    print(f"{chrom}: {name} done in {time.time() - started:.0f} s", flush=True)


def hetsnv(chrom: str) -> None:
    d = common.chrom_dir(chrom)
    d.mkdir(parents=True, exist_ok=True)
    if common.EXTENSION:
        extractor = common._module("extract_fragments.py")
        if common.HPRC_TEST or common.PLATFORM == "k3_hg005":
            vcf = common.HPRC / common.SAMPLE / f"{common.SAMPLE}.dip.vcf.gz"
            bed = common.HPRC / common.SAMPLE / f"{common.SAMPLE}.dip.bed"
        elif common.PLATFORM in ("cal_r10", "cal_v5", "hprc_HG002", "w4_grch38", "w4_chm13"):
            vcf = common.GIAB_V5 / f"HG002_{common.V5_BUILD}_v5.0q_smvar.vcf.gz"
            bed = common.GIAB_V5 / f"HG002_{common.V5_BUILD}_v5.0q_smvar.benchmark.bed"
        elif common.PLATFORM == "q3_q100":
            vcf, bed = common.Q100_VCF, common.Q100_BED
        elif common.PLATFORM in ("hg005_ont", "k1_hg005"):
            vcf = common.GIAB / "HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
            bed = common.GIAB / "HG005_GRCh38_1_22_v4.2.1_benchmark.bed"
        else:
            vcf = common.PEDIGREE / "CEPH1463.GRCh38.family-truthset.ov.vcf.gz"
            bed = common.PEDIGREE / "hq_regions_final.bed"
            if not bed.exists():
                with gzip.open(common.PEDIGREE / "hq_regions_final.bed.gz", "rb") as source, \
                        open(bed.with_suffix(".part"), "wb") as target:
                    shutil.copyfileobj(source, target)
                bed.with_suffix(".part").rename(bed)
        extractor.het_snvs(str(vcf), common.VCF_SAMPLE, chrom, str(bed),
                           str(common.hetsnv_path(chrom)), indels=True, unphase=True)
    elif common.PLATFORM == "hg001_ont":
        extractor = common._module("extract_fragments.py")
        bed = common.PLATINUM / "ConfidentRegions.bed"
        if not bed.exists():
            with gzip.open(common.PLATINUM / "ConfidentRegions.bed.gz", "rb") as source, \
                    open(bed.with_suffix(".part"), "wb") as target:
                shutil.copyfileobj(source, target)
            bed.with_suffix(".part").rename(bed)
        extractor.het_snvs(str(common.PLATINUM / "NA12878.vcf.gz"), "NA12878", chrom,
                           str(bed), str(common.hetsnv_path(chrom)))
    else:
        source = common.ROOT / "data" / "genome_ont" / chrom / "hetsnv.vcf.gz"
        shutil.copy(source, common.hetsnv_path(chrom))
        shutil.copy(str(source) + ".tbi", str(common.hetsnv_path(chrom)) + ".tbi")
    with gzip.open(common.hetsnv_path(chrom), "rt") as source, open(d / "hetsnv.vcf", "w") as target:
        shutil.copyfileobj(source, target)
    for leftover in d.glob("*.input.vcf"):
        leftover.unlink()


def aligned_depth(bam: str, chrom: str, start: int, end: int) -> float:
    import pysam

    with pysam.AlignmentFile(bam) as handle:
        bases = 0
        for read in handle.fetch(chrom, start, end):
            if read.is_secondary or read.is_unmapped:
                continue
            bases += max(0, min(read.reference_end, end) - max(read.reference_start, start))
    return bases / (end - start)


def median_depth(bam: str, chrom: str, window: int = 1_000_000) -> float:
    """Median over `window`-bp windows with any aligned base of non-secondary aligned depth."""
    import statistics

    import pysam

    with pysam.AlignmentFile(bam) as handle:
        length = handle.get_reference_length(chrom)
        bases = [0] * (length // window + 1)
        for read in handle.fetch(chrom):
            if read.is_secondary or read.is_unmapped:
                continue
            start, end = read.reference_start, read.reference_end
            while start < end:
                w = start // window
                stop = min(end, (w + 1) * window)
                bases[w] += stop - start
                start = stop
    covered = [b / window for b in bases[:-1] if b > 0]
    return statistics.median(covered)


def index_dir(d: Path, k: int) -> Path:
    """Where htslib caches the index of the `k`-th remote BAM, which it names by basename:
    two flowcells both called `calls.sorted.bam` must not share one cached index."""
    if len(common.BAM_URLS) == 1:
        return d
    out = d / f"remote{k}"
    out.mkdir(exist_ok=True)
    return out


def estimate(chrom: str) -> None:
    import statistics

    d = common.chrom_dir(chrom)
    here = os.getcwd()
    length = LENGTH_MB[chrom] * 1_000_000
    windows = [0.0] * 10
    try:
        for k, url in enumerate(common.BAM_URLS):  # depths of several BAMs add
            os.chdir(index_dir(d, k))  # htslib caches the remote index in the working directory
            for w in range(10):
                windows[w] += aligned_depth(url, chrom, int((w + 0.5) * length / 10) - 250_000,
                                            int((w + 0.5) * length / 10) + 250_000)
    finally:
        os.chdir(here)
    covered = [w for w in windows if w > 0]
    (d / "estimate.json").write_text(json.dumps({"chrom": chrom, "windows": windows,
                                                 "median_depth": statistics.median(covered)}))


def slice30(chrom: str) -> None:
    import pysam

    d = common.chrom_dir(chrom)
    estimated = json.loads((d / "estimate.json").read_text())["median_depth"]
    # Amendment 1, addendum: the seed-0 remote draw overshoots, so aim above 30 and trim locally.
    fraction = 36 / estimated
    part = d / "30.part.bam"
    drawn = d / "30.drawn.bam"
    if not drawn.exists():  # the rename below happens only after a complete draw
        draw(common.BAM_URLS, chrom, fraction, d, part)
        part.rename(drawn)
    pysam.index(str(drawn))
    measured = median_depth(str(drawn), chrom)
    if measured > 30:
        trim = d / "30.part.bam"
        pysam.view("-b", "--output-fmt-option", "level=1", "--subsample", f"{30 / measured:.6f}",
                   "--subsample-seed", "3", "-o", str(trim), str(drawn), catch_stdout=False)
        trim.rename(d / "30.bam")
        for path in (drawn, Path(str(drawn) + ".bai")):
            path.unlink()
    else:
        drawn.rename(d / "30.bam")
        Path(str(drawn) + ".bai").unlink()
    pysam.index(str(d / "30.bam"))


def draw(urls, chrom: str, fraction: float, d: Path, out: Path) -> None:
    """The seed-0 draw of `chrom` from each BAM in `urls`, merged into `out` when there are several."""
    import pysam

    args = ["--subsample", f"{fraction:.6f}", "--subsample-seed", "0"] if fraction < 1 else []
    parts = [out] if len(urls) == 1 else [d / f"30.part{k}.bam" for k in range(len(urls))]
    here = os.getcwd()
    try:
        for k, (url, part) in enumerate(zip(urls, parts)):
            os.chdir(index_dir(d, k))
            pysam.view("-b", "--output-fmt-option", "level=1", *args, "-o", str(part), url, chrom,
                       catch_stdout=False)
    finally:
        os.chdir(here)
    if len(parts) > 1:
        pysam.merge("-f", "--output-fmt-option", "level=1", str(out), *map(str, parts),
                    catch_stdout=False)
        for part in parts:
            part.unlink()


def subsample(chrom: str) -> None:
    import pysam

    d = common.chrom_dir(chrom)
    measured = median_depth(str(d / "30.bam"), chrom)
    for depth, seed in SEEDS.items():
        fraction = min(1.0, float(depth) / measured)
        part = d / f"{depth}.part.bam"
        pysam.view("-b", "--output-fmt-option", "level=1", "--subsample", f"{fraction:.6f}",
                   "--subsample-seed", seed, "-o", str(part), str(d / "30.bam"),
                   catch_stdout=False)
        part.rename(d / f"{depth}.bam")
        pysam.index(str(d / f"{depth}.bam"))


def measure(chrom: str) -> None:
    d = common.chrom_dir(chrom)
    record = json.loads((d / "estimate.json").read_text())
    for depth in common.DEPTHS:
        record[f"depth_{depth}"] = median_depth(str(d / f"{depth}.bam"), chrom)
    (d / "depth.json").write_text(json.dumps(record))


def run(cmd: list, log: Path) -> None:
    with open(log, "w") as handle:
        subprocess.run(["nice", "-n", "15", "/usr/bin/time", "-v", *cmd], stdout=handle,
                       stderr=subprocess.STDOUT, check=True)


def phase(chrom: str, depth: str) -> None:
    d = common.chrom_dir(chrom)
    ext = d / "external"
    ext.mkdir(exist_ok=True)
    ref = str(common.REFERENCE / f"{chrom}.fa")
    bam = str(d / f"{depth}.bam")
    jobs = {
        "fragments": [sys.executable, str(HERE / "extract_fragments.py"), "--bam", bam,
                      "--vcf", str(common.hetsnv_path(chrom)), "--reference", ref,
                      "--chromosome", chrom, "--sample", common.VCF_SAMPLE,
                      "--out", str(common.fragments_path(chrom, depth, "realign")),
                      "--ignore-read-groups"],
        "whatshap": [WHATSHAP, "phase", "--reference", ref, "--ignore-read-groups",
                     "-o", str(ext / f"whatshap_{depth}.vcf"), str(common.hetsnv_path(chrom)), bam],
        "longphase": [LONGPHASE, "phase", "-s", str(d / "hetsnv.vcf"), "-b", bam, "-r", ref,
                      "--pb" if common.PLATFORM.startswith(("hprc_", "phifi_", "w4_", "q3_", "k1_", "k3_"))
                      else "--ont", "-t", "1",
                      "-o", str(ext / f"longphase_{depth}")],
    }
    if common.EXTENSION:
        jobs["fragments"].append("--indels")
        jobs["longphase"].append("--indels")
    if common.PLATFORM == "ont_hac":
        # Amendment 2: HapCUT2 on HG002 hac only, with the calibration runs' commands.
        home = Path.home() / ".local" / "bin"
        frag = str(ext / f"hapcut2_{depth}.frag")
        jobs["hapcut2"] = [
            [str(home / "extractHAIRS"), "--bam", bam, "--VCF", str(d / "hetsnv.vcf"), "--ont",
             "1", "--ref", ref, "--out", frag],
            [str(home / "HAPCUT2"), "--fragments", frag, "--VCF", str(d / "hetsnv.vcf"),
             "--output", str(ext / f"hapcut2_{depth}.blocks")],
            ["rm", "-f", frag],
        ]
    todo = {n: c for n, c in jobs.items() if not done(chrom, f"{n}_{depth}").exists()}

    def one(item):
        name, cmd = item
        commands = cmd if isinstance(cmd[0], list) else [cmd]

        def action():
            for k, command in enumerate(commands):
                run(command, ext / f"{name}_{depth}.{k}.log")

        step(chrom, f"{name}_{depth}", action)

    with ThreadPoolExecutor(3) as pool:
        list(pool.map(one, todo.items()))
    for name in ("whatshap", "longphase"):
        vcf = ext / f"{name}_{depth}.vcf"
        if vcf.exists():
            subprocess.run(["gzip", "-f", str(vcf)], check=True)


def retried(action, attempts: int = 6, wait: int = 120) -> None:
    """Run a step that reads the remote BAM, retrying after dropped connections."""
    for attempt in range(1, attempts + 1):
        try:
            action()
            return
        except (OSError, ConnectionError) as error:
            if attempt == attempts:
                raise
            print(f"network error ({error}); retry {attempt} in {wait} s", flush=True)
            time.sleep(wait)


def cleanup(chrom: str) -> None:
    d = common.chrom_dir(chrom)
    for pattern in ("*.bam", "*.bai", "*.csi", "remote*/*.bai", "remote*/*.csi"):
        for path in d.glob(pattern):
            path.unlink()


RESERVATIONS = common.ROOT / "data" / "reservations"


def reservation(chrom: str) -> Path:
    return RESERVATIONS / f"{common.PLATFORM}__{chrom}.gb"


def _unwritten(entry: Path) -> float:
    """GB a reservation may still add: its size less the BAMs already in its directory."""
    try:
        fields = entry.read_text().split("\t")
        size = float(fields[0])
        if len(fields) < 2:
            return size
        written = sum(p.stat().st_size for p in Path(fields[1].strip()).glob("*.bam")) / 1e9
    except FileNotFoundError:  # another run released it, or cleaned up, mid-read
        return 0.0
    return max(0.0, size - written)


def reserve(chrom: str) -> None:
    """Wait until the disk holds this chromosome's slices beside every other reservation.

    Another reservation counts only for what it has not written yet, since
    what it has written is already missing from the free space.
    """
    RESERVATIONS.mkdir(parents=True, exist_ok=True)
    own = reservation(chrom)
    size = LENGTH_MB[chrom] * PEAK_GB_PER_MB
    drawn = done(chrom, "slice30").exists()
    while True:
        with open(RESERVATIONS / ".lock", "w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            others = sum(_unwritten(f) for f in RESERVATIONS.glob("*.gb") if f != own)
            free = shutil.disk_usage(RESERVATIONS).free / 1e9
            # A chromosome already drawn holds its disk; it only re-registers.
            if drawn or free - others >= size + FREE_MARGIN_GB:
                own.write_text(f"{size:.1f}\t{common.chrom_dir(chrom)}\n")
                return
        print(f"{chrom}: waiting for {size + FREE_MARGIN_GB:.0f} GB beyond {others:.0f} GB "
              f"still to be written by other runs ({free:.0f} GB free)", flush=True)
        time.sleep(120)


def available_gb() -> float:
    with open("/proc/meminfo") as handle:
        for line in handle:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1e6
    return 0.0


def phase_slot():
    """Block until one of `PHASE_SLOTS` lock files is free and memory allows; return it, held."""
    RESERVATIONS.mkdir(parents=True, exist_ok=True)
    while True:
        for k in range(PHASE_SLOTS):
            if k >= 2 and available_gb() < SLOT_FREE_GB:  # two slots always; more if memory allows
                break
            handle = open(RESERVATIONS / f".phase{k}.lock", "w")
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                return handle
            except BlockingIOError:
                handle.close()
        time.sleep(30)


def process(chrom: str, attempts: int = 3) -> bool:
    """Every step of one chromosome, tried `attempts` times; False if every attempt failed.

    After the last failure the chromosome's reads are deleted and its reservation
    released, so the other chromosomes are not left waiting for its disk.
    """
    for attempt in range(1, attempts + 1):
        if attempt > 1:
            time.sleep(300)
        if attempt_once(chrom):
            return True
        print(f"{chrom}: attempt {attempt} of {attempts} failed", flush=True)
    cleanup(chrom)
    for name in ("slice30", "subsample", "measure"):  # their outputs are gone; redraw next time
        done(chrom, name).unlink(missing_ok=True)
    reservation(chrom).unlink(missing_ok=True)
    return False


def attempt_once(chrom: str) -> bool:
    import traceback

    try:
        step(chrom, "hetsnv", lambda: hetsnv(chrom))
        step(chrom, "estimate", lambda: retried(lambda: estimate(chrom)))
        reserve(chrom)
        step(chrom, "slice30", lambda: retried(lambda: slice30(chrom)))
        step(chrom, "subsample", lambda: subsample(chrom))
        step(chrom, "measure", lambda: measure(chrom))
        with phase_slot():
            for depth in common.DEPTHS:
                phase(chrom, depth)
        step(chrom, "cleanup", lambda: cleanup(chrom))
        reservation(chrom).unlink(missing_ok=True)
        return True
    except Exception:
        print(f"{chrom}: failed\n{traceback.format_exc()}", flush=True)
        return False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chromosomes", nargs="+", default=list(common.CONFIRMATORY))
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if common.PLATFORM not in ("hg001_ont", "ont_hac") and not common.EXTENSION:
        raise SystemExit("set PHASEBUDGET_PLATFORM to hg001_ont, ont_hac, cal_r10 or ped_<sample>")
    todo = [c for c in sorted(args.chromosomes, key=lambda c: LENGTH_MB[c])
            if not done(c, "cleanup").exists()]
    # Processes, not threads: the remote steps chdir so htslib caches each index per chromosome.
    with ProcessPoolExecutor(args.workers) as pool:
        failed = [c for c, ok in zip(todo, pool.map(process, todo)) if not ok]
    if failed:
        raise SystemExit(f"failed: {' '.join(failed)}")


if __name__ == "__main__":
    main()
