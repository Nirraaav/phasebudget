"""Switch error inside GIAB genome stratifications, against GIAB v5.0q.

A mock review of the selection paper asked where the errors
fall. This script scores stored genome-run phasings inside the v5.0q
benchmark regions (`genome_common.v5_truth`) and, within those, inside each
GIAB GRCh38 stratification below. A pair of consecutive compared variants
belongs to a stratum when both variants lie inside it; the stratum "all" is
every pair. Per stratum it prints, for each pipeline, the comparable pairs,
the switch errors and their rate per 10,000; and for the first two pipelines
the paired counts on common pairs (`paired.common_pairs`) with the exact
sign test and the 1 Mb window sign-flip test of `permutation.py`.

The BEDs are GIAB genome stratifications v3.6, GRCh38, downloaded to
`data/strata/` and checked against these md5 sums (from
`GRCh38-genome-stratifications-md5s.txt` in the same directory):

  base: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/
  72ec3f7db18ef9318c64b1df271703fe  SegmentalDuplications/GRCh38_segdups.bed.gz
  8a78766bb99abf1d25a5b307d49a0324  Mappability/GRCh38_lowmappabilityall.bed.gz
  05a302a437da3f3514106d3eedd62696  LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz
  562ebfd3b3994ae50c3d162d96c66e94  Union/GRCh38_notinalldifficultregions.bed.gz

The platform comes from PHASEBUDGET_PLATFORM, as in `genome_common.py`.
The analysis is post hoc: no plan registered it.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402
import permutation  # noqa: E402
import score_genome  # noqa: E402

from phasebudget.metrics import _compared_members  # noqa: E402
from phasebudget.paired import common_pairs, sign_test, window_discordance  # noqa: E402

BASE = ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
        "genome-stratifications/v3.6/GRCh38@all/")
STRATA = (  # (name, path under BASE, md5)
    ("segdups", "SegmentalDuplications/GRCh38_segdups.bed.gz",
     "72ec3f7db18ef9318c64b1df271703fe"),
    ("lowmap", "Mappability/GRCh38_lowmappabilityall.bed.gz",
     "8a78766bb99abf1d25a5b307d49a0324"),
    ("TR+homopol", "LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz",
     "05a302a437da3f3514106d3eedd62696"),
    ("easy", "Union/GRCh38_notinalldifficultregions.bed.gz",
     "562ebfd3b3994ae50c3d162d96c66e94"),
)
DIR = common.ROOT / "data" / "strata"
MIN_FREE = 25 * 2 ** 30
WINDOW = 1_000_000
RESAMPLES = 20000
SEED = 0


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch() -> dict[str, Path]:
    """Download any missing BED, verify every md5, return `{name: path}`."""
    DIR.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name, rel, md5 in STRATA:
        path = DIR / Path(rel).name
        if not path.exists():
            if shutil.disk_usage(DIR).free < MIN_FREE:
                raise SystemExit("under 25 GB free on disk; not downloading")
            subprocess.run(["curl", "-sSfL", "-o", str(path) + ".part", BASE + rel], check=True)
            Path(str(path) + ".part").rename(path)
        if _md5(path) != md5:
            raise SystemExit(f"{path}: md5 mismatch")
        paths[name] = path
    return paths


def load_bed(path: Path, chroms) -> dict[str, tuple[list[int], list[int]]]:
    """Merged `(starts, ends)` per chromosome, 0-based half-open."""
    raw: dict[str, list[tuple[int, int]]] = {c: [] for c in chroms}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            f = line.split("\t", 3)
            if f[0] in raw:
                raw[f[0]].append((int(f[1]), int(f[2])))
    merged = {}
    for chrom, intervals in raw.items():
        starts, ends = [], []
        for s, e in sorted(intervals):
            if ends and s <= ends[-1]:
                ends[-1] = max(ends[-1], e)
            else:
                starts.append(s)
                ends.append(e)
        merged[chrom] = (starts, ends)
    return merged


def inside(bed: tuple[list[int], list[int]], positions) -> list[bool]:
    starts, ends = bed
    out = []
    for p in positions:
        k = bisect.bisect_right(starts, p) - 1
        out.append(k >= 0 and p < ends[k])
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--depth", default="full")
    parser.add_argument("--cap", type=int, default=3)
    parser.add_argument("--pipelines", nargs="+", default=["whatshap", "optimal"],
                        help="the first two are also compared paired, first against second; "
                             "tool names such as 'whatshap phase' read the tool outputs")
    parser.add_argument("--chroms", nargs="+", default=list(common.CONFIRMATORY))
    args = parser.parse_args()
    names = args.pipelines
    first, second = names[0], names[1] if len(names) > 1 else None
    beds = {name: load_bed(path, args.chroms) for name, path in fetch().items()}
    strata = ("all",) + tuple(beds)

    per = {s: {n: [0, 0] for n in names} for s in strata}  # [pairs, switches]
    paired = {s: {"pairs": 0, "first": 0, "second": 0, "both": 0, "windows": []} for s in strata}
    sites = {s: 0 for s in strata}
    objective = {n: [0.0, 0.0] for n in names}
    for chrom in args.chroms:
        variants = score_genome.score.variant_table(common.hetsnv_path(chrom))
        positions = [p - 1 for p, _, _ in variants]
        truth = common.v5_truth(chrom, variants)
        member = {"all": [True] * len(positions)}
        member.update({s: inside(beds[s][chrom], positions) for s in beds})
        for s in strata:
            sites[s] += sum(t is not None and m for t, m in zip(truth, member[s]))
        phasings = {}
        for n in names:
            # Tool outputs (`whatshap phase`, `LongPhase`, `HapCUT2`) are uncapped runs.
            kind = "tool" if n in score_genome.TOOLS else "lib"
            phasing, record = score_genome.load_phasing(kind, chrom, args.depth, n, args.cap)
            if len(phasing) != len(positions):
                raise ValueError(f"{chrom} {n}: wrong variant count")
            phasings[n] = phasing
            if record.get("objective") is not None:
                objective[n][0] += record["objective"]
                objective[n][1] += record["optimum"]
            for members in _compared_members(truth, phasing):
                for left, right in zip(members, members[1:]):
                    wrong = (phasing.alleles[left] ^ phasing.alleles[right]) != \
                        (truth[left] ^ truth[right])
                    for s in strata:
                        if member[s][left] and member[s][right]:
                            per[s][n][0] += 1
                            per[s][n][1] += wrong
        if second is None:
            continue
        a, b = phasings[first], phasings[second]
        usable = [j for j in range(len(truth)) if truth[j] is not None
                  and a.alleles[j] is not None and b.alleles[j] is not None]
        previous = dict(zip(usable[1:], usable))
        rows = common_pairs(truth, a, b)
        for s in strata:
            kept = [r for r in rows if member[s][previous[r[0]]] and member[s][r[0]]]
            q = paired[s]
            q["pairs"] += len(kept)
            q["first"] += sum(1 for _, x, y in kept if x and not y)
            q["second"] += sum(1 for _, x, y in kept if y and not x)
            q["both"] += sum(1 for _, x, y in kept if x and y)
            q["windows"].extend(window_discordance(kept, positions, WINDOW))

    print("switch error by GIAB v3.6 GRCh38 stratification, inside the v5.0q benchmark "
          "regions (post hoc)")
    print(f"platform {common.PLATFORM}; depth {args.depth}; cap {args.cap}; weights logodds; "
          f"chromosomes {', '.join(args.chroms)}; pipelines {', '.join(names)}")
    print("a pair is in a stratum when both its variants are; windows "
          f"{WINDOW // 1_000_000} Mb; sign-flip test {RESAMPLES} resamples, seed {SEED}, "
          f"exact at <= {permutation.EXACT_LIMIT} non-zero windows")
    print("objective / certified optimum, pooled: " + ", ".join(
        f"{n} {o / m:.4f}" if m else f"{n} n/a" for n, (o, m) in objective.items()))
    print()
    head = f"{'stratum':>11} {'v5 sites':>9}" + "".join(
        f" | {n[:12]:>12} {'pairs':>8} {'sw':>5} {'/10k':>6}" for n in names)
    print(head)
    print("-" * len(head))
    for s in strata:
        line = f"{s:>11} {sites[s]:>9}"
        for n in names:
            pairs, sw = per[s][n]
            rate = f"{1e4 * sw / pairs:6.2f}" if pairs else f"{'n/a':>6}"
            line += f" | {'':>12} {pairs:>8} {sw:>5} {rate}"
        print(line)
    if second is None:
        return
    print()
    print(f"paired on common pairs: first = {first}, second = {second}; first/second = pairs "
          "where only that one is wrong")
    head = (f"{'stratum':>11} {'pairs':>8} {'first':>6} {'second':>6} {'both':>5} "
            f"{'sign p':>9} {'windows':>7} {'flip p':>9} {'exact':>5}")
    print(head)
    print("-" * len(head))
    for s in strata:
        q = paired[s]
        flip = permutation.window_sign_flip_test(q["windows"], RESAMPLES, SEED)
        print(f"{s:>11} {q['pairs']:>8} {q['first']:>6} {q['second']:>6} {q['both']:>5} "
              f"{sign_test(q['first'], q['first'] + q['second']):>9.2e} {flip.windows:>7} "
              f"{flip.p_value:>9.2e} {'yes' if flip.exact else 'no':>5}")


if __name__ == "__main__":
    main()
