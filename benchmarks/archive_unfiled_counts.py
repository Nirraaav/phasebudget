"""Recompute, with a printed header, paper numbers that had no archived run file.

Each subcommand reads existing inputs and prints one table; nothing is written
to data/. Mock review 4, major issue 4, lists the numbers.

- `truth-audit`: on the 21 confirmatory chromosomes, how often the
  transmission truth and GIAB v5.0q agree per site, and what v5.0q says at
  the sites the hifiasm filter removed (paper, Data and truths; S5).
- `uncalled`: at one cap and depth, the variants inside a selected read's
  span that no selected read calls, and the boundaries the phasing joins
  with no read calling both variants, as the library reported them before
  commit `33b92bb` (S4). Selection and the DP are unchanged since that
  commit, so the pre-fix phasing is the current DP solution with blocks
  taken over covered variants.
- `gappy`: the share of each selection's reads with gap fraction at or above
  a threshold on one chromosome, depth and cap (paper, Mechanism).
- `arms`: the vote model's rate for each arm, the mean of f(i) over scored
  boundaries with i >= 1, against the observed rate, from an archived
  `benchmark_spanning.py --count informative` output (S9, Table arms).
- `depth`: the measured slice depth of one chromosome (paper, Data).

Usage, from the repository root:
`PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/archive_unfiled_counts.py truth-audit`
`PHASEBUDGET_PLATFORM=ont ... archive_unfiled_counts.py gappy --chrom chr20`
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402

from phasebudget.dp import CoverageDP  # noqa: E402
from phasebudget.selection import read_weights, select_optimal  # noqa: E402

pipelines = common._module("hg002_pipelines.py")
extractor = common._module("extract_fragments.py")
score = common._module("score_hg002.py")


def selections(chrom, depth, cap):
    matrix, positions, _, _ = extractor.load_fragments(
        str(common.fragments_path(chrom, depth, "realign")))
    matrix = extractor.phred_to_log_odds(matrix)
    return matrix, {
        "whatshap": pipelines.whatshap_selection(matrix, cap, positions).selected,
        "optimal": select_optimal(matrix, cap, read_weights(matrix)).selected,
    }


def audit_one(chrom):
    variants = score.variant_table(common.hetsnv_path(chrom))
    _, transmission, hifiasm = common.load_truths(chrom, [p - 1 for p, _, _ in variants])
    v5 = common.v5_truth(chrom, variants)
    both = agree = removed = removed_v5 = removed_agree = 0
    for t, h, v in zip(transmission, hifiasm, v5):
        if t is not None and v is not None:
            both += 1
            agree += t == v
        if t is not None and h is not None and t != h:
            removed += 1
            if v is not None:
                removed_v5 += 1
                removed_agree += t == v
    return chrom, (both, agree, removed, removed_v5, removed_agree)


def uncalled_one(job):
    chrom, depth, cap = job
    matrix, chosen = selections(chrom, depth, cap)
    truth = common.v5_truth(chrom, score.variant_table(common.hetsnv_path(chrom)))
    out = {}
    for name, selected in chosen.items():
        subset = matrix.subset(selected)
        haplotype = CoverageDP(subset).solve().haplotype
        coverage = subset.coverage()
        called = [False] * subset.n_variants
        both = [0] * (subset.n_variants - 1)
        for f in subset:
            for k, a in enumerate(f.alleles):
                if a is not None:
                    called[f.start + k] = True
                    if k + 1 < len(f.alleles) and f.alleles[k + 1] is not None:
                        both[f.start + k] += 1
        spanned = subset.spanned_boundaries()
        block, current, nxt = [None] * subset.n_variants, None, 0  # blocks() before 33b92bb
        for j in range(subset.n_variants):
            if coverage[j] == 0:
                current = None
                continue
            if current is None or j == 0 or not spanned[j - 1]:
                current, nxt = nxt, nxt + 1
            block[j] = current
        uncalled = [coverage[j] > 0 and not called[j] for j in range(subset.n_variants)]
        row = [sum(uncalled), sum(u and truth[j] is not None for j, u in enumerate(uncalled)),
               0, 0, 0, 0]
        for b in range(subset.n_variants - 1):
            if both[b] or block[b] is None or block[b] != block[b + 1]:
                continue
            if truth[b] is None or truth[b + 1] is None:
                continue
            wrong = (haplotype[b] ^ haplotype[b + 1]) != (truth[b] ^ truth[b + 1])
            k = 2 if uncalled[b] or uncalled[b + 1] else 4
            row[k] += 1
            row[k + 1] += wrong
        out[name] = row
    return out


def gappy_one(job):
    chrom, depth, cap, threshold = job
    matrix, chosen = selections(chrom, depth, cap)
    fragments = {f.read_id: f for f in matrix}
    out = {}
    for name, selected in chosen.items():
        gaps = [1 - sum(a is not None for a in fragments[r].alleles) / len(fragments[r].alleles)
                for r in selected]
        out[name] = (sum(g >= threshold for g in gaps), len(gaps))
    return out


def vote_f(i: int, q: float) -> float:
    total = 0.0
    for k in range(i + 1):
        p = math.comb(i, k) * q ** k * (1 - q) ** (i - k)
        total += p if 2 * k > i else p / 2 if 2 * k == i else 0.0
    return total


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("what", choices=("truth-audit", "uncalled", "gappy", "arms", "depth"))
    parser.add_argument("--chrom", default="chr20")
    parser.add_argument("--depth", default="full")
    parser.add_argument("--cap", type=int, default=3)
    parser.add_argument("--threshold", type=float, default=0.1)
    parser.add_argument("--spanning", help="`arms`: a benchmark_spanning.py output")
    parser.add_argument("--e", type=float, help="`arms`: per-allele error of the vote model")
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    head = f"archive_unfiled_counts.py {args.what}; platform {common.PLATFORM}"

    if args.what == "truth-audit":
        print(f"{head}; {len(common.CONFIRMATORY)} confirmatory chromosomes; truths: transmission "
              "(unfiltered), hifiasm phase transfer, GIAB v5.0q inside its benchmark BED")
        with ProcessPoolExecutor(args.jobs) as pool:
            rows = dict(pool.map(audit_one, common.CONFIRMATORY))
        total = [sum(r[i] for r in rows.values()) for i in range(5)]
        for chrom, r in list(rows.items()) + [("total", total)]:
            print(f"{chrom:>6}: both phase {r[0]:>8}, paternal = v5 first allele {r[1]:>8}; "
                  f"filter removed {r[2]:>5}, v5 phases {r[3]:>5}, agrees with transmission {r[4]:>5}")
    elif args.what == "uncalled":
        print(f"{head}; depth {args.depth}, cap {args.cap}, {len(common.CONFIRMATORY)} confirmatory "
              "chromosomes, truth GIAB v5.0q; phasing as reported before 33b92bb")
        jobs = [(c, args.depth, args.cap) for c in common.CONFIRMATORY]
        pooled: dict = {}
        with ProcessPoolExecutor(args.jobs) as pool:
            for out in pool.map(uncalled_one, jobs):
                for name, row in out.items():
                    acc = pooled.setdefault(name, [0] * 6)
                    for i, x in enumerate(row):
                        acc[i] += x
        print("selection: uncalled variants inside blocks (all, v5-phased); boundaries joined with "
              "no read calling both: touching an uncalled variant (scored, wrong), both called "
              "(scored, wrong)")
        for name, r in sorted(pooled.items(), reverse=True):
            print(f"{name:>8}: {r[0]:>5} {r[1]:>5} | {r[2]:>5} {r[3]:>5} | {r[4]:>5} {r[5]:>5}")
    elif args.what == "gappy":
        print(f"{head}; {args.chrom}, depth {args.depth}, cap {args.cap}; share of selected reads "
              f"with gap fraction >= {args.threshold}")
        out = gappy_one((args.chrom, args.depth, args.cap, args.threshold))
        for name, (k, n) in sorted(out.items(), reverse=True):
            print(f"{name:>8}: {k} of {n} = {100 * k / n:.2f} per cent")
    elif args.what == "arms":
        q = 2 * args.e * (1 - args.e)
        print(f"archive_unfiled_counts.py arms; from {args.spanning}; e {args.e}, q {q:.6f}; per 10,000 scored boundaries "
              "with i >= 1; model is the mean of f(i)")
        pooled = {}
        line = re.compile(r"^\s*(\w+)\s+(\d+)\s+(\w+)\s+(\d+)\s+\d+\s+[\d.]+\s+(\d+)\s+(\d+)\s")
        for text in open(args.spanning):
            m = line.match(text)
            if not m or int(m.group(4)) == 0:
                continue
            depth, cap, name, i, scored, wrong = m.groups()
            acc = pooled.setdefault((depth, int(cap), name), [0, 0, 0.0])
            acc[0] += int(scored)
            acc[1] += int(wrong)
            acc[2] += int(scored) * vote_f(int(i), q)
        for (depth, cap, name), (scored, wrong, model) in sorted(pooled.items()):
            obs, mod = 1e4 * wrong / scored, 1e4 * model / scored
            print(f"depth {depth:>4} cap {cap:>2} {name:>8}: observed {obs:6.2f}, model {mod:6.2f}, "
                  f"ratio {mod / obs:.2f}")
    else:
        record = json.loads((common.chrom_dir(args.chrom) / "depth.json").read_text())
        print(f"{head}; {args.chrom}; from "
              f"{(common.chrom_dir(args.chrom) / 'depth.json').relative_to(common.ROOT)}, "
              "written by genome_run.py")
        print(f"records {record['records']}, mean aligned length {record['mean_aligned_length']:.0f}, "
              f"chromosome length {record['length']}, depth {record['depth']:.2f}")


if __name__ == "__main__":
    main()
