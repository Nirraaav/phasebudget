"""Gates of Z1 and Z2, the October directions memo (development repository) (round 2), on data already seen.

Z1: per phaser, the switch error on consecutive pairs, the gene-weighted
relation error (track D's population, genes with at least 100 relations),
and the split of wrong links into flips (two wrong links around one variant)
and long switches. Z2: on consecutive pairs both phasers phase, the error
among pairs where the two phasers agree, against each phaser alone.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/phaser_error_types.py`
"""

from __future__ import annotations

import bisect
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import chain_relations as cr  # noqa: E402
import confidence_score as cs  # noqa: E402
import genome_common as common  # noqa: E402

from phasebudget.confidence import chain_counts  # noqa: E402

MIN_RELATIONS = 100


def links(rows):
    """Wrong flags of consecutive labelled variants in one block, as runs per block."""
    for k in range(len(rows) - 1):
        a, b = rows[k], rows[k + 1]
        if a[1] == b[1] and a[3] is not None and a[2] is not None and b[2] is not None:
            yield k, a[2] != b[2]


def main() -> None:
    chroms = list(common.CONFIRMATORY)
    print(f"Z1 and Z2 gates (development; data seen): {cr.DEPTH}x; chains of track D; genes "
          f"{cr.GENCODE.name} protein-coding with >= {MIN_RELATIONS} relations; flips = two wrong "
          "links around one variant; chromosomes " + " ".join(chroms))
    gene_map = cr.genes(chroms)
    for platform in cr.DATASETS:
        print(f"\n{platform}")
        print("  phaser      links      wrong  switch/10k  flips  long  long share  "
              "gene-weighted relation error  genes")
        for phaser in cs.PHASERS:
            chains = cr.load_chains(platform, phaser, chroms)
            n_links = wrong = flips = 0
            shares = []
            for chrom, rows in chains.items():
                bad = dict(links(rows))
                n_links += len(bad)
                wrong += sum(bad.values())
                flips += sum(2 for k in bad if bad[k] and bad.get(k + 1))
                pos = [r[0] for r in rows]
                for start, end, _ in gene_map[chrom]:
                    lo, hi = bisect.bisect_left(pos, start), bisect.bisect_right(pos, end)
                    sub = rows[lo:hi]
                    total = [0, 0]
                    start_run = 0
                    for k in range(1, len(sub) + 1):  # maximal runs of one block
                        if k == len(sub) or sub[k][1] != sub[k - 1][1] or sub[k - 1][3] is None:
                            run = sub[start_run:k]
                            if len(run) > 1:
                                c = chain_counts([True] * (len(run) - 1), [r[2] for r in run],
                                                 [0] * len(run))
                                total[0] += c.pairs
                                total[1] += c.wrong
                            start_run = k
                    if total[0] >= MIN_RELATIONS:
                        shares.append(total[1] / total[0])
            long = wrong - flips
            print(f"  {phaser:<10} {n_links:>9}  {wrong:>6}  {1e4 * wrong / n_links:>10.2f}  "
                  f"{flips:>5}  {long:>5}  {long / max(1, wrong):>10.3f}  "
                  f"{sum(shares) / len(shares):>28.2e}  {len(shares):>5}")
        rows = {ph: {(c, p1, p2): w for c, p1, k, w, p2, _ in cs.load_rows(platform, cr.DEPTH, ph,
                                                                        chroms)}
                for ph in cs.PHASERS}
        shared = rows[cs.PHASERS[0]].keys() & rows[cs.PHASERS[1]].keys()
        agree = [key for key in shared if rows["whatshap"][key] == rows["longphase"][key]]
        err = {ph: sum(rows[ph][key] for key in shared) / len(shared) for ph in cs.PHASERS}
        agree_err = sum(rows["whatshap"][key] for key in agree) / len(agree)
        best = min(err.values())
        print(f"  Z2: shared pairs {len(shared)}, agree {len(agree)} ({len(agree) / len(shared):.4f}); "
              f"error whatshap {err['whatshap']:.2e}, LongPhase {err['longphase']:.2e}, agreed "
              f"{agree_err:.2e} ({agree_err / best:.3f} of the better)")


if __name__ == "__main__":
    main()
