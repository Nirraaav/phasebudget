"""Post hoc, after scoring X8: what are the wrong trusted pairs of V2's two violating cells?

V2 failed in whatshap's `lowmap` group at beta 3e-3 (30x and 10x). For each
wrong pair the thresholds trust there, it asks whether LongPhase agrees with
whatshap, and whether Platinum Genomes 2017's NA12878 phase (an older,
Illumina-based pedigree truth) sides with whatshap or with the Platinum
Pedigree truth the test scored against. It also counts the violating units
per chromosome. It changes no verdict.

Usage: `PHASEBUDGET_PLATFORM=hg001_v5 PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/t3_lowmap_posthoc.py`
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import confidence_score as cs  # noqa: E402
import stratified_score as ss  # noqa: E402

from phasebudget.confidence import learn_then_test_by_group  # noqa: E402

BETA = 3e-3


def main() -> None:
    from stratified_gate import strata

    import genome_common as common

    assert common.PLATFORM == "hg001_v5", "run with PHASEBUDGET_PLATFORM=hg001_v5"
    chroms = list(common.CONFIRMATORY)
    stratum = strata(chroms)
    sc = common._module("score_hg002.py")
    print(f"post hoc (after X8's scoring): wrong trusted pairs in whatshap lowmap, beta {BETA:g}")
    for depth in ("30", "10"):
        cal = ss.grouped("cal_v5", depth, "whatshap", chroms, stratum)
        t = learn_then_test_by_group(cal, BETA, cs.DELTA)["lowmap"]
        ind = ss.indels("hg001_v5", chroms)
        rows = cs.load_rows("hg001_v5", depth, "whatshap", chroms)
        lp = {(c, a, b): w for c, a, k, w, b, _ in cs.load_rows("hg001_v5", depth, "longphase",
                                                                 chroms)}
        bad = [(c, a, k, b) for c, a, k, w, b, _ in rows
               if w and k >= t and (c, a) not in ind and (c, b) not in ind
               and stratum(c, a, b) == "lowmap"]
        platinum, chrom_units = {}, Counter()
        for chrom in sorted({c for c, *_ in bad}):
            variants = sc.variant_table(common.hetsnv_path(chrom))
            for v, p in zip(variants, common.platinum_truth(chrom, variants)):
                if p is not None:
                    platinum[(chrom, v[0])] = p
        pedigree = {}
        for chrom in sorted({c for c, *_ in bad}):
            variants = sc.variant_table(common.hetsnv_path(chrom))
            for v, p in zip(variants, common.pedigree_truth(chrom, variants)):
                if p is not None:
                    pedigree[(chrom, v[0])] = p
        tally = Counter()
        for c, a, k, b in bad:
            chrom_units[(c, a // cs.UNIT)] += 1
            same_lp = "LongPhase wrong too" if lp.get((c, a, b)) else (
                "LongPhase right" if (c, a, b) in lp else "LongPhase does not phase it")
            if (c, a) in platinum and (c, b) in platinum:
                plat_rel = platinum[(c, a)] == platinum[(c, b)]
                ped_rel = pedigree[(c, a)] == pedigree[(c, b)]
                side = "Platinum 2017 agrees with Pedigree" if plat_rel == ped_rel \
                    else "Platinum 2017 sides with whatshap"
            else:
                side = "Platinum 2017 does not phase it"
            tally[(same_lp, side)] += 1
        print(f"\n{depth}x: threshold {t}; {len(bad)} wrong trusted lowmap pairs in "
              f"{len(chrom_units)} units; read support (Wilson, reads) median "
              f"{sorted(k for _, _, k, _ in bad)[len(bad) // 2] if bad else '-'}")
        for (x, y), n in tally.most_common():
            print(f"  {n:>4}  {x}; {y}")
        per_chrom = Counter(c for c, _ in chrom_units)
        print("  units per chromosome: " + ", ".join(f"{c} {n}" for c, n in sorted(
            per_chrom.items(), key=lambda kv: -kv[1])))


if __name__ == "__main__":
    main()
