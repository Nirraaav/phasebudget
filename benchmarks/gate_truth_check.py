"""Follow-up gate G-C (`docs/followup/error-controlled-phasing.md`): are the remaining HiFi errors real?

For every consecutive pair that stock `whatshap phase` (cap 15, the HiFi
tool-level run) gets wrong against GIAB v5.0q, this prints:
- the read support margin for the called relation (as in
  `gate_evidence_score.py`);
- whether the transmission truth agrees with v5.0q on that pair;
- whether HiPhase and LongPhase call the same relation as WhatsHap;
- the GIAB v3.6 strata the pair lies in.

A pair that every tool and most reads call one way, against v5.0q, points
at the truth rather than at the phaser. The analysis is post hoc and for
development only.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/gate_truth_check.py`
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402
from gate_evidence_score import pair_scores  # noqa: E402

extractor = common._module("extract_fragments.py")
score = common._module("score_hg002.py")
stratify = common._module("stratify.py")
TOOLS = common.ROOT / "data" / "tools"


def existing(path: Path) -> Path:
    gz = path.with_name(path.name + ".gz")
    return path if path.exists() or not gz.exists() else gz


def relation(phasing, j, k):
    a, b = phasing.alleles[j], phasing.alleles[k]
    if a is None or b is None or phasing.blocks[j] != phasing.blocks[k]:
        return None
    return a == b


def main() -> None:
    beds = {name: stratify.load_bed(path, common.CONFIRMATORY)
            for name, path in stratify.fetch().items()}
    print("gate G-C: stock whatshap phase cap 15 errors on HG002 HiFi, full depth, against v5.0q")
    print("chrom      pos1   margin support  transmission  hiphase  longphase  strata")
    totals = {"n": 0, "trans_agrees_v5": 0, "trans_disagrees": 0, "hiphase_same": 0,
              "longphase_same": 0, "margin_pos": 0}
    for chrom in common.CONFIRMATORY:
        variants = score.variant_table(common.hetsnv_path(chrom))
        positions = [p - 1 for p, _, _ in variants]
        v5 = common.v5_truth(chrom, variants)
        transmission = common.load_truths(chrom, positions)[1]
        d = TOOLS / chrom
        wh = score.whatshap_phasing(existing(d / "stock_c15.vcf"), variants)
        hp = score.whatshap_phasing(existing(d / "hiphase.vcf.gz"), variants)
        lp = score.whatshap_phasing(existing(d / "longphase.vcf"), variants)
        matrix, _, _, _ = extractor.load_fragments(
            str(common.fragments_path(chrom, "full", "realign")))
        rows = pair_scores(matrix.fragments, wh.alleles, wh.blocks, v5)
        phased = [j for j, a in enumerate(wh.alleles) if a is not None and v5[j] is not None]
        nxt = [(a, b) for a, b in zip(phased, phased[1:]) if wh.blocks[a] == wh.blocks[b]]
        member = {s: stratify.inside(beds[s][chrom], positions) for s in beds}
        for (j, k), (margin, support, wrong) in zip(nxt, rows):
            if not wrong:
                continue
            called = relation(wh, j, k)
            t = None if transmission[j] is None or transmission[k] is None else \
                (transmission[j] == transmission[k]) == (v5[j] == v5[k])
            h, lo = relation(hp, j, k), relation(lp, j, k)
            strata = ",".join(s for s in beds if member[s][j] and member[s][k]) or "-"
            fmt = lambda x: "n/a" if x is None else ("same" if x == called else "differs")
            print(f"{chrom:<6} {variants[j][0]:>10} {margin:>6} {support:>7}  "
                  f"{'n/a' if t is None else ('agrees' if t else 'DIFFERS'):>12}  "
                  f"{fmt(h):>7}  {fmt(lo):>9}  {strata}")
            totals["n"] += 1
            totals["trans_agrees_v5"] += t is True
            totals["trans_disagrees"] += t is False
            totals["hiphase_same"] += h == called
            totals["longphase_same"] += lo == called
            totals["margin_pos"] += margin > 0
    print(f"\nwrong pairs {totals['n']}: transmission agrees with v5.0q at "
          f"{totals['trans_agrees_v5']}, differs at {totals['trans_disagrees']}; HiPhase calls the "
          f"same (wrong) relation at {totals['hiphase_same']}, LongPhase at "
          f"{totals['longphase_same']}; read margin above 0 for the called relation at "
          f"{totals['margin_pos']}")


if __name__ == "__main__":
    main()
