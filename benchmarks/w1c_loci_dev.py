"""Development check for a W1c registration: do HG002 HiFi error loci predict errors elsewhere?

Bins of BIN bases holding either site of a wrong consecutive pair in `hprc_HG002` at 30x (either
phaser, v5.0q truth) form a list; on other development platforms at 30x, the wrong-pair rate of
pairs touching a listed bin is compared with the rest. No W1 or W1b test genome is read.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/w1c_loci_dev.py [BIN]`
"""
import sys
from collections import Counter
sys.path.insert(0, "benchmarks")
import confidence_score as cs
import genome_common as common
from stratified_gate import strata

CH = list(common.CONFIRMATORY)
BIN = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000
st = strata(CH)
print(f"W1c development: error-locus list from hprc_HG002 30x, both phasers, bins of {BIN} bp, "
      f"chromosomes {CH}; scored at 30x on hg001_v5, ped_NA12879, ped_NA12885, cal_v5", flush=True)
bl = set()
for ph in cs.PHASERS:
    rows = cs.load_rows("hprc_HG002", "30", ph, CH)
    w = [r for r in rows if r[3]]
    print("HG002 hifi", ph, "pairs", len(rows), "wrong", len(w), Counter(st(r[0], r[1], r[4]) for r in w).most_common())
    for c, a, k, wr, b, _ in w:
        bl.add((c, a // BIN)); bl.add((c, b // BIN))
print("blacklist bins", len(bl), "Mb", len(bl) * BIN / 1e6)
for plat in ("hg001_v5", "ped_NA12879", "ped_NA12885", "cal_v5"):
    for ph in cs.PHASERS:
        for d in ("30",):
            rows = cs.load_rows(plat, d, ph, CH)
            tin = [r for r in rows if (r[0], r[1] // BIN) in bl or (r[0], r[4] // BIN) in bl]
            win = sum(r[3] for r in tin); wall = sum(r[3] for r in rows)
            rin = win / max(1, len(tin)); rout = (wall - win) / max(1, len(rows) - len(tin))
            hard = [r for r in rows if st(r[0], r[1], r[4]) != "easy"]
            hin = [r for r in hard if (r[0], r[1] // BIN) in bl or (r[0], r[4] // BIN) in bl]
            print(f"{plat} {ph} {d}: pairs {len(rows)} touch {len(tin)/len(rows):.4f} wrong {wall} "
                  f"in {win} ({win/max(1,wall):.3f}) rate in {rin:.2e} out {rout:.2e} ratio {rin/max(rout,1e-12):.1f}; "
                  f"hard pairs touch {len(hin)/max(1,len(hard)):.3f} hard wrong caught {sum(r[3] for r in hin)}/{sum(r[3] for r in hard)}")
