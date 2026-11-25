"""Exploratory: allele error of reads by gap fraction and by number of calls (one chromosome)."""
import sys
sys.path.insert(0, "benchmarks")
import genome_common as common
ex = common._module("extract_fragments.py"); sc = common._module("score_hg002.py")
chrom = sys.argv[1]
m, pos, _, _ = ex.load_fragments(str(common.fragments_path(chrom, "full", "realign")))
var = sc.variant_table(common.hetsnv_path(chrom))
truth = common.v5_truth(chrom, var) if common.SAMPLE == "HG002" else common.load_truths(chrom, pos)[1]
rows = []
for f in m.fragments:
    calls = [(f.start + k, a) for k, a in enumerate(f.alleles) if a is not None and truth[f.start + k] is not None]
    span = len(f.alleles); called = sum(a is not None for a in f.alleles)
    if len(calls) < 3: continue
    agree = sum(a == truth[j] for j, a in calls)
    if 2 * agree == len(calls): continue
    err = min(agree, len(calls) - agree)
    rows.append((1 - called / span, len(calls), err))
tot = sum(r[1] for r in rows); E = sum(r[2] for r in rows)
print(f"{chrom} {common.PLATFORM}: reads {len(rows)}, e {E/tot:.5f}")
zero = sum(1 for r in rows if r[2] == 0)
import math
# expected reads with zero errors under independence
exp0 = sum((1 - E/tot) ** r[1] for r in rows)
print(f"reads with no error: observed {zero}, expected under independence {exp0:.0f}")
top = sorted(rows, key=lambda r: -r[2] / r[1])
k = len(rows) // 20
print(f"worst 5% of reads hold {sum(r[2] for r in top[:k]) / E:.2f} of errors")
for lo, hi in [(0, .02), (.02, .05), (.05, .1), (.1, .2), (.2, .4), (.4, 1.01)]:
    sel = [r for r in rows if lo <= r[0] < hi]
    if sel:
        print(f"gap fraction [{lo:.2f},{hi:.2f}): reads {len(sel):>6}, e {sum(r[2] for r in sel)/sum(r[1] for r in sel):.5f}")
for lo, hi in [(3, 10), (10, 30), (30, 100), (100, 300), (300, 10**6)]:
    sel = [r for r in rows if lo <= r[1] < hi]
    if sel:
        print(f"calls [{lo},{hi}): reads {len(sel):>6}, e {sum(r[2] for r in sel)/sum(r[1] for r in sel):.5f}")
