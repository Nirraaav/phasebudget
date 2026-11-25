"""Exploratory: the library against `whatshap phase` on one chromosome: wMEC cost, unphased sites, tied support. Usage: explore_ties.py CHROM DEPTH."""
import gzip, json, sys
sys.path.insert(0, "benchmarks")
import genome_common as common
ex = common._module("extract_fragments.py"); sc = common._module("score_hg002.py")
from phasebudget.selection import read_weights
chrom, depth = sys.argv[1], sys.argv[2]
m, pos, _, _ = ex.load_fragments(str(common.fragments_path(chrom, depth, "realign")))
m = ex.phred_to_log_odds(m)
var = sc.variant_table(common.hetsnv_path(chrom)); truth = common.v5_truth(chrom, var)
pl = common._module("hg002_pipelines.py")
sel = pl.whatshap_selection(m, 15, pos).selected
r = json.load(gzip.open(common.chrom_dir(chrom) / "runs" / f"lib_{depth}_whatshap_c15.json.gz", "rt"))
lib = (r["alleles"], r["blocks"])
d = common.chrom_dir(chrom) / "external"
vcf = d / f"whatshap_{depth}.vcf"
if not vcf.exists(): vcf = d / f"whatshap_{depth}.vcf.gz"
tp = sc.whatshap_phasing(vcf, var)
tool = ([-1 if a is None else a for a in tp.alleles], [-1 if b is None else b for b in tp.blocks])
def cost(ph, reads):
    a, b = ph; total = 0.0; unph = 0
    for f in reads:
        groups = {}
        for k, al in enumerate(f.alleles):
            if al is None: continue
            j = f.start + k
            if a[j] < 0: unph += 1; continue
            g = groups.setdefault(b[j], [0.0, 0.0]); w = f.weights[k]
            g[0] += w * (al != a[j]); g[1] += w * (al == a[j])
        total += sum(min(x, y) for x, y in groups.values())
    return total, unph
def switches(ph):
    a, b = ph; n = w = 0
    for j in range(len(a) - 1):
        if a[j] >= 0 and a[j+1] >= 0 and b[j] >= 0 and b[j] == b[j+1] and truth[j] is not None and truth[j+1] is not None:
            n += 1; w += (a[j] == a[j+1]) != (truth[j] == truth[j+1])
    return w, n
reads = [m.fragments[k] for k in sel]
print(f"{chrom} {depth}x: selected {len(reads)} of {len(m.fragments)} reads")
for name, ph in (("library", lib), ("whatshap phase", tool)):
    c, u = cost(ph, reads); w, n = switches(ph)
    print(f"{name:>15}: MEC on selected reads {c:.1f} (alleles at unphased sites {u}); switches {w}/{n}")
# Equal terms: drop sites the tool leaves unphased from both; count blocks.
keep = [a >= 0 for a in tool[0]]
lib2 = ([x if k else -1 for x, k in zip(*lib[0:1], keep)], lib[1])
lib2 = ([x if k else -1 for x, k in zip(lib[0], keep)], lib[1])
for name, ph in (("library, tool sites", lib2), ("whatshap phase", tool)):
    c, u = cost(ph, reads)
    nb = len({b for a, b in zip(*ph) if a >= 0})
    print(f"{name:>20}: MEC {c:.1f}, blocks {nb}")
# Same block partition: score the library's alleles under the tool's blocks.
mixed = (lib2[0], tool[1])
print(f"{'library alleles, tool blocks':>20}: MEC {cost(mixed, reads)[0]:.1f}")
print("switches, library on tool sites:", switches(lib2), "; library on tool sites and tool blocks:", switches(mixed))
diff = sum(1 for x, y in zip(lib2[0], tool[0]) if x >= 0 and y >= 0)
print("sites phased by both:", diff)
# Per-site support under the library phasing: reads take the orientation that fits their block best.
a, b = lib
support = {}
for f in reads:
    groups = {}
    for k, al in enumerate(f.alleles):
        if al is None or a[f.start + k] < 0: continue
        j = f.start + k; g = groups.setdefault(b[j], [])
        g.append((j, al))
    for g in groups.values():
        agree = sum(al == a[j] for j, al in g); flip = agree < len(g) - agree
        for j, al in g:
            s = support.setdefault(j, [0, 0]); s[(al == a[j]) == flip] += 1   # [for, against]
from collections import Counter
dropped = [j for j in range(len(a)) if a[j] >= 0 and tool[0][j] < 0]
print("dropped sites:", len(dropped))
print("support [for, against] at dropped sites:", Counter(tuple(support.get(j, [0, 0])) for j in dropped).most_common(8))
kept = [j for j in range(len(a)) if a[j] >= 0 and tool[0][j] >= 0]
tie_kept = sum(1 for j in kept if support.get(j, [0, 0])[0] <= support.get(j, [0, 0])[1])
print(f"kept sites with for <= against: {tie_kept} of {len(kept)}")
wrong_at = lambda j: truth[j] is not None
