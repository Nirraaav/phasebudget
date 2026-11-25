"""Round 8, L3 (development): parental labels decide disputes between a phaser and its key.

`docs/design/l3-outside-adjudication.md` states the computation and the rules. Per genome
(HG002 `q3_q100` against Q100, HG005 `k3_hg005` against the HPRC release 2 dipcall), phaser,
side of the GIAB BED and group, it counts pairs by whether the phaser disputes the key and
which side the parents' `agree` label takes. Rows are Q3's and K3b's caches.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/l3_outside_adjudication.py`
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

PHASERS = ("whatshap", "longphase")
GROUPS = ("indel", "segdup", "lowmap", "tandem", "difficult", "easy")
CHROMS = tuple(f"chr{i}" for i in range(1, 23) if i != 20)
GENOMES = {  # name: (platform, rows path pattern, BED, key)
    "HG002": ("q3_q100", "data/r6_q3/rows/{phaser}_30_{chrom}.json",
              "data/giab_v5/HG002_GRCh38_v5.0q_smvar.benchmark.bed", "Q100 dipcall"),
    "HG005": ("k3_hg005", "data/k2/rows/k3_hg005_{phaser}_30_{chrom}.json",
              "data/giab_hg005/HG005_GRCh38_1_22_v4.2.1_benchmark.bed",
              "HPRC release 2 dipcall"),
}
A_IN_MIN, P_MIN, MIN_LABELLED = 0.9, 0.5, 20


# ---- pure logic, tested in tests/test_l3_outside_adjudication.py ------------------

def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """95% Wilson interval of `k / n`; (0, 1) when `n` is 0."""
    if n == 0:
        return 0.0, 1.0
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(0.0, centre - half), min(1.0, centre + half)


def corrected(f: float, e: float) -> float:
    """Phaser-error share `p` from `f = p(1 - e) + (1 - p)e`, clipped to [0, 1]."""
    if e >= 0.5:
        return float("nan")
    return min(1.0, max(0.0, (f - e) / (1 - 2 * e)))


def phaser_share(k_key: int, k_ph: int, e_wrong: int, e_n: int) -> tuple[float, float, float]:
    """`(p, lower, upper)`: corrected phaser share of disputes and its 95% interval, from
    `k_key` labels siding with the key, `k_ph` with the phaser, and `e_wrong` of `e_n`
    labelled undisputed pairs whose label differs from both."""
    n = k_key + k_ph
    if n == 0:
        return float("nan"), 0.0, 1.0
    e = e_wrong / e_n if e_n else 0.0
    f_lo, f_hi = wilson(k_key, n)
    e_lo, e_hi = wilson(e_wrong, e_n) if e_n else (0.0, 0.0)
    # the low end takes whichever bound of e lowers p (e_lo when f > 1/2), the high end the other
    lows = [corrected(f_lo, x) for x in (e_lo, e_hi)]
    highs = [corrected(f_hi, x) for x in (e_lo, e_hi)]
    return corrected(k_key / n, e), min(lows), max(highs)


# ---- counts in a platform's environment -------------------------------------------

def counts(name: str) -> dict:
    """`{"phaser|side|group|disputed|label": pairs}`, label one of none, key, phaser."""
    import genome_common as common
    import q3_outside_benchmark as q3
    import w5_called as wc
    import w5_noise as wn
    from stratified_gate import strata

    platform, pattern, bed, _ = GENOMES[name]
    assert common.PLATFORM == platform and common.SAMPLE == name
    sc = common._module("score_hg002.py")
    chroms = list(CHROMS)
    stratum = strata(chroms)
    inside = q3.bed_lookup(ROOT / bed, chroms)
    parents = {r: (wc.calls(fa, set(chroms)), wc.calls(mo, set(chroms)))
               for r, (fa, mo) in wc.REPLICATES.items()}
    out: Counter = Counter()
    for chrom in chroms:
        variants = sc.variant_table(common.hetsnv_path(chrom))
        key = {v[0]: t for v, t in zip(variants, common.primary_truth(chrom, variants))
               if t is not None}
        per = [wn.labelled(chrom, variants, *parents[r]) for r in wc.REPLICATES]
        label = {}
        for pos1, _, _ in variants:
            kept = wn.FILTERS["agree"](*(p.get(pos1) for p in per))
            if kept:
                label[pos1] = kept[0]
        for phaser in PHASERS:
            rows = json.loads((ROOT / pattern.format(phaser=phaser, chrom=chrom)).read_text())
            for p1, p2, _, _, _, wrong, indel, _ in rows:
                key_rel = key[p1] == key[p2]
                side = "inside" if inside(chrom, p1) and inside(chrom, p2) else "outside"
                group = "indel" if indel else stratum(chrom, p1, p2)
                if p1 in label and p2 in label:
                    lab_rel = label[p1] == label[p2]
                    verdict = "key" if lab_rel == key_rel else "phaser"
                else:
                    verdict = "none"
                out[f"{phaser}|{side}|{group}|{int(bool(wrong))}|{verdict}"] += 1
        print(f"  {name} {chrom} done", flush=True)
    return dict(out)


def load_counts(name: str) -> dict:
    path = ROOT / "data" / "l3" / f"{name}_counts.json"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_name(path.name + ".part")
        env = dict(os.environ, PHASEBUDGET_PLATFORM=GENOMES[name][0],
                   PYTHONPATH=f"{ROOT / 'src'}:{HERE}")
        subprocess.run([sys.executable, __file__, "--counts", name, str(part)], env=env,
                       check=True)
        part.replace(path)
    return json.loads(path.read_text())


# ---- report -----------------------------------------------------------------------

def tally(c: dict, phaser: str, side: str, groups) -> dict:
    t = Counter()
    for k, n in c.items():
        ph, sd, g, disputed, verdict = k.split("|")
        if ph == phaser and sd == side and g in groups:
            t[(int(disputed), verdict)] += n
    return t


def line(t: dict) -> tuple[str, tuple]:
    pairs = sum(t.values())
    disputed = sum(n for (d, _), n in t.items() if d)
    k_key, k_ph = t[(1, "key")], t[(1, "phaser")]
    e_wrong, e_n = t[(0, "phaser")], t[(0, "key")] + t[(0, "phaser")]
    p, lo, hi = phaser_share(k_key, k_ph, e_wrong, e_n)
    text = (f"pairs {pairs:>8} disputed {disputed:>6} labelled {e_n / max(1, pairs - disputed):.3f}"
            f"/{(k_key + k_ph) / max(1, disputed):.3f} e {e_wrong}/{e_n} = "
            f"{e_wrong / max(1, e_n):.2e} verdicts key {k_key:>5} phaser {k_ph:>5} "
            f"p {p:.3f} [{lo:.3f}, {hi:.3f}]")
    return text, (disputed, k_key, k_ph, p, lo, hi)


def main() -> None:
    if len(sys.argv) > 1 and sys.argv[1] == "--counts":
        Path(sys.argv[3]).write_text(json.dumps(counts(sys.argv[2])))
        return
    started = time.time()
    print("L3 (development, docs/design/l3-outside-adjudication.md): 30x, whatshap 2.8 and "
          "LongPhase 2.0.2 as stored; " + "; ".join(
              f"{g}: platform {p}, key {k}, inside = both variants in {Path(b).name}"
              for g, (p, _, b, k) in GENOMES.items())
          + "; labels: parents' ONT wf-human-variation sup calls, both flowcell replicates, "
          "`agree` filter, paternal allele by w5_called.paternal; pairs = Q3 and K3b cached "
          "consecutive pairs the key phases; groups indel first, then GIAB v3.6 stratum; "
          f"chromosomes {' '.join(CHROMS)}; p = (f - e)/(1 - 2e), 95% Wilson bounds with e at "
          f"its own bound; L3-0 a_in >= {A_IN_MIN}; L3-1/L3-2 lower bound of p >= {P_MIN}; "
          f"reweighting counts groups with < {MIN_LABELLED} labelled disputes at 0", flush=True)
    with ThreadPoolExecutor(2) as pool:
        data = dict(zip(GENOMES, pool.map(load_counts, GENOMES)))
    print("\ncolumns: pairs, disputed (phaser differs from key), labelled share of undisputed / "
          "disputed, e = label differs from both where phaser and key agree, verdicts on "
          "labelled disputes, corrected phaser-error share p of disputes [95%]")
    results = {}
    for name, c in data.items():
        for phaser in PHASERS:
            for side in ("inside", "outside"):
                print(f"\n{name} {phaser} {side}")
                for label, groups in (("all", GROUPS), ("SNV only", GROUPS[1:]),
                                      *((g, (g,)) for g in GROUPS)):
                    text, r = line(tally(c, phaser, side, groups))
                    results[(name, phaser, side, label)] = r
                    print(f"  {label:<9} {text}")
                # worst-case reweighting by group
                total = sum(results[(name, phaser, side, g)][0] for g in GROUPS)
                worst = sum(r[0] * (r[3] if r[1] + r[2] >= MIN_LABELLED else 0.0)
                            for g in GROUPS for r in [results[(name, phaser, side, g)]])
                print(f"  reweighted by group, unlabelled groups at 0: {worst / max(1, total):.3f}")
    print()
    ok0 = True
    for phaser in PHASERS:
        _, k_key, k_ph, *_ = results[("HG002", phaser, "inside", "all")]
        a_in = k_key / max(1, k_key + k_ph)
        lo, hi = wilson(k_key, k_key + k_ph)
        ok0 &= a_in >= A_IN_MIN
        print(f"L3-0 {phaser}: a_in {k_key}/{k_key + k_ph} = {a_in:.4f} [{lo:.4f}, {hi:.4f}]")
    print(f"L3-0: {'holds' if ok0 else 'fails'}")
    for rule, name in (("L3-1", "HG005"), ("L3-2", "HG002")):
        ok = True
        for phaser in PHASERS:
            _, _, _, p, lo, hi = results[(name, phaser, "outside", "all")]
            ok &= lo >= P_MIN
            print(f"{rule} {name} {phaser} outside: p {p:.3f} [{lo:.3f}, {hi:.3f}]")
        verdict = ("not decidable (L3-0 fails)" if not ok0 else
                   "supported" if ok else "not supported")
        print(f"{rule}: {verdict}")
    print(f"runtime {time.time() - started:.0f} s")


if __name__ == "__main__":
    main()
