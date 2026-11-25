"""Round 5, R1 (development): error rates of phase truths without a gold standard.

`docs/design/r1-truth-audit.md` states the model, the raters, the fits and
the gate. Per genome, every rater's phase is read once per chromosome (in a
subprocess with that genome's platform) and cached; the fits run on counts
of agreement patterns per group.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/r1_truth_audit.py`
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

ROOT = HERE.parent
CACHE = ROOT / "data" / "confidence_rows"
CHROMS = tuple(f"chr{i}" for i in range(1, 23) if i != 20)
# Per genome (its platform): rater -> ("truth", loader name) or ("phaser", platform of the reads).
RATERS = {
    "hg001_v5": {"PED": ("truth", "pedigree"), "PG17": ("truth", "platinum"),
                 "KG": ("truth", "kg"), "ONT": ("phaser", "hg001_v5"),
                 "HIFI": ("phaser", "phifi_NA12878")},
    "cal_v5": {"V5": ("truth", "v5"), "TRIO": ("truth", "trio"), "HIFIASM": ("truth", "hifiasm"),
               "ONT": ("phaser", "cal_v5"), "HIFI": ("phaser", "hprc_HG002")},
}
GENOME = {"hg001_v5": "NA12878", "cal_v5": "HG002"}
MAIN = {"hg001_v5": ("PED", "PG17", "KG", "ONT"), "cal_v5": ("V5", "TRIO", "HIFIASM", "ONT")}
CHECK = "HIFI"
HARD = ("lowmap", "segdup", "tandem")
BLOCK = 1_000_000
BOOTSTRAPS = 100


def fit(counts: dict, iterations: int = 5000, tol: float = 1e-13) -> tuple:
    """Two-class latent model on `{pattern: count}`, a pattern a tuple of 0/1 (1 = the rater says
    the two variants are on the same haplotype). Returns `(pi, errors, loglik)`: `pi` the share
    of relations that are `same`, `errors[r]` the chance rater `r` reports the wrong relation."""
    patterns = [p for p, n in counts.items() if n]
    if not patterns:
        raise ValueError("no relations")
    m = len(patterns[0])
    total = sum(counts[p] for p in patterns)
    pi, err = 0.5, [0.05] * m
    previous = -math.inf
    for _ in range(iterations):
        post = {p: posterior(p, pi, err) for p in patterns}
        pi = sum(counts[p] * post[p] for p in patterns) / total
        err = [sum(counts[p] * (post[p] * (p[r] != 1) + (1 - post[p]) * (p[r] != 0))
                   for p in patterns) / total for r in range(m)]
        err = [min(max(e, 1e-12), 0.5 - 1e-9) for e in err]
        loglik = sum(counts[p] * math.log(likelihood(p, pi, err)) for p in patterns)
        if loglik - previous < tol * total:
            break
        previous = loglik
    return pi, tuple(err), loglik


def likelihood(pattern: tuple, pi: float, err) -> float:
    same = diff = 1.0
    for r, x in enumerate(pattern):
        same *= (1 - err[r]) if x == 1 else err[r]
        diff *= (1 - err[r]) if x == 0 else err[r]
    return pi * same + (1 - pi) * diff


def posterior(pattern: tuple, pi: float, err) -> float:
    """P(the relation is `same` | pattern)."""
    same = diff = 1.0
    for r, x in enumerate(pattern):
        same *= (1 - err[r]) if x == 1 else err[r]
        diff *= (1 - err[r]) if x == 0 else err[r]
    return pi * same / (pi * same + (1 - pi) * diff)


def g2(counts: dict, pi: float, err) -> tuple[float, int]:
    """Likelihood-ratio statistic of the observed pattern counts against the fit, and its degrees
    of freedom (patterns minus one minus parameters)."""
    total = sum(counts.values())
    stat = sum(2 * n * math.log(n / (total * likelihood(p, pi, err))) for p, n in counts.items()
               if n)
    m = len(err)
    return stat, 2 ** m - 1 - (m + 1)


def pair_outliers(counts: dict, pi: float, err) -> list:
    """Per pair of raters: relations where those two agree with each other and disagree with
    every other rater, observed and expected under the fit."""
    total = sum(counts.values())
    m = len(err)
    out = []
    for a in range(m):
        for b in range(a + 1, m):
            observed = expected = 0.0
            for z in (0, 1):
                p = tuple(z if r in (a, b) else 1 - z for r in range(m))
                observed += counts.get(p, 0)
                expected += total * likelihood(p, pi, err)
            out.append((a, b, int(observed), expected))
    return out


def read_table(platform: str, chrom: str) -> dict:
    """`{"variants": [(pos1, ref, alt)], rater: [value per variant]}` in `platform`'s environment;
    a truth's value is its first allele, a phaser's `[allele, block]`, None where it has none."""
    import pysam

    import genome_common as common
    from phase_confidence import phasing_from_vcf

    sc = common._module("score_hg002.py")
    variants = sc.variant_table(common.hetsnv_path(chrom))
    positions = [p - 1 for p, _, _ in variants]

    def patmat(path):
        index = {v: j for j, v in enumerate(variants)}
        out = [None] * len(variants)
        with pysam.TabixFile(str(path)) as tabix:
            for line in tabix.fetch(chrom):
                f = line.rstrip("\n").split("\t")
                j = index.get((int(f[1]), f[3], f[4]))
                fields = dict(zip(f[8].split(":"), f[9].split(":")))
                if j is not None and fields.get("GT") in ("0|1", "1|0") \
                        and fields.get("PS") == "PATMAT":
                    out[j] = int(fields["GT"][0])
        return out

    loaders = {
        "pedigree": lambda: common.pedigree_truth(chrom, variants),
        "platinum": lambda: common.platinum_truth(chrom, variants),
        "kg": lambda: common.kg_transmission(chrom, variants, "NA12878"),
        "v5": lambda: common.v5_truth(chrom, variants),
        # The per-chromosome trio files have no tabix index: read them whole.
        "trio": lambda: sc.parental_truth(common.GIAB / "truth" / f"HG002_{chrom}_trio_truth.vcf.gz",
                                          variants),
        "hifiasm": lambda: patmat(common.GIAB / "HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11"
                                  "_phasetransfer.vcf.gz"),
    }
    table = {"variants": [list(v) for v in variants]}
    for name, (kind, source) in RATERS[platform].items():
        if kind == "truth":
            table[name] = loaders[source]()
            continue
        vcf = common.ROOT / "data" / f"ext_{source}" / chrom / "external" / "whatshap_30.vcf.gz"
        if not vcf.exists():
            table[name] = None
            continue
        phasing = phasing_from_vcf(str(vcf), common.VCF_SAMPLE, chrom, positions)
        table[name] = [None if a is None else [a, b]
                       for a, b in zip(phasing.alleles, phasing.blocks)]
    return table


def load_table(platform: str, chrom: str) -> dict:
    path = CACHE / f"r1_{platform}_{chrom}.json"
    # A missing HiFi phasing is not cached as missing for good: retry once its run has finished.
    if path.exists():
        table = json.loads(path.read_text())
        if all(table.get(r) is not None for r in RATERS[platform]):
            return table
        source = RATERS[platform][CHECK][1]
        if not (ROOT / "data" / f"ext_{source}" / chrom / "external" / "whatshap_30.vcf.gz"
                ).exists():
            return table
    env = dict(os.environ, PHASEBUDGET_PLATFORM=platform, PYTHONPATH=f"{ROOT / 'src'}:{HERE}")
    subprocess.run([sys.executable, __file__, "--table", platform, chrom, str(path)], env=env,
                   check=True)
    return json.loads(path.read_text())


def relation(value_a, value_b) -> int | None:
    """1 if the rater puts the two variants on the same haplotype, 0 if not, None if it does not
    relate them (a value missing, or a phaser's two blocks differ)."""
    if value_a is None or value_b is None:
        return None
    if isinstance(value_a, list):
        if value_a[1] is None or value_a[1] != value_b[1]:
            return None
        return int(value_a[0] == value_b[0])
    return int(value_a == value_b)


def relations(tables: dict, raters, stratum, extra=None) -> list:
    """`(chrom, pos1 a, pos1 b, group, block, pattern, extra relation)` for consecutive variants
    that every rater in `raters` relates; `extra` is one more rater read but not required."""
    out = []
    for chrom, table in tables.items():
        variants = table["variants"]
        keep = [j for j in range(len(variants)) if all(table[r] is not None and
                                                       table[r][j] is not None for r in raters)]
        for i, j in zip(keep, keep[1:]):
            pattern = tuple(relation(table[r][i], table[r][j]) for r in raters)
            if None in pattern:
                continue
            (a, ra, aa), (b, rb, ab) = variants[i], variants[j]
            indel = any(len(x) != 1 for x in (ra, aa, rb, ab))
            group = "indel" if indel else stratum(chrom, a, b)
            more = None
            if extra is not None and table.get(extra) is not None:
                more = relation(table[extra][i], table[extra][j])
            out.append((chrom, a, b, group, (chrom, a // BLOCK), pattern, more))
    return out


def counts_of(rows) -> Counter:
    return Counter(r[5] for r in rows)


def bootstrap(rows, seed: int = 0) -> list:
    """`BOOTSTRAPS` fits on 1 Mb blocks resampled with replacement."""
    blocks = defaultdict(Counter)
    for r in rows:
        blocks[r[4]][r[5]] += 1
    keys = sorted(blocks)
    rng = random.Random(seed)
    fits = []
    for _ in range(BOOTSTRAPS):
        total = Counter()
        for k in rng.choices(keys, k=len(keys)):
            total.update(blocks[k])
        fits.append(fit(total)[1])
    return fits


def quantile(values, q: float) -> float:
    values = sorted(values)
    return values[min(len(values) - 1, max(0, int(round(q * (len(values) - 1)))))]


def audit(platform: str, stratum) -> tuple[list, dict]:
    """Report lines and gate parts for one genome."""
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(6) as pool:
        tables = dict(zip(CHROMS, pool.map(lambda c: load_table(platform, c), CHROMS)))
    raters = MAIN[platform]
    have_check = all(t.get(CHECK) is not None for t in tables.values())
    lines = [f"\n{GENOME[platform]} ({platform}): main fit raters {', '.join(raters)}; check "
             f"rater {CHECK} {'present' if have_check else 'NOT YET AVAILABLE'}"]
    rows = relations(tables, raters, stratum, CHECK if have_check else None)
    groups = sorted({r[3] for r in rows})
    gate = {"stable": True, "found": [], "check": [0, 0]}
    lines.append("group       relations  pi      " + "  ".join(f"{r:>24}" for r in raters)
                 + "   G2 (df)  odd/even")
    for g in groups:
        mine = [r for r in rows if r[3] == g]
        counts = counts_of(mine)
        pi, err, _ = fit(counts)
        boots = bootstrap(mine)
        odd = fit(counts_of([r for r in mine if int(r[0][3:]) % 2]))[1]
        even = fit(counts_of([r for r in mine if not int(r[0][3:]) % 2]))[1]
        stat, df = g2(counts, pi, err)
        cells = []
        for k, name in enumerate(raters):
            lo, hi = quantile([b[k] for b in boots], 0.025), quantile([b[k] for b in boots], 0.975)
            cells.append(f"{err[k]:.2e} [{lo:.1e},{hi:.1e}]")
            truth = RATERS[platform][name][0] == "truth"
            if truth and err[k] >= 1e-4:
                if max(odd[k], even[k]) > 2 * max(min(odd[k], even[k]), 1e-12):
                    gate["stable"] = False
            if truth and g in HARD and err[k] >= 1e-3 and lo >= 1e-4:
                gate["found"].append((name, g, err[k], lo))
        ratio = "; ".join(f"{raters[k]} {odd[k]:.1e}/{even[k]:.1e}" for k in range(len(raters)))
        lines.append(f"{g:<10} {len(mine):>10}  {pi:.3f}  " + "  ".join(f"{c:>24}" for c in cells)
                     + f"   {stat:.1f} ({df})  {ratio}")
        for a, b, observed, expected in pair_outliers(counts, pi, err):
            if observed > 2 * expected + 5:
                lines.append(f"      {raters[a]} and {raters[b]} agree against the rest "
                             f"{observed} times, {expected:.1f} expected under independence")
        if have_check and g in HARD:
            for r in mine:
                p_same = posterior(r[5], pi, err)
                call = 1 if p_same >= 0.5 else 0
                wrong_truth = [name for k, name in enumerate(raters)
                               if RATERS[platform][name][0] == "truth"
                               and (p_same if r[5][k] == 0 else 1 - p_same) >= 0.99]
                if wrong_truth and r[6] is not None:
                    gate["check"][0] += 1
                    gate["check"][1] += r[6] == call
    lines.append("leave-one-out (group: rater error rates without the named rater)")
    for drop in raters:
        rest = tuple(r for r in raters if r != drop)
        sub = relations(tables, rest, stratum)
        parts = []
        for g in groups:
            mine = counts_of([r for r in sub if r[3] == g])
            if mine:
                err = fit(mine)[1]
                parts.append(f"{g} " + ",".join(f"{n}={e:.1e}" for n, e in zip(rest, err)))
        lines.append(f"  without {drop}: " + "; ".join(parts))
    return lines, gate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", nargs=3, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.table:
        platform, chrom, path = args.table
        Path(path + ".part").write_text(json.dumps(read_table(platform, chrom)))
        Path(path + ".part").rename(path)
        return
    from stratified_gate import strata

    print(f"R1 (development, docs/design/r1-truth-audit.md): latent class fit per group of the "
          f"relation between consecutive variants every rater relates; chromosomes "
          f"{' '.join(CHROMS)}; {BOOTSTRAPS} bootstraps of {BLOCK // 1_000_000} Mb blocks; error "
          "rates with 95% intervals; G2 against the fit; odd/even fits", flush=True)
    stratum = strata(list(CHROMS))
    gates = {}
    for platform in RATERS:
        lines, gates[platform] = audit(platform, stratum)
        print("\n".join(lines), flush=True)
    gate = gates["hg001_v5"]
    checked, agreed = gate["check"]
    c_ok = checked >= 20 and agreed >= 0.75 * checked
    print(f"\nR1 gate (NA12878): (a) stable {'ok' if gate['stable'] else 'FAIL'}; (b) truth error "
          f">= 1e-3 in a hard group {'ok' if gate['found'] else 'FAIL'} "
          f"({', '.join(f'{n} {g} {e:.1e} (2.5% {lo:.1e})' for n, g, e, lo in gate['found'])}); "
          f"(c) HiFi agrees with the fit on {agreed} of {checked} flagged relations "
          f"{'ok' if c_ok else ('PENDING' if checked == 0 else 'FAIL')}")
    print(f"R1 gate: {'passed' if gate['stable'] and gate['found'] and c_ok else 'not passed'}"
          + (" (part c waits for phifi_NA12878)" if checked == 0 else ""))


if __name__ == "__main__":
    main()
