"""Round 4, W2 gate (development): parent of origin from imprinted DMRs, with and without certificates.

Long-read methylation can tell which haplotype is maternal: at a known
imprinted DMR one parent's allele is methylated. Tools that assign parent of
origin this way (Akbari et al., listed in the October directions memo (development repository))
carry the label from the DMR to every variant of the phased block, so a phase
error between the DMR and a variant gives that variant the wrong parent. This
script measures that on HG002 (`cal_v5`, ONT R10.4.1 sup v5, the third
test's calibration reads, 30-fold), against v5.0q oriented by the trio:

1. draw every read of both flowcells over the region (methylation needs the
   depth; haplotagging needs only the phased VCFs, not the phasers' reads);
2. haplotag it three ways with `whatshap haplotag`: by the truth (v5.0q
   phase at the phaser-input variants, one block, paternal first), by
   `whatshap phase`'s output, and by that output split into certified blocks
   (a block ends wherever a consecutive pair is not trusted by the per-group
   thresholds at beta 1e-3, calibrated on the even chromosomes at 30-fold);
3. `modkit pileup` (CpG, strands combined, 5mC only) per haplotype;
4. per known imprinted DMR (eLife 2022;11:e77898, supplementary file 1, hg38),
   the methylation of each haplotype, and, for the phaser and certified
   tagging, the parent assigned to every variant of the block holding the DMR.

Gate: at least 3 known DMRs on the chromosome show the listed parent's
haplotype methylated at least 0.2 above the other under truth tagging; the
parent assigned by phaser blocks is wrong for a measurable share of
variants, and the certified blocks' share is lower.

The gate passed on chr15:20-30 Mb on 15 October. The scale-up runs every
autosome with a listed DMR (383 DMRs) and draws only the reads that overlap
a DMR, since the blocks come from the stored phaser VCFs and only the
methylation needs reads. It reports two views of the assignment: every DMR
whose methylation passes the call rule under that tagging (what a tool
given the DMR list would do), and only the DMRs that the truth tagging
recovers in HG002 (so a wrong variant is a phase error, not a DMR that is
not imprinted in this genome). A variant takes the call of its nearest DMR.

Usage: `PHASEBUDGET_PLATFORM=cal_v5 PYTHONPATH=src:benchmarks .venv/bin/python
benchmarks/w2_pofo.py [--chromosomes chr15 ...]`
"""

from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import sys
import time
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import genome_common as common  # noqa: E402

OUT = common.ROOT / "data" / "w2"
DMRS = OUT / "elife-77898-supp1-v1.xlsx"  # md5 fe1c6e3e2fda4ebf0f2af84851892513
BETA, DELTA, UNIT = 1e-3, 0.05, 100_000
MIN_DIFF, MIN_CPGS, MIN_COVERAGE = 0.2, 5, 3
WHATSHAP = str(Path(sys.executable).with_name("whatshap"))
MODKIT = str(Path.home() / ".local" / "bin" / "modkit")


def known_dmrs(chrom: str) -> list:
    """`[(start, end, name, methylated parent)]` of supplementary file 1 on `chrom`."""
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    z = zipfile.ZipFile(DMRS)
    shared = ["".join(t.text or "" for t in si.iter(ns + "t"))
              for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(ns + "si")]
    out = []
    for row in ET.fromstring(z.read("xl/worksheets/sheet1.xml")).iter(ns + "row"):
        vals = []
        for c in row.findall(ns + "c"):
            v = c.find(ns + "v")
            vals.append("" if v is None else shared[int(v.text)] if c.get("t") == "s" else v.text)
        if len(vals) >= 5 and vals[0] == chrom:
            out.append((int(vals[1]), int(vals[2]), vals[3], vals[4]))
    return out


def regions(chrom: str) -> list:
    """`chrom:start-end` for the listed DMRs on `chrom`, overlapping ones merged."""
    merged: list = []
    for s, e in sorted((s, e) for s, e, _, _ in known_dmrs(chrom)):
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [f"{chrom}:{s}-{e}" for s, e in merged]


def draw(chrom: str, spans: list, d: Path) -> Path:
    """Every read of both flowcells overlapping `spans`, merged: methylation needs depth, not
    the phasers' 30-fold slice, and haplotagging needs only the phased VCFs."""
    import pysam

    out = d / "reads.bam"
    if out.exists():
        return out
    d.mkdir(parents=True, exist_ok=True)
    import os

    parts, here = [], os.getcwd()
    try:
        for k, url in enumerate(common.BAM_URLS):
            part = d / f"part{k}.bam"
            sub = d / f"remote{k}"  # htslib caches each remote index by basename, in the cwd
            sub.mkdir(exist_ok=True)
            os.chdir(sub)
            for attempt in range(6):  # S3 reads fail now and then: "retrieval of region #0 failed"
                try:
                    pysam.view("-b", "-M", "-o", str(part), url, *spans, catch_stdout=False)
                    break
                except pysam.utils.SamtoolsError:
                    if attempt == 5:
                        raise
                    time.sleep(60)
            parts.append(part)
    finally:
        os.chdir(here)
    pysam.merge("-f", str(d / "reads.part.bam"), *map(str, parts), catch_stdout=False)
    (d / "reads.part.bam").rename(out)
    pysam.index(str(out))
    for part in parts:
        part.unlink()
    return out


def truth_vcf(chrom: str, d: Path) -> tuple[Path, dict]:
    """Phaser-input variants with v5.0q's phase, paternal first, in one block; `{pos1: pat}`."""
    sc = common._module("score_hg002.py")
    variants = sc.variant_table(common.hetsnv_path(chrom))
    v5 = common.v5_truth(chrom, variants)
    trio = sc.parental_truth(common.truth_path(chrom), variants)
    both = [(a, b) for a, b in zip(v5, trio) if a is not None and b is not None]
    same = sum(a == b for a, b in both)
    flip = same < len(both) / 2  # v5.0q's allele order relative to the trio's paternal allele
    pat = {v[0]: (t if not flip else 1 - t) for v, t in zip(variants, v5) if t is not None}
    out = d / "truth.vcf"
    with gzip.open(common.hetsnv_path(chrom), "rt") as source, open(out, "w") as target:
        for line in source:
            if line.startswith("##"):
                target.write(line)
                continue
            if line.startswith("#"):
                target.write('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n')
                target.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            p = pat.get(int(f[1]))
            if p is not None:
                f[8], f[9] = "GT:PS", f"{p}|{1 - p}:1"
            else:
                f[8], f[9] = "GT", "0/1"
            target.write("\t".join(f) + "\n")
    print(f"{chrom}: v5.0q against trio transmission, {same} of {len(both)} sites agree "
          f"({'flipped' if flip else 'kept'}); {len(pat)} truth-phased variants", flush=True)
    return out, pat


def certified_vcf(chrom: str, d: Path) -> tuple[Path, Path]:
    """The phaser VCF, and a copy whose blocks end at every untrusted consecutive pair."""
    import confidence_score as cs
    import stratified_score as ss
    from stratified_gate import strata

    from phasebudget.confidence import learn_then_test_by_group

    even = [c for c in common.CONFIRMATORY if int(c[3:]) % 2 == 0 and c != chrom]
    stratum = strata([f"chr{i}" for i in range(1, 23)])  # chr20 too: W2 uses it for GNAS
    cal = ss.grouped(common.PLATFORM, "30", "whatshap", even, stratum)
    thresholds = learn_then_test_by_group(cal, BETA, DELTA)
    ind = ss.indels(common.PLATFORM, [chrom])
    breaks = set()
    for c, a, key, wrong, b, _ in cs.load_rows(common.PLATFORM, "30", "whatshap", [chrom]):
        group = "indel" if (c, a) in ind or (c, b) in ind else stratum(c, a, b)
        t = thresholds.get(group)
        if t is None or key < t:
            breaks.add(b)
    # Pairs the truth does not phase are not scored; a block also ends at them.
    phaser = common.chrom_dir(chrom) / "external" / "whatshap_30.vcf.gz"
    scored = {b for c, a, key, wrong, b, _ in cs.load_rows(common.PLATFORM, "30", "whatshap",
                                                               [chrom])}
    out = d / "certified.vcf"
    block, last_ps = None, None
    with gzip.open(phaser, "rt") as source, open(out, "w") as target:
        for line in source:
            if line.startswith("#"):
                target.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            keys, vals = f[8].split(":"), f[9].split(":")
            fields = dict(zip(keys, vals))
            pos = int(f[1])
            if "|" in fields.get("GT", "") and "PS" in fields:
                if fields["PS"] != last_ps or pos in breaks or pos not in scored:
                    block = pos
                last_ps = fields["PS"]
                fields["PS"] = str(block)
                f[9] = ":".join(fields.get(k, ".") for k in keys)
            target.write("\t".join(f) + "\n")
    print(f"{chrom}: thresholds at beta {BETA:g} from {len(even)} even chromosomes "
          f"({json.dumps({g: t for g, t in thresholds.items()})}); {len(breaks)} untrusted links",
          flush=True)
    return phaser, out


def indexed(vcf: Path, out: Path) -> Path:
    """A bgzipped, tabix-indexed copy of `vcf` (plain or gzip), which `whatshap haplotag` needs."""
    import pysam

    if not Path(str(out) + ".tbi").exists():
        plain = out.with_suffix("")
        with (gzip.open(vcf, "rb") if vcf.suffix == ".gz" else open(vcf, "rb")) as source, \
                open(plain, "wb") as target:
            target.write(source.read())
        pysam.tabix_index(str(plain), preset="vcf", force=True)
    return out


def tag_and_pileup(name: str, vcf: Path, bam: Path, chrom: str, d: Path) -> dict:
    """`{(hp, cpg position0): (modified, valid)}` after haplotagging `bam` with `vcf`."""
    ref = str(common.REFERENCE / f"{chrom}.fa")
    tagged = d / f"{name}.bam"
    vcf = indexed(vcf, d / f"{name}.tagging.vcf.gz")
    if not tagged.exists():
        subprocess.run([WHATSHAP, "haplotag", "--ignore-read-groups", "--reference", ref,
                        "--output-threads", "4", "-o", str(tagged), str(vcf), str(bam)],
                       check=True, stdout=subprocess.DEVNULL, stderr=open(d / f"{name}.log", "w"))
        subprocess.run([sys.executable, "-c", "import pysam, sys; pysam.index(sys.argv[1])",
                        str(tagged)], check=True)
    pile = d / f"pileup_{name}"
    if not pile.exists():
        subprocess.run([MODKIT, "pileup", str(tagged), str(pile), "--phased",
                        "--cpg", "--reference", ref, "--combine-strands", "--modified-bases", "5mC",
                        "--threads", "8", "--region", chrom],
                       check=True, stdout=subprocess.DEVNULL,
                       stderr=open(d / f"pileup_{name}.log", "w"))
    out = {}
    for hp in (1, 2):
        with open(pile / f"hp{hp}.bedmethyl") as handle:
            for line in handle:
                f = line.split()
                if f[3] != "m":
                    continue
                valid, modified = int(f[9]), int(f[11])
                out[(hp, int(f[1]))] = (modified, valid)
    return out


def dmr_methylation(pile: dict, start: int, end: int) -> tuple:
    """Per haplotype: (methylated share over CpGs with coverage, CpGs counted)."""
    res = []
    for hp in (1, 2):
        m = v = n = 0
        for (h, p), (mod, valid) in pile.items():
            if h == hp and start <= p < end and valid >= MIN_COVERAGE:
                m, v, n = m + mod, v + valid, n + 1
        res.append((m / v if v else None, n))
    return tuple(res)


def blocks(vcf: Path) -> list:
    """`[(pos1, ps, first allele)]` of phased calls."""
    out = []
    with open(vcf) if vcf.suffix != ".gz" else gzip.open(vcf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            fields = dict(zip(f[8].split(":"), f[9].split(":")))
            gt = fields.get("GT", "")
            if gt in ("0|1", "1|0") and "PS" in fields:
                out.append((int(f[1]), fields["PS"], int(gt[0])))
    return out


def chromosome(chrom: str) -> tuple:
    """`(text, calls, DMRs with coverage, recovered DMRs, truth-phased variants)` for `chrom`.
    A call is `(tagging, dmr, recovered, pos1, assigned paternal allele, truth's paternal
    allele, distance to the DMR's middle)`."""
    import bisect

    d = OUT / common.PLATFORM / f"{chrom}_dmrs"
    bam = draw(chrom, regions(chrom), d)
    truth, pat = truth_vcf(chrom, d)
    phaser, certified = certified_vcf(chrom, d)
    piles = {name: tag_and_pileup(name, vcf, bam, chrom, d)
             for name, vcf in (("truth", truth), ("phaser", phaser), ("certified", certified))}
    sets = {}
    for name, vcf in (("phaser", phaser), ("certified", certified)):
        calls = blocks(vcf)
        by_ps: dict = {}
        for p, ps, first in calls:
            by_ps.setdefault(ps, []).append((p, first))
        sets[name] = ([p for p, _, _ in calls], calls, by_ps)
    text = [f"\n{chrom}: DMRs; methylation of HP1 / HP2 (CpGs) per tagging; under truth tagging "
            "HP1 is paternal"]
    out, covered, recovered = [], 0, 0
    for start, end, name, parent in known_dmrs(chrom):
        mid = (start + end) // 2
        line = f"  {name[:26]:<26} {start:>10}-{end:<10} {parent:<9}"
        ok = False
        for tagging, pile in piles.items():
            (m1, n1), (m2, n2) = dmr_methylation(pile, start, end)
            line += (f" {tagging} " + (f"{m1:.2f}/{m2:.2f} ({min(n1, n2)})"
                                       if m1 is not None and m2 is not None else "-"))
            if m1 is None or m2 is None or min(n1, n2) < MIN_CPGS:
                continue
            if tagging == "truth":
                covered += 1
                ok = (m2 - m1 if parent == "Maternal" else m1 - m2) >= MIN_DIFF
                recovered += ok
                line += " recovered" if ok else ""
            elif abs(m1 - m2) >= MIN_DIFF:
                # The haplotype the listed parent methylates is that parent's.
                hp_parent = 1 if m1 > m2 else 2
                paternal_hp = hp_parent if parent == "Paternal" else 3 - hp_parent
                positions, calls, by_ps = sets[tagging]
                k = bisect.bisect_left(positions, mid)
                near = [calls[i] for i in (k - 1, k) if 0 <= i < len(calls)]
                if not near:
                    continue
                ps = min(near, key=lambda c: abs(c[0] - mid))[1]
                for p, first in by_ps[ps]:
                    if p in pat:  # HP1 carries each variant's first allele in whatshap haplotag
                        out.append((tagging, f"{chrom}:{start}", ok, p,
                                    first if paternal_hp == 1 else 1 - first, pat[p], abs(p - mid)))
        text.append(line)
    text.append(f"  {chrom}: {covered} DMRs with coverage, {recovered} recovered under truth "
                f"tagging; {len(pat)} truth-phased variants")
    print(f"  {chrom} done", flush=True)
    return "\n".join(text), out, covered, recovered, len(pat)


def main() -> None:
    from concurrent.futures import ProcessPoolExecutor

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chromosomes", nargs="+",
                        default=[f"chr{i}" for i in range(1, 23)
                                 if common.hetsnv_path(f"chr{i}").exists()])
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    print(f"W2 scale-up (development): platform {common.PLATFORM}, chromosomes "
          f"{' '.join(args.chromosomes)}, beta {BETA:g}; DMR call: listed parent's haplotype "
          f"methylated at least {MIN_DIFF} above the other over at least {MIN_CPGS} CpGs with "
          f"coverage {MIN_COVERAGE} or more per haplotype", flush=True)
    import confidence_score as cs

    # Every odd chromosome calibrates on the same even ones: fill that cache once, not in
    # parallel processes that would write it at the same time.
    cs.load_rows(common.PLATFORM, "30", "whatshap",
                 [c for c in common.CONFIRMATORY if int(c[3:]) % 2 == 0])
    with ProcessPoolExecutor(args.workers) as pool:
        results = list(pool.map(chromosome, args.chromosomes))
    calls = []
    for text, out, *_ in results:
        print(text)
        calls += out
    covered = sum(r[2] for r in results)
    recovered = sum(r[3] for r in results)
    phased = sum(r[4] for r in results)
    print(f"\nall chromosomes: {covered} DMRs with coverage, {recovered} recovered under truth "
          f"tagging; {phased} truth-phased variants")
    bins = ((0, 100_000, "under 100 kb"), (100_000, 1_000_000, "100 kb to 1 Mb"),
            (1_000_000, 10**10, "1 Mb or more"))
    print("\nparent of origin per variant, from the nearest DMR whose call passes the rule: "
          "view, tagging, DMRs used, variants assigned (share of truth-phased variants), wrong, "
          "rate; then by distance from the DMR")
    for view, keep in (("every DMR", lambda c: True), ("recovered DMRs only", lambda c: c[2])):
        best = {}
        for tagging in ("phaser", "certified"):
            best[tagging] = {}
            for c in sorted((c for c in calls if c[0] == tagging and keep(c)), key=lambda c: c[6]):
                best[tagging].setdefault(c[3], c)
            mine = list(best[tagging].values())
            wrong = sum(c[4] != c[5] for c in mine)
            dmrs = len({c[1] for c in mine})
            bad = len({c[1] for c in mine if c[4] != c[5]})
            print(f"  {view:<20} {tagging:<10} DMRs {dmrs:>4} variants {len(mine):>8} "
                  f"({len(mine) / max(1, phased):.4f}) wrong {wrong:>6} "
                  f"({wrong / max(1, len(mine)):.2e}); DMRs with a wrong variant {bad}")
            for lo, hi, label in bins:
                part = [c for c in mine if lo <= c[6] < hi]
                w = sum(c[4] != c[5] for c in part)
                print(f"      {label:<15} {len(part):>8} wrong {w:>6} ({w / max(1, len(part)):.2e})")
        # The certified blocks are short, so their variants sit close to the DMR, where the
        # phaser errs less: compare on the same variants.
        same = [best["phaser"][p] for p in best["certified"] if p in best["phaser"]]
        w = sum(c[4] != c[5] for c in same)
        print(f"  {view:<20} phaser on the certified variants: {len(same)} of "
              f"{len(best['certified'])}, wrong {w} ({w / max(1, len(same)):.2e}), DMRs with a "
              f"wrong variant {len({c[1] for c in same if c[4] != c[5]})}")


if __name__ == "__main__":
    main()
