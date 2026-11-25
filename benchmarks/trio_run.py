"""The registered HG002 trio run on real reads (`docs/design/trio-run-plan.md`).

Per chromosome:
1. the child's heterozygous SNVs of the genome-wide HiFi run, kept where
   both parents lie inside their v4.2.1 BEDs with a compatible call;
2. each parent's reads streamed for that chromosome, extracted at the kept
   sites with WhatsHap's reader, and deleted;
3. the child's full-depth fragments projected onto the kept sites;
4. every allocation selected, phased with WhatsHap's `PedigreeDPTable`, and
   the child's alleles and blocks written with the DP's time.

`--score` then scores the written runs against GIAB v5.0q and decides P1 to
P3. Usage:
`PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/trio_run.py --chromosomes chr20`
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import genome_common as common  # noqa: E402

from phasebudget._types import Fragment, Phasing  # noqa: E402
from phasebudget.bootstrap import statistic_interval  # noqa: E402
from phasebudget.fragments import FragmentMatrix  # noqa: E402
from phasebudget.metrics import block_lengths, comparable_pairs, n50, switch_positions  # noqa: E402
from phasebudget.paired import common_pairs, sign_test, window_discordance  # noqa: E402
from phasebudget.selection import read_weights, select_optimal  # noqa: E402

extractor = common._module("extract_fragments.py")
score = common._module("score_hg002.py")
pedigree = common._module("benchmark_pedigree.py")
pipelines = common._module("hg002_pipelines.py")

TRIO = common.ROOT / "data" / "trio"
GIAB_DATA = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio"
PARENTS = {
    "father": ("HG003", f"{GIAB_DATA}/HG003_NA24149_father/PacBio_CCS_15kb_20kb_chemistry2/GRCh38/"
                        "HG003.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam"),
    "mother": ("HG004", f"{GIAB_DATA}/HG004_NA24143_mother/PacBio_CCS_15kb_20kb_chemistry2/GRCh38/"
                        "HG004.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam"),
}
CONFIRMATORY = ("chr16", "chr18", "chr21", "chr22")
RUNS = (("thirds", 15), ("thirds", 9), ("shared", 9), ("shared", 15),
        # Amendment 2: WhatsHap's own selector per member, as its pedigree mode deploys it.
        ("whatshap_thirds", 15), ("whatshap_thirds", 9))
MEMBERS = ("child", "father", "mother")
WINDOW = 1_000_000


def child_paths(chrom):
    """`(variant VCF, full-depth fragments)` of HG002; chromosome 20 is the pilot's."""
    if chrom == "chr20":
        pilot = common.ROOT / "data" / "hg002"
        return pilot / "HG002_chr20_hetsnv_bed.vcf", pilot / "HG002_chr20_fragments.json.gz"
    return common.hetsnv_path(chrom), common.fragments_path(chrom, "full", "realign")


def parent_genotypes(chrom, variants, sample):
    """Genotype `(a, b)` per child variant, or None where the site is dropped."""
    split = common.GIAB / "split"
    starts, ends = [], []
    with open(split / f"{sample}_{chrom}_v4.2.1.bed") as handle:
        for line in handle:
            f = line.split("\t")
            starts.append(int(f[1]))
            ends.append(int(f[2]))
    calls = {}
    with gzip.open(split / f"{sample}_{chrom}_v4.2.1.vcf.gz", "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            calls[int(f[1])] = (f[3], f[4], f[9].split(":")[0].replace("|", "/"))
    import bisect

    out = []
    for pos1, ref, alt in variants:
        k = bisect.bisect_right(starts, pos1 - 1) - 1
        if k < 0 or pos1 - 1 >= ends[k]:
            out.append(None)
            continue
        if pos1 not in calls:
            out.append((0, 0))
            continue
        cref, calt, gt = calls[pos1]
        if (cref, calt) != (ref, alt) or gt not in ("0/1", "1/0", "1/1", "0/0"):
            out.append(None)
            continue
        out.append(tuple(sorted(int(x) for x in gt.split("/"))))
    return out


def project(matrix, keep):
    """The matrix restricted to the columns in `keep`, renumbered in order."""
    new = {old: i for i, old in enumerate(keep)}
    fragments = []
    for f in matrix:
        calls = [(new[j], a, w) for j, a, w in f.observed() if j in new]
        if not calls:
            continue
        start, last = calls[0][0], calls[-1][0]
        alleles = [None] * (last - start + 1)
        weights = [0.0] * (last - start + 1)
        for j, a, w in calls:
            alleles[j - start] = a
            weights[j - start] = w
        fragments.append(Fragment(len(fragments), start, tuple(alleles), tuple(weights)))
    return FragmentMatrix(fragments, len(keep))


def prepare(chrom):
    """Kept sites, genotypes and the three members' matrices; parents cached on disk."""
    d = TRIO / chrom
    d.mkdir(parents=True, exist_ok=True)
    variant_vcf, child_fragments = child_paths(chrom)
    variants = score.variant_table(variant_vcf)
    genotypes = {"father": parent_genotypes(chrom, variants, "HG003"),
                 "mother": parent_genotypes(chrom, variants, "HG004")}
    # Amendment 1: a heterozygous child cannot have two parents homozygous for one allele.
    keep = [j for j in range(len(variants))
            if genotypes["father"][j] is not None and genotypes["mother"][j] is not None
            and not (genotypes["father"][j] == genotypes["mother"][j]
                     and genotypes["father"][j][0] == genotypes["father"][j][1])]
    kept_variants = [variants[j] for j in keep]
    vcf = d / "kept.vcf.gz"
    if not vcf.exists():
        import pysam

        wanted = set(kept_variants)
        source = pysam.VariantFile(str(variant_vcf))
        with pysam.VariantFile(str(vcf), "wz", header=source.header) as target:
            for record in source:
                if (record.pos, record.ref, record.alts[0]) in wanted:
                    target.write(record)
        pysam.tabix_index(str(vcf), preset="vcf", force=True)
    child, positions, _, _ = extractor.load_fragments(str(child_fragments))
    matrices = {"child": extractor.phred_to_log_odds(project(child, keep))}
    kept_positions = [positions[j] for j in keep]
    for member, (sample, url) in PARENTS.items():
        path = d / f"{sample}_fragments.json.gz"
        if not path.exists():
            extract_parent(chrom, sample, url, vcf, path)
        matrix, parent_positions, _, _ = extractor.load_fragments(str(path))
        if parent_positions != kept_positions:
            raise ValueError(f"{chrom} {sample}: extracted sites differ from the kept sites")
        matrices[member] = extractor.phred_to_log_odds(matrix)
    genomes = {"child": [(0, 1)] * len(keep)}
    for member in ("father", "mother"):
        genomes[member] = [genotypes[member][j] for j in keep]
    return kept_variants, kept_positions, genomes, matrices, len(variants) - len(keep)


def extract_parent(chrom, sample, url, vcf, out):
    import pysam
    import shutil

    d = out.parent
    import time

    # The plan requires 20 GB free before a download; wait for it rather than fail, since
    # the ONT run frees its slices as it goes.
    while shutil.disk_usage(d).free < 20e9:
        time.sleep(300)
    here = os.getcwd()
    os.chdir(d)  # htslib caches the remote index in the working directory
    try:
        bam = d / f"{sample}.bam"
        pysam.view("-b", "-o", str(bam), url, chrom, catch_stdout=False)
        pysam.index(str(bam))
        matrix, positions, names = extractor.read_fragments(
            str(bam), str(vcf), str(common.REFERENCE / f"{chrom}.fa"), chrom, "HG002", 20,
            ignore_read_groups=True)
        extractor.save_fragments(str(out), matrix, positions, names,
                                 {"sample": sample, "bam": url, "chromosome": chrom,
                                  "reference": f"{chrom}.fa", "mapq": 20,
                                  "ignore_read_groups": True})
    finally:
        os.chdir(here)
        for path in d.glob(f"{sample}.bam*"):
            path.unlink()
        for path in d.glob("*.bai"):
            path.unlink()


def select(allocation, total, matrices, positions):
    """As in `benchmark_trio_budget.select`, for the thirds and shared arms."""
    if allocation == "whatshap_thirds":
        per = total // 3
        return {m: matrices[m].subset(
                    pipelines.whatshap_selection(matrices[m], per, positions).selected)
                for m in MEMBERS}
    if allocation == "thirds":
        per = total // 3
        return {m: matrices[m].subset(select_optimal(matrices[m], per,
                                                     read_weights(matrices[m])).selected)
                for m in MEMBERS}
    fragments, origin = [], []
    for m in MEMBERS:
        for f in matrices[m]:
            fragments.append(Fragment(len(fragments), f.start, f.alleles, f.weights))
            origin.append((m, f.read_id))
    union = FragmentMatrix(fragments, matrices["child"].n_variants)
    chosen = select_optimal(union, total, read_weights(union)).selected
    picked = {m: [] for m in MEMBERS}
    for index in chosen:
        m, read_id = origin[index]
        picked[m].append(read_id)
    return {m: matrices[m].subset(picked[m]) for m in MEMBERS}


def phase(selected, genomes, positions):
    """The child's alleles and blocks from WhatsHap's pedigree DP, and the DP time.

    As `whatshap phase` does (amendment 3), the DP, the genotypes and the
    master block cover only the positions some selected read calls; every
    other kept site is left unphased.
    """
    from whatshap.cli.phase import find_components
    from whatshap.core import (Genotype, NumericSampleIds, Pedigree, PedigreeDPTable,
                               ReadSet)
    from whatshap.pedigree import UniformRecombinationCostComputer

    ids = NumericSampleIds()
    readset = ReadSet()
    for member in MEMBERS:
        pedigree._add_reads(readset, selected[member], positions, ids[member], member[0])
    readset.sort()
    accessible = sorted(readset.get_positions())
    index = {p: j for j, p in enumerate(positions)}
    ped = Pedigree(ids)
    for member in MEMBERS:
        ped.add_individual(member, [Genotype(list(genomes[member][index[p]])) for p in accessible])
    ped.add_relationship(father_id="father", mother_id="mother", child_id="child")
    costs = UniformRecombinationCostComputer.uniform_recombination_map(
        pedigree.RECOMBINATION_RATE, accessible)
    started = time.perf_counter()
    table = PedigreeDPTable(readset, costs, ped, False, accessible)
    superreads, _ = table.get_super_reads()
    seconds = time.perf_counter() - started
    child = superreads[MEMBERS.index("child")]
    alleles = {variant.position: variant.allele for variant in child[0]}
    # A set, as in WhatsHap: a site homozygous in both parents appears once.
    homozygous = sorted({p for member in ("father", "mother") for p in accessible
                         if genomes[member][index[p]][0] == genomes[member][index[p]][1]})
    components = find_components(accessible, readset, master_block=homozygous)
    return ([alleles.get(p, -1) for p in positions],
            [components.get(p, -1) for p in positions], seconds)


def run_chromosome(chrom):
    kept_variants, positions, genomes, matrices, dropped = prepare(chrom)
    out = TRIO / chrom / "runs.json.gz"
    record = {"chrom": chrom, "variants": [list(v) for v in kept_variants], "dropped": dropped,
              "reads": {m: len(matrices[m]) for m in MEMBERS}, "runs": {}}
    for allocation, total in RUNS:
        selected = select(allocation, total, matrices, positions)
        coverage = [0] * len(positions)
        for m in MEMBERS:
            for j, c in enumerate(selected[m].coverage()):
                coverage[j] += c
        alleles, blocks, seconds = phase(selected, genomes, positions)
        record["runs"][f"{allocation}_{total}"] = {
            "alleles": alleles, "blocks": blocks, "seconds": round(seconds, 2),
            "max_active": max(coverage, default=0),
            "selected": {m: len(selected[m]) for m in MEMBERS}}
        print(f"{chrom} {allocation} {total}: DP {seconds:.1f} s, max active "
              f"{max(coverage, default=0)}", flush=True)
    with gzip.open(out, "wt") as handle:
        json.dump(record, handle)


def load_run(chrom, name):
    with gzip.open(TRIO / chrom / "runs.json.gz", "rt") as handle:
        record = json.load(handle)
    variants = [tuple(v) for v in record["variants"]]
    run = record["runs"][name]
    phasing = Phasing(tuple(None if a < 0 else a for a in run["alleles"]),
                      tuple(None if b < 0 else b for b in run["blocks"]))
    return variants, phasing, run, record


def decide(chromosomes):
    names = [f"{a}_{t}" for a, t in RUNS]
    pooled = {n: {"switches": 0, "pairs": 0, "blocks": [], "seconds": 0.0, "max_active": 0}
              for n in names}
    tests = {("shared_9", "thirds_15"): [0, 0, 0, []], ("shared_9", "thirds_9"): [0, 0, 0, []],
             ("thirds_15", "whatshap_thirds_15"): [0, 0, 0, []]}
    dropped = 0
    for chrom in chromosomes:
        variants, _, _, record = load_run(chrom, names[0])
        dropped += record["dropped"]
        truth = common.v5_truth(chrom, variants)
        positions = [p - 1 for p, _, _ in variants]
        phasings = {}
        for n in names:
            _, phasing, run, _ = load_run(chrom, n)
            phasings[n] = phasing
            p = pooled[n]
            p["switches"] += len(switch_positions(truth, phasing))
            p["pairs"] += comparable_pairs(phasing, truth)
            p["blocks"].extend(block_lengths(phasing, positions))
            p["seconds"] += run["seconds"]
            p["max_active"] = max(p["max_active"], run["max_active"])
        for (a, b), q in tests.items():
            rows = common_pairs(truth, phasings[a], phasings[b])
            q[0] += len(rows)
            q[1] += sum(1 for _, x, y in rows if x and not y)
            q[2] += sum(1 for _, x, y in rows if y and not x)
            q[3].extend(window_discordance(rows, positions, WINDOW))
    print(f"trio run, chromosomes {', '.join(chromosomes)}; child sites dropped for parents "
          f"outside their BED or with incompatible calls: {dropped}; truth GIAB v5.0q")
    print(f"{'run':>10} {'switches':>9} {'pairs':>9} {'per 10k':>8} {'N50 kb':>8} "
          f"{'DP s':>9} {'max active':>10}")
    for n, p in pooled.items():
        print(f"{n:>10} {p['switches']:>9} {p['pairs']:>9} "
              f"{1e4 * p['switches'] / max(p['pairs'], 1):>8.3f} {n50(p['blocks']) / 1e3:>8.0f} "
              f"{p['seconds']:>9.1f} {p['max_active']:>10}")

    def difference(windows):
        def statistic(picks):
            pairs = sum(windows[k][0] for k in picks)
            return None if pairs == 0 else \
                1e4 * (sum(windows[k][1] for k in picks) - sum(windows[k][2] for k in picks)) / pairs
        return statistic_interval(statistic, len(windows), resamples=2000, seed=0)

    q = tests[("shared_9", "thirds_15")]
    interval = difference(q[3])
    ratio = pooled["shared_9"]["seconds"] / pooled["thirds_15"]["seconds"]
    p1 = interval.high < 0.5 and ratio <= 0.1
    print(f"P1 shared 9 vs thirds 15: discordant {q[1]}/{q[2]} of {q[0]} common pairs; "
          f"difference per 10k {interval.estimate:+.3f} [{interval.low:+.3f}, "
          f"{interval.high:+.3f}]; DP time ratio {ratio:.3f} -> "
          f"{'supported' if p1 else 'not supported'} (rule: upper end below +0.5 and ratio at "
          "most 0.1)")
    q = tests[("shared_9", "thirds_9")]
    pv = sign_test(q[1], q[1] + q[2])
    p2 = pv < 0.05 and q[1] < q[2]
    print(f"P2 shared 9 vs thirds 9: discordant {q[1]}/{q[2]}, p {pv:.3g} -> "
          f"{'supported' if p2 else 'not supported'} (rule: significant for shared at 0.05)")
    q = tests[("thirds_15", "whatshap_thirds_15")]
    pv = sign_test(q[1], q[1] + q[2])
    p4 = pv < 0.05 and q[1] < q[2]
    print(f"P4 exact thirds 15 vs WhatsHap thirds 15: discordant {q[1]}/{q[2]} of {q[0]} common "
          f"pairs, p {pv:.3g} -> {'supported' if p4 else 'not supported'} (rule: significant for "
          "exact selection at 0.05)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chromosomes", nargs="+", default=list(CONFIRMATORY))
    parser.add_argument("--score", action="store_true")
    args = parser.parse_args()
    if args.score:
        decide(args.chromosomes)
        return
    for chrom in args.chromosomes:
        if (TRIO / chrom / "runs.json.gz").exists():
            print(f"{chrom}: done")
            continue
        run_chromosome(chrom)


if __name__ == "__main__":
    main()
