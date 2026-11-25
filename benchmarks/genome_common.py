"""Shared paths and truth loading for the genome-wide run (`docs/design/genome-run-plan.md`)."""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path
from typing import Optional

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "genome"
GIAB = ROOT / "data" / "giab_wg"
REFERENCE = ROOT / "data" / "reference"
BAM_URL = ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
           "HG002_NA24385_son/PacBio_SequelII_CCS_11kb/HG002_GRCh38/"
           "HG002_GRCh38.haplotag.10x.bam")
CONFIRMATORY = tuple(f"chr{i}" for i in range(1, 23) if i != 20)
ODD = tuple(c for c in CONFIRMATORY if int(c[3:]) % 2 == 1)
EVEN = tuple(c for c in CONFIRMATORY if int(c[3:]) % 2 == 0)
DEPTHS = ("full", "20", "15", "10")
MODES = ("realign", "basequality")
SAMPLE = "HG002"
# `PHASEBUDGET_PLATFORM=ont` switches every script to the ONT run
# (`docs/design/ont-run-plan.md`); the variants and truths are the same.
# `hg005` switches to the HG005 replication (`docs/design/hg005-run-plan.md`).
PLATFORM = os.environ.get("PHASEBUDGET_PLATFORM", "hifi")
if PLATFORM == "hg005":
    DATA = ROOT / "data" / "genome_hg005"
    GIAB = ROOT / "data" / "giab_hg005"
    BAM_URL = ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/"
               "HG005_NA24631_son/PacBio_SequelII_CCS_11kb/HG005_GRCh38/"
               "HG005_GRCh38.haplotag.10x.bam")
    DEPTHS = ("full", "10")
    MODES = ("realign",)
    SAMPLE = "HG005"
elif PLATFORM == "ont":
    DATA = ROOT / "data" / "genome_ont"
    BAM_URL = ("https://ont-open-data.s3.amazonaws.com/giab_2023.05/analysis/variant_calling/"
               "hg002_sup_60x/hg002.haplotagged.bam")
    DEPTHS = ("full", "30", "20", "10")
    MODES = ("realign",)
elif PLATFORM in ("hg001_ont", "ont_hac"):
    # The test data of `docs/design/confidence-run-plan.md`: new reads, never used before it.
    run = {"hg001_ont": ("conf_hg001_ont", "hg001_sup_60x/hg001.haplotagged.bam", "NA12878"),
           "ont_hac": ("conf_hg002_hac", "hg002_hac_60x/hg002.haplotagged.bam", "HG002")}
    folder, path, SAMPLE = run[PLATFORM]
    DATA = ROOT / "data" / folder
    BAM_URL = ("https://ont-open-data.s3.amazonaws.com/giab_2023.05/analysis/variant_calling/"
               + path)
    DEPTHS = ("30", "20", "10")
    MODES = ("realign",)
elif PLATFORM == "hg005_ont":
    # Track B's R10 test genome: both giab_2025.01 sup flowcells, drawn from and merged.
    if not (ROOT / "docs" / "design" / "stratified-run-plan.md").exists():
        raise SystemExit("HG005 is a test genome: register the stratified run first")
    SAMPLE = "HG005"
    GIAB = ROOT / "data" / "giab_hg005"
    BAM_URLS = tuple("https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/"
                     f"HG005/{flowcell}/calls.sorted.bam" for flowcell in ("PAW87816", "PAW88001"))
    BAM_URL = BAM_URLS[0]
    DATA = ROOT / "data" / "ext_hg005_ont"
    DEPTHS = ("30", "20", "10")
    MODES = ("realign",)
elif PLATFORM in ("cal_v5", "hg001_v5"):
    # `docs/design/third-test-plan.md`: one basecaller (sup v5.0.0) for calibration and test.
    if PLATFORM == "hg001_v5" and not (ROOT / "docs" / "design" / "third-test-plan.md").exists():
        raise SystemExit("HG001 v5.0.0 is a test genome: register the third test first")
    SAMPLE, folder, cells = {"cal_v5": ("HG002", "HG002", ("PAW70337", "PAW71238")),
                             "hg001_v5": ("NA12878", "HG001", ("PAW79146", "PAW81754"))}[PLATFORM]
    BAM_URLS = tuple("https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/"
                     f"{folder}/{flowcell}/calls.sorted.bam" for flowcell in cells)
    BAM_URL = BAM_URLS[0]
    DATA = ROOT / "data" / f"ext_{PLATFORM}"
    DEPTHS = ("30", "20", "10")
    MODES = ("realign",)
elif PLATFORM == "cal_r10" or PLATFORM.startswith("ped_"):
    # Track B of the extension plan (development repository): SNVs and indels, phased by both tools.
    # Calibration genomes only, until `docs/design/stratified-run-plan.md` is committed.
    if PLATFORM == "cal_r10":
        SAMPLE = "HG002"
        BAM_URL = ("https://ont-open-data.s3.amazonaws.com/giab_2023.05/analysis/"
                   "variant_calling/hg002_sup_60x/hg002.haplotagged.bam")
    else:
        SAMPLE = PLATFORM[4:]
        if SAMPLE != "NA12878" and not (ROOT / "docs" / "design" / "stratified-run-plan.md").exists():
            raise SystemExit(f"{SAMPLE} is a test genome: register the stratified run first")
        BAM_URL = ("https://platinum-pedigree-data.s3.amazonaws.com/data/ont/mapped/GRCh38/"
                   f"{SAMPLE}.minimap2.bam")
    DATA = ROOT / "data" / f"ext_{PLATFORM}"
    DEPTHS = ("30", "20", "10")
    MODES = ("realign",)
elif PLATFORM.startswith("hprc_"):
    # Round 4, W1 (`docs/design/w1-ancestry-plan.md`): HPRC HiFi reads aligned with winnowmap.
    # HG002 calibrates against v5.0q; every other genome is a test genome scored against
    # dipcall on its own release 2 assemblies (`scripts/dipcall_truth.sh`).
    SAMPLE = PLATFORM[5:]
    HPRC_BUCKET = "https://human-pangenomics.s3.amazonaws.com/working/"
    if SAMPLE == "HG002":
        BAM_URL = (HPRC_BUCKET + "HPRC_PLUS/HG002/analysis/aligned_reads/hifi/GRCh38/"
                   "HG002_aligned_GRCh38_winnowmap.sorted.bam")
    else:
        if not (ROOT / "docs" / "design" / "w1-ancestry-plan.md").exists():
            raise SystemExit(f"{SAMPLE} is a W1 test genome: register W1 first")
        with open(ROOT / "docs" / "design" / "w1-candidates.tsv") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            rows = {r[0]: dict(zip(header, r)) for r in
                    (line.rstrip("\n").split("\t") for line in handle)}
        BAM_URL = rows[SAMPLE]["bam_url"]
        VCF_SAMPLE = "syndip"  # dipcall's sample name
    DATA = ROOT / "data" / f"ext_{PLATFORM}"
    DEPTHS = ("30", "20", "10")
    MODES = ("realign",)
elif PLATFORM.startswith("phifi_"):
    # W1b (`docs/design/w1b-pooled-plan.md`): Platinum Pedigree HiFi, scored against its truth.
    SAMPLE = PLATFORM[6:]
    if not (ROOT / "docs" / "design" / "w1b-pooled-plan.md").exists():
        raise SystemExit(f"{SAMPLE} HiFi is a W1b test genome: register W1b first")
    if SAMPLE == "NA12877" and not (ROOT / "docs" / "design" / "k1b-agreement-plan.md").exists():
        raise SystemExit("NA12877 HiFi is a K1b test genome: register K1b first")
    BAM_URL = ("https://platinum-pedigree-data.s3.amazonaws.com/data/hifi/mapped/GRCh38/"
               f"{SAMPLE}.GRCh38.haplotagged.bam")
    DATA = ROOT / "data" / f"ext_{PLATFORM}"
    DEPTHS = ("30", "20", "10")
    MODES = ("realign",)
elif PLATFORM == "k1_hg005":
    # K1b (`docs/design/k1b-agreement-plan.md`): HG005 HPRC HiFi (winnowmap), GIAB v4.2.1
    # heterozygous calls, scored against trio transmission as in the second test.
    if not (ROOT / "docs" / "design" / "k1b-agreement-plan.md").exists():
        raise SystemExit("HG005 HPRC HiFi is a K1b test genome: register K1b first")
    SAMPLE = "HG005"
    GIAB = ROOT / "data" / "giab_hg005"
    BAM_URL = ("https://human-pangenomics.s3.amazonaws.com/working/HPRC_PLUS/HG005/analysis/"
               "aligned_reads/hifi/GRCh38/HG005_aligned_GRCh38_winnowmap.sorted.bam")
    DATA = ROOT / "data" / "ext_k1_hg005"
    DEPTHS = ("30",)
    MODES = ("realign",)
elif PLATFORM == "k3_hg005":
    # K3b (`docs/design/k3b-hg005-plan.md`): the same HG005 HPRC HiFi reads, phased on the HPRC
    # release 2 dipcall variants and scored against their phase.
    if not (ROOT / "docs" / "design" / "k3b-hg005-plan.md").exists():
        raise SystemExit("k3_hg005 is a K3b test platform: register K3b first")
    SAMPLE = "HG005"
    GIAB = ROOT / "data" / "giab_hg005"
    BAM_URL = ("https://human-pangenomics.s3.amazonaws.com/working/HPRC_PLUS/HG005/analysis/"
               "aligned_reads/hifi/GRCh38/HG005_aligned_GRCh38_winnowmap.sorted.bam")
    DATA = ROOT / "data" / "ext_k3_hg005"
    DEPTHS = ("30",)
    MODES = ("realign",)
    VCF_SAMPLE = "syndip"  # dipcall's sample name
elif PLATFORM in ("w4_grch38", "w4_chm13"):
    # W4 (`docs/design/w4-reference-gate.md`): HG002 HPRC HiFi realigned with minimap2 to each
    # reference (`scripts/realign_hifi.sh`), scored against v5.0q on that reference.
    BAM_URL = str(ROOT / "data" / "w4" / f"HG002_{PLATFORM[3:]}_mm2.bam")
    if PLATFORM == "w4_chm13":
        REFERENCE = ROOT / "data" / "reference_chm13"
    DATA = ROOT / "data" / f"ext_{PLATFORM}"
    DEPTHS = ("30",)
    MODES = ("realign",)
elif PLATFORM == "q3_q100":
    # Q3 (`docs/design/q3-outside-benchmark.md`): the W4 GRCh38 minimap2 alignments, phased on
    # the Q100 dipcall heterozygous variants inside the Q100 dip BED and scored against Q100.
    BAM_URL = str(ROOT / "data" / "w4" / "HG002_grch38_mm2.bam")
    DATA = ROOT / "data" / "ext_q3_q100"
    DEPTHS = ("30",)
    MODES = ("realign",)
    VCF_SAMPLE = "syndip"  # dipcall's sample name
elif PLATFORM not in ("hifi", "hg005"):
    raise ValueError(f"unknown PHASEBUDGET_PLATFORM {PLATFORM}")
EXTENSION = (PLATFORM in ("cal_r10", "hg005_ont", "cal_v5", "hg001_v5", "w4_grch38", "w4_chm13",
                          "q3_q100", "k1_hg005", "k3_hg005")
             or PLATFORM.startswith(("ped_", "hprc_", "phifi_")))
V5_BUILD = "CHM13v2.0" if PLATFORM == "w4_chm13" else "GRCh38"  # v5.0q's file prefix
HPRC_TEST = PLATFORM.startswith("hprc_") and SAMPLE != "HG002"
if "VCF_SAMPLE" not in globals():
    VCF_SAMPLE = SAMPLE  # the sample column of the phaser input
if "BAM_URLS" not in globals():
    BAM_URLS = (BAM_URL,)
ONT_PLATFORMS = ("ont", "hg001_ont", "ont_hac")


def _module(name):
    path = Path(__file__).with_name(name)
    spec = importlib.util.spec_from_file_location(name[:-3], path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def chrom_dir(chrom: str) -> Path:
    return DATA / chrom


def fragments_path(chrom: str, depth: str, mode: str) -> Path:
    return chrom_dir(chrom) / f"fragments_{depth}_{mode}.json.gz"


def hetsnv_path(chrom: str) -> Path:
    return chrom_dir(chrom) / "hetsnv.vcf.gz"


def truth_path(chrom: str) -> Path:
    return GIAB / "truth" / f"{SAMPLE}_{chrom}_trio_truth.vcf.gz"


def hifiasm_path(chrom: str) -> Path:
    return GIAB / "split" / f"{SAMPLE}_{chrom}_hifiasm_phasetransfer.vcf.gz"


def load_truths(chrom: str, positions) -> tuple[list[Optional[int]], list[Optional[int]],
                                                  list[Optional[int]]]:
    """`(filtered, transmission, hifiasm)` paternal alleles aligned to `positions` (0-based).

    `filtered` is the registered primary truth: transmission minus every site
    where the hifiasm PATMAT call disagrees with it.
    """
    score = _module("score_hg002.py")
    variants = score.variant_table(hetsnv_path(chrom))
    if [p - 1 for p, _, _ in variants] != list(positions):
        raise ValueError(f"{chrom}: variant table and fragments disagree")
    transmission = score.parental_truth(truth_path(chrom), variants)
    hifiasm = score.parental_truth(hifiasm_path(chrom), variants)
    filtered = [t if (t is None or h is None or t == h) else None
                for t, h in zip(transmission, hifiasm)]
    return filtered, transmission, hifiasm


GIAB_V5 = ROOT / "data" / "giab_v5"


def v5_truth(chrom: str, variants) -> list[Optional[int]]:
    """First allele of GIAB v5.0q's phased call per `(pos1, ref, alt)`, None elsewhere.

    Amendment 2 of `docs/design/genome-run-plan.md`: only calls inside the
    v5.0q benchmark BED count, and only `0|1` or `1|0` at the same alleles.
    """
    import bisect
    import pysam

    if SAMPLE != "HG002":
        return [None] * len(variants)  # v5.0q exists for HG002 only
    starts, ends = [], []
    with open(GIAB_V5 / f"HG002_{V5_BUILD}_v5.0q_smvar.benchmark.bed") as handle:
        for line in handle:
            f = line.split("\t")
            if f[0] == chrom:
                starts.append(int(f[1]))
                ends.append(int(f[2]))
    calls = {}
    with pysam.TabixFile(str(GIAB_V5 / f"HG002_{V5_BUILD}_v5.0q_smvar.vcf.gz")) as tabix:
        for line in tabix.fetch(chrom):
            f = line.split("\t")
            gt = f[9].split(":")[0]
            if gt in ("0|1", "1|0"):
                calls[(int(f[1]), f[3], f[4])] = int(gt[0])
    truth: list[Optional[int]] = []
    for pos1, ref, alt in variants:
        k = bisect.bisect_right(starts, pos1 - 1) - 1
        inside = k >= 0 and pos1 - 1 < ends[k]
        truth.append(calls.get((pos1, ref, alt)) if inside else None)
    return truth


PLATINUM = ROOT / "data" / "platinum"


def platinum_truth(chrom: str, variants) -> list[Optional[int]]:
    """First allele of Platinum Genomes NA12878's phased call per `(pos1, ref, alt)`.

    Only calls inside `ConfidentRegions.bed.gz` count, and only `0|1` or `1|0`
    at the same alleles, by the rules `v5_truth` uses.
    """
    import bisect
    import pysam

    starts, ends = [], []
    with pysam.TabixFile(str(PLATINUM / "ConfidentRegions.bed.gz")) as tabix:
        for line in tabix.fetch(chrom):
            f = line.split("\t")
            starts.append(int(f[1]))
            ends.append(int(f[2]))
    calls = {}
    with pysam.TabixFile(str(PLATINUM / "NA12878.vcf.gz")) as tabix:
        for line in tabix.fetch(chrom):
            f = line.split("\t")
            gt = f[9].split(":")[0]
            if gt in ("0|1", "1|0"):
                calls[(int(f[1]), f[3], f[4])] = int(gt[0])
    truth: list[Optional[int]] = []
    for pos1, ref, alt in variants:
        k = bisect.bisect_right(starts, pos1 - 1) - 1
        inside = k >= 0 and pos1 - 1 < ends[k]
        truth.append(calls.get((pos1, ref, alt)) if inside else None)
    return truth


PEDIGREE = ROOT / "data" / "pedigree_truth"


def pedigree_truth(chrom: str, variants, sample: Optional[str] = None) -> list[Optional[int]]:
    """First allele of the Platinum Pedigree truthset's phased call for `sample` (default
    `SAMPLE`), per `(pos1, ref, alt)`.

    Only `0|1` or `1|0` calls at the same alleles inside `hq_regions_final.bed.gz` count.
    """
    import bisect
    import gzip

    import pysam

    starts, ends = [], []
    with gzip.open(PEDIGREE / "hq_regions_final.bed.gz", "rt") as handle:
        for line in handle:
            f = line.split("\t")
            if f[0] == chrom:
                starts.append(int(f[1]))
                ends.append(int(f[2]))
    calls = {}
    with pysam.TabixFile(str(PEDIGREE / "CEPH1463.GRCh38.family-truthset.ov.vcf.gz")) as tabix:
        column = 9 + tabix.header[-1].split("\t")[9:].index(sample or SAMPLE)
        for line in tabix.fetch(chrom):
            f = line.split("\t")
            gt = f[column].split(":")[0]
            if gt in ("0|1", "1|0"):
                calls[(int(f[1]), f[3], f[4])] = int(gt[0])
    truth: list[Optional[int]] = []
    for pos1, ref, alt in variants:
        k = bisect.bisect_right(starts, pos1 - 1) - 1
        inside = k >= 0 and pos1 - 1 < ends[k]
        truth.append(calls.get((pos1, ref, alt)) if inside else None)
    return truth


HPRC = ROOT / "data" / "hprc"
Q100_VCF = GIAB_V5 / "GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.vcf.gz"  # not indexed
Q100_BED = GIAB_V5 / "GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.bed"


def dipcall_truth(chrom: str, variants, vcf: Optional[Path] = None,
                  bed: Optional[Path] = None) -> list[Optional[int]]:
    """First allele of dipcall's call on `SAMPLE`'s assemblies per `(pos1, ref, alt)`.

    The first allele is the paternal assembly's. Only `0|1` or `1|0` calls at
    the same alleles inside dipcall's BED count, by the rules `v5_truth` uses.
    `vcf` and `bed` default to `SAMPLE`'s HPRC files; a VCF without a tabix
    index is streamed.
    """
    import bisect
    import gzip

    import pysam

    vcf = vcf or HPRC / SAMPLE / f"{SAMPLE}.dip.vcf.gz"
    starts, ends = [], []
    with open(bed or HPRC / SAMPLE / f"{SAMPLE}.dip.bed") as handle:
        for line in handle:
            f = line.split("\t")
            if f[0] == chrom:
                starts.append(int(f[1]))
                ends.append(int(f[2]))
    calls = {}
    if Path(str(vcf) + ".tbi").exists():
        with pysam.TabixFile(str(vcf)) as tabix:
            lines = list(tabix.fetch(chrom))
    else:
        with gzip.open(vcf, "rt") as handle:
            lines = [line for line in handle if line.startswith(chrom + "\t")]
    for line in lines:
        f = line.split("\t")
        gt = f[9].split(":")[0]
        if gt in ("0|1", "1|0"):
            calls[(int(f[1]), f[3], f[4])] = int(gt[0])
    truth: list[Optional[int]] = []
    for pos1, ref, alt in variants:
        k = bisect.bisect_right(starts, pos1 - 1) - 1
        inside = k >= 0 and pos1 - 1 < ends[k]
        truth.append(calls.get((pos1, ref, alt)) if inside else None)
    return truth


KG_PANEL = ROOT / "data" / "kg_panel"


def kg_parents(sample: str) -> Optional[tuple[str, str]]:
    """`(father, mother)` of `sample` in the 1000 Genomes 3,202-sample pedigree, if both listed."""
    with open(KG_PANEL / "20130606_g1k_3202_samples_ped_population.txt") as handle:
        for line in handle:
            f = line.split()
            if f[1] == sample and f[2] != "0" and f[3] != "0":
                return f[2], f[3]
    return None


def transmitted(child: Optional[int], father: Optional[int], mother: Optional[int]
                ) -> Optional[int]:
    """The paternal allele of a heterozygous child, from the trio's alt-allele counts.

    None when a count is missing, the child is not heterozygous, both parents
    are heterozygous, or no assignment is Mendelian.
    """
    if child != 1 or father is None or mother is None:
        return None
    options = [p for p in (0, 1)
               if p in (father // 2, (father + 1) // 2) and 1 - p in (mother // 2, (mother + 1) // 2)]
    return options[0] if len(options) == 1 else None


def _alt_count(field: str) -> Optional[int]:
    alleles = field.split(":")[0].replace("|", "/").split("/")
    return sum(int(a) for a in alleles) if len(alleles) == 2 and set(alleles) <= {"0", "1"} else None


def kg_transmission(chrom: str, variants, sample: Optional[str] = None) -> list[Optional[int]]:
    """Paternal allele per `(pos1, ref, alt)` by transmission in the 1000 Genomes panel.

    Reads the trio subset of the high-coverage phased panel (`scripts/kg_trios.py`), using
    its genotypes, not its phase. A variant with no record at the same alleles gets None.
    A trio's labels on a chromosome are computed once and kept in `data/kg_panel/transmission/`,
    since every depth and phaser of a genome asks for them.
    """
    import json

    import pysam

    sample = sample or SAMPLE
    cache = KG_PANEL / "transmission" / f"{sample}_{chrom}.json"
    if cache.exists():
        labels = {(p, r, a): x for p, r, a, x in json.loads(cache.read_text())}
    else:
        father, mother = kg_parents(sample)
        labels = {}
        with pysam.TabixFile(str(KG_PANEL / "trios" / f"{chrom}.vcf.gz")) as tabix:
            names = tabix.header[-1].split("\t")[9:]
            columns = [9 + names.index(s) for s in (sample, father, mother)]
            for line in tabix.fetch(chrom):
                f = line.split("\t")
                x = transmitted(*(_alt_count(f[c]) for c in columns))
                if x is not None:
                    labels[(int(f[1]), f[3], f[4])] = x
        cache.parent.mkdir(exist_ok=True)
        tmp = cache.with_name(f"{cache.name}.{os.getpid()}.tmp")
        tmp.write_text(json.dumps([[*k, x] for k, x in labels.items()]))
        tmp.replace(cache)
    return [labels.get(tuple(v)) for v in variants]


def consensus_truth(chrom: str, variants) -> list[Optional[int]]:
    """Dipcall's first allele where 1000 Genomes transmission gives the same paternal allele
    (`docs/design/w1-ancestry-plan.md`, amendment 1); None elsewhere."""
    return [d if d is not None and d == k else None
            for d, k in zip(dipcall_truth(chrom, variants), kg_transmission(chrom, variants))]


def truth_at(chrom: str, variants, positions) -> list[Optional[int]]:
    """`primary_truth` aligned to the extractor's 0-based `positions`, matched by position.

    The extractor can drop input variants (it does with indels), so the
    truth list of `variants` must not be indexed by fragment column.
    """
    by_position = {v[0]: t for v, t in zip(variants, primary_truth(chrom, variants))}
    return [by_position.get(p + 1) for p in positions]


def primary_truth(chrom: str, variants) -> list[Optional[int]]:
    """The truth each dataset is scored against: v5.0q for HG002, Platinum for HG001.

    In track B's pedigree arm it is the Platinum Pedigree truthset for every member.
    """
    if PLATFORM.startswith(("ped_", "phifi_")) or PLATFORM == "hg001_v5":
        return pedigree_truth(chrom, variants)
    if PLATFORM == "q3_q100":
        return dipcall_truth(chrom, variants, Q100_VCF, Q100_BED)
    if PLATFORM == "k3_hg005":
        return dipcall_truth(chrom, variants)
    if PLATFORM.startswith("hprc_") and os.environ.get("PHASEBUDGET_TRUTH") == "dipcall":
        return dipcall_truth(chrom, variants)  # W1's A6, descriptive
    if HPRC_TEST:
        return consensus_truth(chrom, variants)
    if SAMPLE == "NA12878":
        return platinum_truth(chrom, variants)
    if SAMPLE == "HG002":
        return v5_truth(chrom, variants)
    if SAMPLE == "HG005":
        return load_truths(chrom, [p - 1 for p, _, _ in variants])[1]  # transmission
    raise ValueError(f"no primary truth for {SAMPLE} here")
