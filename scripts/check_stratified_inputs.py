"""Check that every file `benchmarks/stratified_score.py` reads exists and decodes.

Run before the registered scoring. No truth is compared with any phasing and no
count of phased sites is printed, so the check reveals nothing about S1 to S8.
With `--third` it checks the inputs of `benchmarks/third_score.py` instead
(`docs/design/third-test-plan.md`), and with `--w1` those of `w1_score.py` and
`w1b_score.py` (`docs/design/w1-ancestry-plan.md`, `docs/design/w1b-pooled-plan.md`).

Usage: `.venv/bin/python scripts/check_stratified_inputs.py [--third | --w1]`
"""
import gzip
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "benchmarks"))
import stratify  # noqa: E402

CHROMS = [f"chr{i}" for i in range(1, 23) if i != 20]
PLATS = {"cal_r10": "HG002", "ped_NA12878": "NA12878", "hg005_ont": "HG005",
         **{f"ped_{s}": s for s in ("NA12877", "NA12879", "NA12881", "NA12882", "NA12885", "NA12886")}}
TESTS = [p for p in PLATS if p not in ("cal_r10", "ped_NA12878")]
if "--third" in sys.argv:
    PLATS, TESTS = {"cal_v5": "HG002", "cal_r10": "HG002", "hg001_v5": "NA12878"}, ["hg001_v5"]
W1 = "--w1" in sys.argv
if W1:
    import genome_common
    from w1b_score import PEDIGREE
    from w1_score import candidates

    W1_SAMPLES = [s for s in candidates() if genome_common.kg_parents(s)]
    PLATS = {"hprc_HG002": "HG002", **{f"hprc_{s}": "syndip" for s in W1_SAMPLES},
             **{f"phifi_{s}": s for s in PEDIGREE}}
    TESTS = [p for p in PLATS if p != "hprc_HG002"]
problems, checked, pending = [], 0, []


def gz_lines(path):
    with gzip.open(path, "rt") as h:  # reading to the end checks every block's CRC
        return h.readlines()


def vcf_ok(path, sample):
    lines = gz_lines(path)
    head = [l for l in lines if l.startswith("#CHROM")]
    if not head or head[0].rstrip("\n").split("\t")[9:] != [sample]:
        return f"sample column {head[0].split()[9:] if head else None}, expected {sample}"
    if not any(not l.startswith("#") for l in lines):
        return "no records"
    if not lines[-1].endswith("\n"):
        return "truncated last line"
    return None


for plat, sample in PLATS.items():
    for chrom in CHROMS:
        d = ROOT / "data" / f"ext_{plat}" / chrom
        if not (ROOT / "data" / f"ext_{plat}" / "done" / f"{chrom}__cleanup.ok").exists():
            pending.append(f"{plat}:{chrom}")
            continue
        try:
            for f in ("hetsnv.vcf.gz", "hetsnv.vcf.gz.tbi"):
                if not (d / f).exists():
                    problems.append(f"{plat} {chrom}: {f} missing")
            err = vcf_ok(d / "hetsnv.vcf.gz", sample)
            if err:
                problems.append(f"{plat} {chrom} hetsnv: {err}")
            for depth in ("30", "20", "10"):
                frag = d / f"fragments_{depth}_realign.json.gz"
                json.loads("".join(gz_lines(frag)))
                for ph in ("whatshap", "longphase"):
                    v = d / "external" / f"{ph}_{depth}.vcf"
                    v = v if v.exists() else v.with_name(v.name + ".gz")
                    if not v.exists():
                        problems.append(f"{plat} {chrom}: {ph}_{depth}.vcf(.gz) missing")
                        continue
                    if v.suffix == ".gz":
                        err = vcf_ok(v, sample)
                    else:
                        lines = v.read_text().splitlines()
                        err = None if any(l and not l.startswith("#") for l in lines) else "no records"
                    if err:
                        problems.append(f"{plat} {chrom} {v.name}: {err}")
            checked += 1
        except Exception as e:  # noqa: BLE001
            problems.append(f"{plat} {chrom}: {type(e).__name__}: {e}")

shared = [ROOT / "data/pedigree_truth" / n for n in
          ("hq_regions_final.bed.gz", "CEPH1463.GRCh38.family-truthset.ov.vcf.gz",
           "CEPH1463.GRCh38.family-truthset.ov.vcf.gz.tbi")]
shared += [ROOT / "data/giab_hg005" / sub / f"HG005_{c}_{kind}.vcf.gz" for c in CHROMS
           for sub, kind in (("truth", "trio_truth"), ("split", "hifiasm_phasetransfer"))]
shared += [ROOT / "data/gencode/gencode.v50.basic.annotation.gtf.gz",
           ROOT / "data/clinvar/gene_specific_summary.txt"]
if W1:
    shared += [ROOT / "data/hprc" / s / f"{s}.dip.{x}" for s in ["HG002", *W1_SAMPLES]
               for x in ("vcf.gz", "vcf.gz.tbi", "bed")]
    shared += [ROOT / "data/kg_panel/trios" / f"{c}.vcf.gz{x}" for c in CHROMS for x in ("", ".tbi")]
    shared += [ROOT / "data/kg_panel/20130606_g1k_3202_samples_ped_population.txt"]
    check = ROOT / "data/kg_panel/md5.check"
    if not check.exists() or any(not line.rstrip().endswith(": OK") for line in check.open()):
        problems.append("panel md5 check missing or not all OK")
for p in shared:
    if not p.exists():
        problems.append(f"missing {p.relative_to(ROOT)}")
for name, rel, md5 in stratify.STRATA:
    p = stratify.DIR / Path(rel).name
    if not p.exists() or stratify._md5(p) != md5:
        problems.append(f"stratum {name}: missing or md5 mismatch")

cache = ROOT / "data" / "confidence_rows"
stale = [f.name for f in cache.iterdir() if any(t in f.name for t in TESTS)]
if stale:
    problems.append(f"test-genome caches already present: {stale}")

print(f"checked {checked} finished chromosomes; pending {len(pending)}: {' '.join(pending)}")
print(f"shared files {len(shared)}, strata {len(stratify.STRATA)}")
print("problems:" if problems else "no problems")
for p in problems:
    print("  " + p)
sys.exit(1 if problems else 0)
