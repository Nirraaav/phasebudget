"""List HPRC genomes that have release 2 parental assemblies and a HiFi GRCh38 BAM.

Round 4, direction W1 (`docs/design/w1-ancestry-plan.md`). For every sample
folder under `working/HPRC/` and `working/HPRC_PLUS/` of the public
`human-pangenomics` bucket, a sample is a candidate when the bucket holds
`<s>_{mat,pat}_hprc_r2_*.fa.gz` under `assemblies/release2/` and
`<s>_aligned_GRCh38_winnowmap.sorted.bam` with its `.bai` under
`analysis/aligned_reads/hifi/GRCh38/`. Population and sex come from the 1000
Genomes 3,202-sample table, then from the phase 3 pedigree file, whose
population codes map to super-populations through the first table; a sample
in neither gets `NA`.

Usage: `python scripts/hprc_candidates.py > docs/design/w1-candidates.tsv`
"""

from __future__ import annotations

import hashlib
import re
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

BUCKET = "https://human-pangenomics.s3.amazonaws.com/"
POPULATIONS = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
               "1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt")
PEDIGREE = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
            "integrated_call_samples_v3.20200731.ALL.ped")
SEX = {"1": "male", "2": "female"}


def fetch(url: str, attempts: int = 5) -> str:
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode()
        except OSError:
            if attempt == attempts - 1:
                raise
            time.sleep(5)


def listing(prefix: str, delimiter: bool = False) -> tuple[list[str], dict[str, int]]:
    """Common prefixes and `{key: size}` under `prefix`, following continuation tokens."""
    prefixes, keys, token = [], {}, None
    while True:
        query = {"list-type": "2", "prefix": prefix}
        if delimiter:
            query["delimiter"] = "/"
        if token:
            query["continuation-token"] = token
        text = fetch(BUCKET + "?" + urllib.parse.urlencode(query))
        prefixes += re.findall(r"<Prefix>([^<]+)</Prefix>", text)[1:]  # the first is `prefix`
        for key, size in re.findall(r"<Key>([^<]+)</Key>.*?<Size>(\d+)</Size>", text):
            keys[key] = int(size)
        match = re.search(r"<NextContinuationToken>([^<]+)<", text)
        if not match:
            return prefixes, keys
        token = match.group(1)


def files(folder: str) -> tuple:
    sample = folder.rstrip("/").rsplit("/", 1)[1]
    _, assemblies = listing(folder + "assemblies/release2/", delimiter=True)
    _, reads = listing(folder + "analysis/aligned_reads/hifi/GRCh38/", delimiter=True)
    pick = {}
    for hap in ("mat", "pat"):
        found = sorted(k for k in assemblies
                       if re.fullmatch(rf".*/{sample}_{hap}_hprc_r2_v[0-9.]+\.fa\.gz", k))
        pick[hap] = found[-1] if found else None  # the latest version when several exist
    bam = f"{folder}analysis/aligned_reads/hifi/GRCh38/{sample}_aligned_GRCh38_winnowmap.sorted.bam"
    ok = pick["mat"] and pick["pat"] and bam in reads and bam + ".bai" in reads
    return sample, folder, pick, bam, reads.get(bam, 0), bool(ok)


def main() -> None:
    rows = [line.split() for line in fetch(POPULATIONS).splitlines()[1:] if line.strip()]
    population = {f[1]: (f[6], f[5], SEX.get(f[4], "NA")) for f in rows}
    superpopulation = {f[5]: f[6] for f in rows}
    for line in fetch(PEDIGREE).splitlines()[1:]:
        f = line.split("\t")
        if len(f) > 6 and f[1] not in population:
            population[f[1]] = (superpopulation.get(f[6], "NA"), f[6], SEX.get(f[4], "NA"))
    folders = []
    for top in ("working/HPRC/", "working/HPRC_PLUS/"):
        folders += listing(top, delimiter=True)[0]
    with ThreadPoolExecutor(16) as pool:
        found = list(pool.map(files, folders))
    print("\t".join(["sample", "folder", "superpopulation", "population", "sex", "bam_gb",
                     "rank", "bam_url", "pat_url", "mat_url"]))
    seen = set()
    for sample, folder, pick, bam, size, ok in sorted(found):
        if not ok or sample in seen:  # a sample in both folders is listed once, from HPRC
            continue
        seen.add(sample)
        superpop, pop, sex = population.get(sample, ("NA", "NA", "NA"))
        rank = hashlib.sha256(sample.encode()).hexdigest()[:16]
        print("\t".join([sample, folder, superpop, pop, sex, f"{size / 1e9:.1f}", rank,
                         BUCKET + bam, BUCKET + pick["pat"], BUCKET + pick["mat"]]))
    print(f"{len(folders)} folders, {len(seen)} candidates", file=sys.stderr)


if __name__ == "__main__":
    main()
