"""Track D development: relations between any two SNVs of one gene, `docs/design/chain-relations.md`.

For HG002 ONT hac and HG001 ONT sup at 30x, `whatshap phase` and LongPhase,
each chromosome is reduced, in a subprocess with the dataset's environment,
to its phased SNVs in position order with block, label (phaser allele XOR
truth allele, None where the truth has no phase) and the frozen score of the
link to the next phased SNV. Thresholds are the registered pooled ones,
Learn then Test on HG002 ONT sup rows. Relations are counted inside each
protein-coding gene of GENCODE v50 basic, within maximal runs of one block,
by `phasebudget.confidence.chain_counts`.

With `--clinical` (Y2's gate, `docs/design/family-and-clinical-gates.md`)
only genes with at least 10 ClinVar pathogenic or likely pathogenic alleles
count, and with `--stratified` a link is trusted when its score reaches its
stratum's threshold from `learn_then_test_by_group` on the same calibration
rows; the 20 genes with the lowest trusted share are listed.

The data were seen before: this is development, not a test.

Usage: `PYTHONPATH=src:benchmarks .venv/bin/python benchmarks/chain_relations.py`
"""

from __future__ import annotations

import bisect
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import confidence_score as cs  # noqa: E402

GENCODE = HERE.parent / "data" / "gencode" / "gencode.v50.basic.annotation.gtf.gz"
CLINVAR = HERE.parent / "data" / "clinvar" / "gene_specific_summary.txt"
MIN_PATHOGENIC = 10
DATASETS = ("ont_hac", "hg001_ont")
DEPTH = "30"


def chains_for(phaser: str, chromosomes) -> dict:
    """`{chrom: [(pos1, block, label, link score or None)]}` in the current platform's environment."""
    import genome_common as common
    from phase_confidence import phasing_from_vcf
    from phasebudget.confidence import pair_support, score

    extractor = common._module("extract_fragments.py")
    sc = common._module("score_hg002.py")
    out = {}
    for chrom in chromosomes:
        variants = sc.variant_table(common.hetsnv_path(chrom))
        matrix, positions, _, _ = extractor.load_fragments(
            str(common.fragments_path(chrom, DEPTH, "realign")))
        truth = common.truth_at(chrom, variants, positions)
        vcf = common.chrom_dir(chrom) / "external" / f"{phaser}_{DEPTH}.vcf"
        vcf = vcf if vcf.exists() else vcf.with_name(vcf.name + ".gz")
        phasing = phasing_from_vcf(str(vcf), common.VCF_SAMPLE, chrom, positions)
        link = {s.first: list(score(s)) for s in pair_support(matrix, phasing)}
        out[chrom] = [(positions[j] + 1, str(phasing.blocks[j]),
                       None if truth[j] is None else a ^ truth[j], link.get(j))
                      for j, a in enumerate(phasing.alleles) if a is not None]
    return out


def load_chains(platform: str, phaser: str, chromosomes) -> dict:
    path = cs.CACHE / f"chains_{platform}_{DEPTH}_{phaser}.json.gz"
    if not path.exists():
        env = dict(os.environ, PHASEBUDGET_PLATFORM=platform,
                   PYTHONPATH=f"{HERE.parent / 'src'}:{HERE}")
        subprocess.run([sys.executable, __file__, "--chains", phaser, str(path), *chromosomes],
                       env=env, check=True)
    with gzip.open(path, "rt") as handle:
        return json.load(handle)


def clinical_symbols() -> set:
    """ClinVar symbols with at least `MIN_PATHOGENIC` pathogenic or likely pathogenic alleles."""
    out = set()
    with open(CLINVAR) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if not line.startswith("#") and f[5] != "-" and int(f[5]) >= MIN_PATHOGENIC:
                out.add(f[0])
    return out


def genes(chromosomes, names=None) -> dict:
    """`{chrom: [(start, end, name)]}`, 1-based inclusive, protein-coding genes, optionally named."""
    import re

    out: dict = {c: [] for c in chromosomes}
    with gzip.open(GENCODE, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] == "gene" and f[0] in out and 'gene_type "protein_coding"' in f[8]:
                name = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
                if names is None or name in names:
                    out[f[0]].append((int(f[3]), int(f[4]), name))
    return out


def count(chain: list, trusted):
    """Totals over the given SNVs: all relations, then those whose chain the truth phases."""
    from phasebudget.confidence import chain_counts

    def runs(rows, cut):
        start = 0
        for k in range(1, len(rows) + 1):
            if k == len(rows) or cut(rows[k - 1], rows[k]):
                yield rows[start:k]
                start = k

    def add(total, rows):
        c = chain_counts([trusted(a, b) for a, b in zip(rows, rows[1:])], [r[2] for r in rows],
                         [r[0] // cs.UNIT for r in rows])
        for i, v in enumerate((c.pairs, c.wrong, c.trusted, c.wrong_trusted, c.trusted_span)):
            total[i] += v

    every, covered = [0] * 5, [0] * 5
    for block in runs(chain, lambda a, b: a[1] != b[1] or a[3] is None):
        add(every, block)
        for run in runs(block, lambda a, b: a[2] is None or b[2] is None):
            if run[0][2] is not None:
                add(covered, run)
    return every, covered


def main() -> None:
    import argparse

    import genome_common as common
    from phasebudget.confidence import (Units, coverage_grid, learn_then_test,
                                        learn_then_test_by_group)

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", action="store_true")
    parser.add_argument("--stratified", action="store_true")
    parser.add_argument("--per-gene", type=Path,
                        help="write {dataset phaser: {gene: [relations, trusted]}} at beta 1e-3")
    args = parser.parse_args()
    dump: dict = {}
    chroms = list(common.CONFIRMATORY)
    names = clinical_symbols() if args.clinical else None
    print(f"track D development, docs/design/chain-relations.md (data seen): {DEPTH}x; relations "
          f"between truth-phased SNVs in one block and one protein-coding gene of {GENCODE.name}"
          + (f", restricted to {len(names)} symbols of {CLINVAR.name} with at least "
             f"{MIN_PATHOGENIC} P/LP alleles" if names else "")
          + f"; {'per-stratum' if args.stratified else 'pooled'} thresholds by Learn then Test on "
          f"HG002 ONT sup rows, units {cs.UNIT // 1000} kb, delta {cs.DELTA}, grid from 0.5 in "
          f"steps of 0.0025; m = bins spanned by the endpoints; chromosomes {' '.join(chroms)}")
    gene_map = genes(chroms, names)
    print(f"protein-coding genes: {sum(len(v) for v in gene_map.values())}")
    stratum = None
    if args.stratified:
        from stratified_gate import strata

        stratum = strata(chroms)
    for phaser in cs.PHASERS:
        cal = cs.load_rows("ont", DEPTH, phaser, chroms)
        if args.stratified:
            rows = [(stratum(c, p1, p2), (c, p1 // cs.UNIT), k, w) for c, p1, k, w, p2, _ in cal]
            thresholds = {b: learn_then_test_by_group(rows, b, cs.DELTA) for b in cs.BETAS}
        else:
            units = Units.from_pairs([((c, p // cs.UNIT), k, w) for c, p, k, w, *_ in cal])
            grid = coverage_grid([r[2] for r in cal], start=0.5, step=0.0025)
            thresholds = {b: learn_then_test(units, grid, b, cs.DELTA) for b in cs.BETAS}
        for platform in DATASETS:
            chains = load_chains(platform, phaser, chroms)
            print(f"\n{platform} {phaser}")
            print("  beta    population  relations  wrong      trusted share  wrong trusted  "
                  "rate      mean bound beta*m")
            for beta, t in thresholds.items():
                if t is None:
                    print(f"  {beta:<7g} no threshold")
                    continue
                every, covered, per_gene = [0] * 5, [0] * 5, {}
                for chrom, rows in chains.items():
                    if args.stratified:
                        def trusted(a, b, chrom=chrom):
                            s = t.get(stratum(chrom, a[0], b[0]))
                            return s is not None and tuple(a[3]) >= s
                    else:
                        def trusted(a, b):
                            return tuple(a[3]) >= t
                    pos = [r[0] for r in rows]
                    for start, end, name in gene_map[chrom]:
                        lo, hi = bisect.bisect_left(pos, start), bisect.bisect_right(pos, end)
                        if hi - lo < 2:
                            continue
                        e, c = count(rows[lo:hi], trusted)
                        every = [a + b for a, b in zip(every, e)]
                        covered = [a + b for a, b in zip(covered, c)]
                        g = per_gene.setdefault(name, [0, 0])
                        g[0] += e[0]
                        g[1] += e[2]
                if args.stratified:
                    print(f"  {beta:<7g} thresholds: "
                          + ", ".join(f"{k} {'none' if v is None else (round(v[0], 6), v[1])}"
                                      for k, v in sorted(t.items())))
                for name, (n, w, tr, wt, span) in (("all", every), ("covered", covered)):
                    print(f"  {beta:<7g} {name:<10}  {n:>9}  {w:>9}  {tr / max(1, n):>13.4f}  "
                          f"{wt:>13}  {wt / max(1, tr):<8.2e}  {beta * span / max(1, tr):.2e}")
                if beta == 1e-3:
                    dump[f"{platform} {phaser}"] = per_gene
                if args.clinical and beta == 1e-3:
                    low = sorted((tr / n, g, n) for g, (n, tr) in per_gene.items() if n >= 100)[:20]
                    print("  lowest trusted share, genes with >= 100 relations: "
                          + ", ".join(f"{g} {share:.2f} ({n})" for share, g, n in low))
    if args.per_gene:
        args.per_gene.write_text(json.dumps(dump))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--chains":
        phaser, path, *chromosomes = sys.argv[2:]
        data = chains_for(phaser, chromosomes)
        with gzip.open(path + ".part", "wt") as handle:
            json.dump(data, handle)
        os.replace(path + ".part", path)
    else:
        main()
