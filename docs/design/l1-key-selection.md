# L1: does the consensus key select easy links? (development, rules stated before scoring)

Written on 22 November 2026, before any of the numbers below were computed.
W1b's scores (development notebook, 19 October) are already known: B1 failed with
163 violating cells of 540 at level 9.26e-5, on five Platinum Pedigree
children scored against the Pedigree truth. This is a development analysis
of outputs that exist, with its decision rules written down first. It is not
a blind test.

## Claim

W1b calibrated on HG002 (v5.0q) and 20 HPRC genomes whose truth is W1
amendment 1's consensus key: dipcall's phase on the genome's own assemblies,
kept only where 1000 Genomes transmission gives the same paternal allele. A
link (consecutive phased pair) enters calibration only if the key labels
both variants. The claim is that the key labels the links a phaser gets
right far more often than the links it gets wrong, in the hard groups. The
pooled thresholds then control error on key-labelled links and are too loose
on all links, which is what a test against the complete Pedigree truth
measures.

The alternatives are dependence between genomes at one locus (the
registered risk) and a sequencing or pipeline batch effect. Step C describes
both; it decides nothing.

## Data

All rows are the cached `confidence_score.load_rows` rows already used by
W1 and W1b (`data/confidence_rows/`), groups by `stratified_gate.strata`
plus `indel`, unit (genome, chromosome, 100 kb).

## Step A (decisive): one genome, same reads, two keys

`hprc_HG002`, 30x, whatshap and LongPhase, 21 confirmatory chromosomes.
Wrong is against v5.0q, as in W1 and W1b. The HPRC-style key for HG002 is
the closest analogue of the consensus key that exists for a genome with no
1000 Genomes parents: dipcall on HG002's HPRC assemblies
(`data/hprc/HG002/HG002.dip.vcf.gz`, inside its BED), kept where GIAB's trio
transmission truth (`data/giab_wg/truth/HG002_<chrom>_trio_truth.vcf.gz`,
PATMAT calls) gives the same paternal allele. A link is key-labelled when
both variants are.

Per group: v5.0q-labelled links, share key-labelled, wrong rate on all
links, on key-labelled links and on the rest. The selection ratio of a set
of groups is the wrong rate on all its v5.0q-labelled links divided by the
wrong rate on its key-labelled links.

**Rule A.** Supported if the selection ratio over `segdup`, `lowmap` and
`tandem` together exceeds 3 for both phasers at 30x. A key-labelled wrong
count of 0 counts as a ratio above 3 only if the all-links count is at least
3. The per-group ratios and the worst one are printed; 20x and 10x are
printed and decide nothing.

Limit stated now: GIAB's transmission truth is restricted to its v4.2.1
benchmark regions, while the 1000 Genomes panel calls genome-wide. The HG002
key may therefore be stricter than the consensus key, which would inflate
the ratio. Step A prints the share labelled by dipcall alone and by
transmission alone, so the two filters can be read apart.

## Step B: the five children, rescored on links a key labels

Same pooled thresholds as `benchmarks/w1b_score.py` (HG002 on v5.0q plus the
20 consensus-truth genomes, per depth 30/20/10, phaser and beta in
{3e-4, 1e-3, 3e-3}), applied to `phifi_NA12879`, `NA12881`, `NA12882`,
`NA12885`, `NA12886`, wrong against the Pedigree truth as in B1. Level
0.05/540 = 9.26e-5, a cell violating when it has test units and p is below
the level, exactly W1b's `decide`.

The children have no 1000 Genomes parents in the panel, but both parents
(NA12877, NA12878) are panel members. The family key is dipcall on the
child's own hifiasm assembly (`data/k1/<child>/dip/`, inside its BED), kept
where transmission from the parents' panel genotypes (child heterozygous,
which every scored variant is) gives the same paternal allele. This is the
consensus key with the child's own short-read genotype replaced by the
truthset's heterozygous call.

Three scorings of the same cells:

1. **All links** (reproduction): must give 163 violating cells. If it does
   not, Step B stops and reports the difference.
2. **Key links**: only links the family key labels.
3. **Random control**: per (child, depth, phaser, group), links kept at
   random with probability equal to that cell's key-labelled share (seed 0),
   so the unit count falls as in scoring 2 without favouring easy links.

**Rule B.** Supported if scoring 2 has at most 32 violating cells (a cut of
more than 80 per cent from 163) and scoring 3 has at least 82 (half of 163
survive). If both scorings cut violations by more than 80 per cent, the drop
is attributed to lost power and Rule B fails.

## Step C (descriptive)

- **Locus sharing.** At 30x, per phaser and hard group: of each child's
  wrong units (by locus, chromosome and 100 kb index), the share also wrong
  in at least one other child, and the share wrong in at least one of the 21
  pool genomes; the same share for pool genomes against the rest of the pool.
- **Batch.** Per group at 30x, wrong rate on key-labelled links: the five
  children (family key labels, Pedigree truth) against the 20 HPRC pool
  genomes (consensus truth). Equal rates on labelled links point away from a
  batch effect.

## Outputs

`benchmarks/l1_key_selection.py`; output
`benchmarks/results/runs/2026-11-22_l1_key_selection.txt`, header printing
every parameter. Key labels are cached in `data/l1/keys/`.

## Amendment 1 (22 November 2026, before any scoring)

A smoke run of the key labels on chr21, which reads no phasing and no wrong
flag, showed that the family key as written cannot work. On NA12879, dipcall
on the child's hifiasm assembly and parental transmission disagree on the
paternal allele at 12,455 of 25,938 variants both label. The K1 assemblies
are single-sample hifiasm assemblies without trio binning, so dipcall's
first allele is hap1 of an arbitrary contig orientation, not the paternal
allele. HG002's HPRC assembly is trio-binned and shows 2 disagreements of
21,213 on the same chromosome, so Step A is unchanged.

The family key is redefined on links. A link is labelled when the child's
assembly relates both variants (the Q6 contig rule, `q6_assembly_phase.relation`
on the assembly values cached by K1b in `data/k1/tables/phifi_<child>_30_<chrom>.json`),
transmission labels both variants, and the two give the same relation (same
haplotype or opposite). For a trio-binned assembly this equals the consensus
key on that link, apart from links across contigs, which the contig rule
leaves unlabelled. Rules A and B are unchanged.
