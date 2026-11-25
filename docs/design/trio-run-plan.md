# Pre-registered plan for the HG002 trio run on real reads

Written on 28 September 2026, before any read of HG003 or HG004 was
downloaded or read. The commit that adds this file fixes that. Material read
before writing:
- the simulated trio budget benchmark (T3, the development notebook of 27 September, and
  OQ-16);
- the genome-wide HiFi and v5.0q results;
- the GIAB file listing of the parents' data, which contains no reads.

This is track T3's real-data step, and next step 4 of the development notebook.

## Amendments

1. **28 September 2026, before any parental read was read.** The variant
   rule also drops sites where both parents are homozygous for the same
   allele, which a heterozygous child contradicts. The pedigree DP cannot
   satisfy them, and the transmission truth excludes them as inconsistent.
   On chromosome 20 there are 34 such sites among the 43,835 where both
   parents have a genotype.
2. **30 September 2026, after the chromosome-20 smoke test started and
   before any parental read of a confirmatory chromosome was read.**
   - *Why.* In pedigree mode WhatsHap 2.8 splits its cap across members:
     `max(1, max_coverage // len(family))`, 5 per member at the default 15
     (`whatshap/cli/phase.py`, line 498). The registered `thirds` arm uses
     exact selection per member, so the deployed configuration was
     missing.
   - *The new arm.* `whatshap_thirds` selects each member's reads with
     WhatsHap's own selector at `C // 3`, and is phased like every other
     arm. It runs at totals 15 and 9.
   - **P4.** `thirds` at 15 (exact selection, 5 per member) against
     `whatshap_thirds` at 15, by the exact sign test on the child's common
     pairs at level 0.05. Supported if significant in the exact
     selection's favour. P1 to P3 are unchanged.

3. **1 October 2026, after chromosome 16 was phased and before any other
   confirmatory chromosome was phased.**
   - *The deviation.* The code review of 1 October found that the driver
     gave the child an allele at every kept site, including sites that no
     selected read covers; these were phased by genotype alone. `whatshap
     phase` restricts the DP, the genotypes and the master block to
     positions that reads cover (WhatsHap 2.8 `cli/phase.py`, line 711).
   - *The fix.* The driver now does the same, as the plan's text already
     said ("as `whatshap phase` does").
   - *The re-run.* Chromosome 16 is re-phased from its cached parental
     fragments, and its earlier phasing is discarded unscored.

## Question

WhatsHap's pedigree DP tracks every active read of every member, so its
state space grows as `2^C` in the total coverage `C` across the trio.
`whatshap phase` splits a cap into equal thirds. On simulation, one exact
selection over the union of the three members' reads under a total of 9
(`shared`) matched thirds at a total of 15 on the child's switch error, with
a DP about 60 times faster. That was not the registered arm; this run is its
confirmatory test on real reads.

## Data

**Child.** HG002: the full-depth `realign` fragments of the genome-wide
HiFi run (`data/genome/chrN/fragments_full_realign.json.gz`), unchanged.

**Parents.** HG003 and HG004: the GIAB PacBio HiFi chemistry-2 alignments
to GRCh38,
`.../HG003_NA24149_father/PacBio_CCS_15kb_20kb_chemistry2/GRCh38/HG003.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam`
and the matching HG004 file. Each is streamed one chromosome at a time and
deleted after extraction. Fragments come from WhatsHap 2.8's reader, with
realignment against the reference, at mapping quality 20, and are converted
to log-odds, as for the child.

**Chromosomes.** 16, 18, 21 and 22, chosen for size, all confirmatory.
Chromosome 20 runs first as a smoke test and enters no test.

**Variants.** The child's v4.2.1 heterozygous SNVs of the genome-wide run,
restricted to sites where both parents lie inside their v4.2.1 benchmark
BEDs. A parent's genotype is its v4.2.1 call, or homozygous reference
inside its BED where it has no call. Sites where a parent's call is not a
biallelic SNV at the same alleles are dropped, and every drop is counted.

## Allocations

For a total `C`, reads selected by exact min-cost flow:
- `thirds`: each member separately at `C // 3`, as in the simulation;
- `shared`: one selection over the union of the three members' reads under
  one cap `C` on their summed span coverage.

The runs are `thirds` at 15 and 9, and `shared` at 9 and 15. Phasing is
WhatsHap's `PedigreeDPTable` with a uniform recombination map of 1.26
cM/Mb, WhatsHap's default, as in `benchmark_pedigree.py`. Allele qualities are the rounded log-odds
weights. Blocks are WhatsHap's `find_components`, with every site at which
some member is homozygous as the master block, as `whatshap phase` does in
pedigree mode.

## Truth and statistics

The child's phase against GIAB v5.0q, assigned as in amendment 2 of
`genome-run-plan.md`. The statistics are those of `genome-run-plan.md` and
its amendment 1: common pairs, the exact sign test, and the difference per
10,000 common pairs with a 1 Mb window bootstrap pooled over the four
chromosomes. The DP's wall time and the largest number of active reads at
one variant are recorded per run.

## Hypotheses and decision rules

**P1, non-inferiority of the shared budget.** `shared` at 9 against
`thirds` at 15. Supported if both hold:
- the upper end of the 95 per cent interval of the difference (`shared` 9
  minus `thirds` 15) is below +0.5 switch errors per 10,000 common pairs;
- the pooled DP time of `shared` at 9 is at most one tenth of `thirds` at
  15.

**P2, shared against thirds at the same total.** `shared` at 9 against
`thirds` at 9: the exact sign test at level 0.05. Supported if significant
in `shared`'s favour.

**P3 (descriptive).** Every run: switch error, flips, phased fraction,
block N50, DP time, and the largest active read count.

**Power.** On HiFi the child's switch error in pedigree mode is likely a
few per 100,000 pairs, since transmission phases most sites. P2 may
therefore have few discordant pairs. The interval in P1 is the informative
result, and a wide interval is reported as undecided, not as support.

## Safety and resources

Each parent's slice of one chromosome is streamed, extracted and deleted
before the next is fetched, with at least 20 GB free before a download. The
four chromosomes total about 300 Mb. At an estimated 0.035 GB per Mb for 35-fold
HiFi and two parents, that is about 21 GB in all, but never more than one slice on
disk at once. The run starts after the ONT run has released its disk.

## Not in this plan

The weighted shared arm, which was worse on simulation; ONT parents; any
chromosome beyond the four.
