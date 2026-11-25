# Q6: phase relations from a single-sample assembly, where read-based phasers fail

Written 17 October 2026, before the assembly ran. Development gate, not
registered. Direction Y6 of the October directions memo (development repository), which was
infeasible on the laptop and is feasible on the server (1 TB of memory).

## Question

Read-based phasers align reads to the reference and then phase. In the hard
groups (`segdup`, `lowmap`, `tandem`, `indel`) their errors are often
supported by the reads themselves (W4), which no read-based score can see. A
diploid assembly places reads by overlaps between reads, not by the
reference. Are the phase relations of a single-sample assembly (HiFi only,
no parents, no Hi-C) right where the read-based phasers are wrong, and does
a hybrid that takes the assembly's relation where the read-based
certificate abstains certify more of the hard groups?

## Material

- Reads: HG002 HPRC HiFi, `data/w4/reads.fq.gz` (all reads, the same reads
  as W4's alignments). The read-based side is W4's `w4_grch38` platform
  (minimap2, 30-fold slice), whatshap and LongPhase as in W4.
- Assembly: hifiasm 0.25.0-r726 (built from source in `~/.local/opt/hifiasm`),
  default options, no trio or Hi-C data, output `hap1` and `hap2` contigs.
  Its contigs carry no phase across contig ends, so a relation is
  assembly-phased only when both variants lie in one `hap1` contig and in
  one `hap2` contig, as aligned by dipcall 0.3 to the GRCh38 no-alt set.
- Truth: GIAB v5.0q (the same as W4). v5.0q comes from the Q100 assembly,
  which used HiFi, ultra-long ONT and the parents, so it may share an
  assembly's errors; the GIAB trio transmission truth
  (`data/giab_wg/truth/`, independent of any assembly) is a second truth on
  the relations it phases.

## Measurements

Consecutive pairs of v5.0q heterozygous variants inside its BED on the 21
confirmatory autosomes. Per group: pairs that each source phases, error of
each source on the pairs all three phase, and error of each read-based
phaser on the pairs the assembly phases and the phaser's certificate
(frozen score, thresholds from W4's GRCh38 rows calibrated per group on odd
chromosomes) does not trust.

## Gate

Fixed before the assembly finishes:

- **(a) the assembly is right where reads mislead:** on even chromosomes,
  among pairs that whatshap phases wrongly and the assembly phases, the
  assembly is right on at least 80 per cent, in each of `segdup`, `lowmap`
  and `tandem` that hold at least 20 such pairs;
- **(b) the assembly is accurate in the hard groups:** on pairs all three
  phase, the assembly's error per pair is at most a third of whatshap's in
  each hard group with at least 20 whatshap errors;
- **(c) the truth does not decide it:** (a) and (b) hold also against the
  trio transmission truth, on the pairs it phases.

Q6 passes when (a), (b) and (c) hold. A pass motivates a registered hybrid
certificate; a fail says that a HiFi-only assembly does not fix the hard
groups, which is worth stating since assemblies are the usual answer.

## Clarifications, 18 October, before scoring

Written after the assembly and dipcall finished and before any relation was
compared with a truth, to fix what the gate left open:

- The pairs are consecutive variants of the read-based input
  (`data/ext_w4_grch38/<chrom>/hetsnv.vcf.gz`, v5.0q's heterozygous calls in
  its BED); a source relates a pair when it phases both variants in one
  block (for the assembly: one `hap1` and one `hap2` contig cover each
  variant's reference span, both in dipcall's BED, `0|1` or `1|0` at the
  same alleles).
- (a) uses even chromosomes as stated; (b) uses all 21, since it involves no
  calibration. (b)'s hard groups are `segdup`, `lowmap`, `tandem` and
  `indel`, the groups the question names.
- A group with fewer than 20 qualifying pairs (a) or 20 whatshap errors (b)
  is not judged; (a) and (b) each need at least one judged group.
- The scorer checks that its own whatshap and LongPhase relations agree
  with W4's cached rows on every pair both hold.

## Code

`benchmarks/q6_assembly_phase.py`; `tests/test_q6_assembly_phase.py` checks
the contig rule on fixed examples.
