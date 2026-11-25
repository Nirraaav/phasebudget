# Pre-registered W1: per-group error control across ancestries, HPRC HiFi

Written on 14 October 2026. The commit that adds this file is the
registration; its push time is the third-party record. Creating this file is
what lets `benchmarks/genome_common.py` open a test genome's reads
(`hprc_<sample>` for any sample other than HG002).

## Why

Every registered test so far put the calibrated thresholds on genomes of
European or East Asian ancestry (HG001, HG005, the CEPH 1463 family). Round 4
of the October directions memo (development repository) lists this as the largest gap a
reader will see. The HPRC holds HiFi reads aligned to GRCh38 and trio-phased
assemblies for 45 genomes, most of them of African or admixed American
ancestry. This test asks whether thresholds calibrated on HG002 keep their
bound on each of them.

## Material seen before writing

Everything in the development notebook up to 14 October. The third test (X8,
`docs/design/third-test-plan.md`) is still running and has not been scored.
Of the HPRC data, only these have been read:

- the bucket listing, file names and sizes (`scripts/hprc_candidates.py`);
- the header and the alignments in one 50 kb window of HG002's HiFi BAM
  (chromosome 21, 20.00 to 20.05 Mb), to check its tags: no read groups,
  winnowmap and samtools in the program lines, no kinetics tags. HG002 is the
  calibration genome;
- the gate's inputs, running since 14 October, 20:40: dipcall on HG002's
  year 1 assemblies and on HG005's release 2 assemblies. The gate compares
  truths with truths and reads no phaser output. When this file was written
  dipcall had not finished, so no gate result existed. HG005 is therefore
  not a test genome.

## What is tested

The procedure is the third test's, unchanged: Learn then Test calibrated
within each group, Vovk's conditional construction applied to Learn then
Test's risk, not claimed as new. The question is whether its guarantee
carries over from HG002 to genomes of other ancestries, sequenced and aligned
by the same consortium pipeline.

## Frozen design

Everything in the third test's "Frozen design" holds unchanged (the score,
the groups `indel`, `segdup`, `lowmap`, `tandem`, `difficult`, `easy`, the
100 kb unit, the grid, `delta = 0.05`, `beta` in {3e-4, 1e-3, 3e-3}, the
population of consecutive pairs of phased variants that the truth phases,
the truth's own heterozygous calls with phase removed as phaser input, the
depths 30, 20 and 10 drawn as there, both phasers, the 21 confirmatory
autosomes), except as listed here.

| | calibration (`hprc_HG002`) | test (`hprc_<sample>`) |
|---|---|---|
| genome | HG002 | the 44 genomes below |
| reads | HPRC HiFi, winnowmap GRCh38 (`working/HPRC_PLUS/HG002/analysis/aligned_reads/hifi/GRCh38/`) | the same pipeline, `bam_url` in `w1-candidates.tsv` |
| phaser input and truth | GIAB v5.0q heterozygous SNVs and indels in its BED | dipcall 0.3 on the sample's release 2 paternal and maternal assemblies against the GRCh38 no-alt analysis set (`scripts/dipcall_truth.sh`); `0|1` and `1|0` calls inside dipcall's BED; the first allele is the paternal assembly's |

- LongPhase runs with `--pb` instead of `--ont`; its other options, and
  `whatshap phase`'s, are unchanged.
- **Test genomes.** Every row of `docs/design/w1-candidates.tsv` except
  HG005: 44 genomes (22 AFR, 16 AMR, 4 EAS, and HG02109 and NA21309, which
  neither 1000 Genomes table lists). The table was made on 14 October by
  `scripts/hprc_candidates.py` at this commit. Genomes are processed in the
  order of its `rank` column (the first 16 hex digits of SHA-256 of the
  sample name), so a run stopped early leaves a subset that did not depend on
  anything about the genomes. If not all 44 finish, the scored genomes are
  those with all 21 chromosomes done, and the count is reported.
- **Depth.** If a genome's BAM holds less than 30-fold, its 30-fold slice is
  every read, and its measured depth is reported with the scores.
- **Shift stated now.** Calibration and test differ in person, ancestry and
  in how the truth was built (v5.0q, assembly-based from several
  technologies, against dipcall on HiFi-based trio assemblies). They share
  the sequencing platform, the aligner and the pipeline.

## Hypotheses

A cell is (genome, depth, phaser, beta, group). In a cell, `n` test units
hold a trusted pair and `k` hold a wrong one; it violates when the exact
binomial `P(X >= k | n, beta)` is below the hypothesis's level. A group that
certifies nothing cannot violate and is reported as abstaining. Levels are
fixed by the counts registered here, whatever number of cells turns out to
have test pairs or genomes finish.

- **A1, the easy group carries over (primary).** 44 genomes x 3 depths x
  2 phasers x 3 betas = 792 cells, level 0.05/792. *Supported* if no cell
  violates.
- **A2, the difficult and indel groups carry over (primary).** Groups
  `segdup`, `lowmap`, `tandem`, `difficult`, `indel`: 3,960 cells, level
  0.05/3,960. *Supported* if no cell violates. The scorer reports how many
  certify something.
- **A3, usefulness (primary).** Per genome, group `easy`, `beta = 1e-3`,
  depths 30 and 20, both phasers: 4 cells with the second test's S4 rule (at
  least 99 per cent of the group's pairs trusted, and wrong trusted pairs at
  most half the wrong pairs when there are at least 20, else under 1 per
  10,000 trusted). *Supported* if at least 90 per cent of scored genomes
  succeed in all 4 cells.
- **A4, ancestry (primary).** Group `easy`, `beta = 1e-3`, depth 30, each
  phaser: the mean certified share of the AFR genomes minus that of the AMR
  and EAS genomes lies within plus or minus 0.02, by two one-sided Welch t
  tests at 0.05 each, genome as the unit. *Supported* if both phasers show
  equivalence. HG02109 and NA21309 are left out of A4.
- **A5 (descriptive).** Per genome and group at `beta = 1e-3`: certified
  share, realized error among trusted pairs, the share of pairs in groups
  that certify nothing, and the number of truth pairs. By super-population:
  heterozygous variants per Mb.
- **A6 (descriptive).** HG002's thresholds calibrated against the dipcall
  truth from the gate instead of v5.0q, applied to the same test cells, to
  show what the truth's construction does to the thresholds. No decision is
  drawn from it.

## Gate and scorer

- **Gate, truth against truth.** Round 4 stated the gate as at most 1
  phase disagreement in 10,000 pairs. The error bound is per 100 kb unit,
  and a unit holds about a hundred pairs, so that rule is too loose; it is
  replaced here, before any gate result was seen, by a rule in the bound's
  own unit. For HG005 (release 2 dipcall) and HG002 (year 1 dipcall), take
  the consecutive pairs of dipcall's heterozygous calls (the W1 phaser
  input) whose two variants the trio transmission truth also phases
  (`data/giab_hg005/truth/` and `data/giab_wg/truth/`, phase set `PATMAT`,
  built from the GIAB v4.2.1 trio benchmarks and independent of any
  assembly). Group them as A1 and A2 do, in 100 kb units. A unit disagrees
  when it holds a pair whose relative phase the two truths disagree on. A
  (group, beta) is *testable* when the one-sided 95 per cent exact binomial
  upper bound on the share of disagreeing units, pooled over the two
  genomes, is at most beta/3. A disagreement rate of beta/3 raises a cell's
  observed rate by at most a third of its bound, which the Bonferroni levels
  of A1 and A2 absorb.
- A1 and A2 count only testable cells, at their registered levels. A3 and
  A4 are scored only if (`easy`, 1e-3) is testable. If no (group, beta) is
  testable, a dated amendment states what replaces the truth, before any
  phaser output of a test genome is compared with a truth.
- The gate also reports, descriptively, the same comparison against v5.0q
  for HG002 and against GIAB's `hifiasm_v11_phasetransfer` phase for HG005,
  and the heterozygous calls each pair of truths does not share.
- `benchmarks/w1_score.py` decides A1 to A4 and reports A5 and A6, reading
  only stored fragments, phaser VCFs and truths. It imports the second
  test's scoring functions from `benchmarks/stratified_score.py` unchanged.
  It is committed, as a dated amendment, with a dry run on `hprc_HG002`
  (odd chromosomes calibrate, even ones stand in for a test genome), before
  any test genome's phasing is compared with its truth.
- Every threshold is scored at full precision from the run file, never from
  printed text. Every deviation is a dated amendment below, stating whether
  results had been seen.

## Software and machine

`benchmarks/confidence_run.py`, `benchmarks/genome_common.py`,
`scripts/dipcall_truth.sh` and `scripts/hprc_candidates.py` at the
registration commit, on the KCDH compute server, with the third test's Python
environment, WhatsHap 2.8, LongPhase 2.0.2, dipcall 0.3 and minimap2 2.31.
The number of genomes, chromosomes and phasing jobs run at once does not
change any output: each draw has a fixed seed, and both phasers run
single-threaded.

Cost, estimated now: the 45 BAMs hold 4.67 TB; the 21 autosomes are about
0.88 of that, about 4.1 TB, or 45 hours at 25 MB/s. dipcall takes about 1 to
2 hours per genome.

## Limits stated now

- The test truths come from assemblies built from the same HiFi reads that
  are phased, with phase from parental k-mers. A read-level error both the
  assembly and the phaser follow is invisible to this truth, so a supported
  A1 or A2 bounds the errors this truth can see.
- No genome of European ancestry other than the calibration genome, none of
  South Asian ancestry, and 4 of East Asian ancestry are in the candidate
  set: A4 compares AFR with AMR and EAS only.
- The cells are correlated (nested depths, shared reads between the
  phasers, shared calibration). Bonferroni levels stay valid under that
  dependence and are conservative.

## Amendments

**1. The test truth becomes a consensus of dipcall and 1000 Genomes
transmission (15 October 2026; the gate's result of 14 October seen; no test
genome's phasing compared with any truth).** The gate found that dipcall on
one genome's assemblies disagrees with trio transmission in about 0.6 per
cent of 100 kb units, so only (`difficult`, 3e-3) was testable
(development notebook, 14 October). Of the two repairs listed there, this amendment
is the truth repair. The calibration repair is a separate registration, W1b
(`docs/design/w1b-pooled-plan.md`), committed with this amendment. At this
commit the runs of HG002 and three test genomes (HG02559, HG00735, HG01243)
are part done, no rows of any test genome have been computed, and no
comparison of dipcall with the 1000 Genomes panel has been run.

- **Truth.** For a test genome whose father and mother are both in the 1000
  Genomes 3,202-sample pedigree (`20130606_g1k_3202_samples_ped_population.txt`,
  md5 bac176a443f941890ac7b51fb1cabe0d), the truth phases a variant when
  dipcall's paternal allele equals the paternal allele that transmission
  gives. Transmission uses the trio's genotypes, not the phase, in the 1000
  Genomes high-coverage phased panel (`20220422_3202_phased_SNV_INDEL_SV`,
  each file checked against the panel's `20220804_manifest.txt`): a record at
  the same position and alleles, the child heterozygous, and exactly one
  Mendelian assignment, so at least one parent homozygous. Anything else
  leaves the variant unphased by the truth. Code: `genome_common.transmitted`
  (tested in `tests/test_kg_transmission.py`), `kg_transmission` and
  `consensus_truth`; `scripts/kg_trios.py` subsets the panel to the trios.
  The population is the registered one restricted to pairs whose two
  variants the consensus phases. The genome-wide HG002 run used the same
  construction (transmission minus the sites where the hifiasm phase
  disagreed).
- **Genomes.** 39 of the 44 test genomes have both parents in the pedigree.
  HG02559, HG02486, HG01123, NA21309 and HG02109 do not and are not scored.
  The registered levels (0.05/792, 0.05/3,960) stay.
- **Gate, replacing the 14 October rule for this truth.** A consensus pair is
  wrong only where the transmission pair is wrong, so transmission's error
  rate bounds the consensus truth's, whatever the relation between dipcall's
  errors and the panel's. NA12878 is the one child here with 1000 Genomes
  parents (NA12891, NA12892) and an independent truth, the Platinum Pedigree
  truthset. Take the consecutive pairs of her Pedigree heterozygous calls
  inside its high-quality regions that both transmission and the Pedigree
  truth phase, grouped as A1 and A2 do. A unit is wrong when it holds a pair
  whose relative phase the two disagree on. A (group, beta) is testable when
  the one-sided 95 per cent exact upper bound on the share of wrong units is
  at most beta/3 (`benchmarks/w1_consensus_gate.py`). A1 and A2 count
  testable cells only; A3 and A4 are scored only if (`easy`, 1e-3) is
  testable, as registered.
- **Descriptive** (`w1_consensus_gate.py --describe`): the consensus of
  dipcall and GIAB's trio transmission against v5.0q for HG002 and against
  GIAB's hifiasm phase transfer for HG005; for each test genome, dipcall
  against 1000 Genomes transmission and the share of dipcall's variants the
  consensus phases. Reported with the scores, not used.
- **Limits stated now.** NA12878 is one genome of European ancestry; the
  panel's genotype accuracy in African and admixed American trios is assumed
  to be no worse. The consensus drops pairs where the two truths disagree,
  and a phaser error at such a pair is not counted. Indels written
  differently in dipcall and the panel get no transmission label, so the
  `indel` group shrinks. If HPRC's trio binning used the same parental short
  reads as the panel (not checked here), dipcall and transmission share an
  input; the gate's bound does not depend on their independence.
- **Order of scoring.** W1b's calibration reads 20 of these test genomes.
  W1 and W1b are scored in one session after their runs end: the input
  check, `w1_score.py`, then `w1b_score.py`, each run once. No W1b result is
  computed before W1's scorer has run.
