# Pre-registered K1b: certificates that require a single-sample assembly to agree

Written on 19 October 2026. The commit that adds this file is the
registration; its push time is the third-party record. Creating this file is
what lets `benchmarks/genome_common.py` open the two test genomes' reads
(`phifi_NA12877`, `k1_hg005`). The author chose on 19 October to register K1b
although K1a did not pass by its letter (development notebook, 19 October).

## Why

K1a (`docs/design/k1-assembly-agreement.md`) found on HG002 HiFi that a
relation trusted only when the read score passes and a HiFi-only hifiasm
assembly phases the pair the same way removes nearly all hard-group errors:
0 wrong of about 53,000 agreeing `segdup` pairs, 1 of about 207,000 `tandem`
pairs, against v5.0q, and 1 of 101,838 agreeing `tandem` pairs against the
GIAB trio truth. `tandem` and `indel` certified at `beta = 1e-3`; `segdup`
could not, because one genome's calibration chromosomes hold about 1,400
`segdup` units and Learn then Test needs 2,995 error-free units at 1e-3. K1b
pools calibration genomes to pass that floor and tests on genomes the
calibration did not see.

## Material seen before writing

Everything in the development notebook up to 19 October, including K1a and K4. The pool's
Pedigree genomes were W1b test genomes on the same HiFi reads (scored 19
October) and second-test genomes on ONT. NA12877's reads of any platform have
not been read by this project. HG005 was a test genome before: HiFi reads
from GIAB's 11 kb Sequel II library in the HG005 replication (H5, 1 October)
and ONT in the second test (X6, 13 October); its HPRC HiFi reads (HudsonAlpha,
2020), used here, have not been read, and no assembly of HG005 has been made.

## Design

Unchanged from W1b unless listed: score, groups `indel`, `segdup`, `lowmap`,
`tandem`, `difficult`, `easy`, the 100 kb unit, the grid,
`learn_then_test_by_group`, `delta = 0.05`, `beta` in {3e-4, 1e-3, 3e-3}, both
phasers (WhatsHap 2.8 defaults; LongPhase 2.0.2 `--pb --indels -t 1`), the 21
confirmatory autosomes.

- **Depth.** 30-fold only, drawn as `confidence_run.py` draws it. The
  assembly uses all of a genome's reads, which a user has only at full depth.
- **Assembly.** hifiasm 0.25.0-r726, default options, `-t 24`, on all of a
  genome's HiFi reads, no parental or Hi-C data; `hap1` and `hap2` contigs
  aligned to the GRCh38 no-alt analysis set by dipcall 0.3 with
  `hs38.PAR.bed`. For Platinum Pedigree members the reads are the primary and
  unmapped records of the provider's GRCh38 BAM (`samtools fastq -F 0x900`,
  `tmp/k1b_asm.sh`); for HG005 the seven HudsonAlpha CCS FASTQ files of HPRC
  (`working/HPRC_PLUS/HG005/raw_data/PacBio_HiFi/*.fastq.gz`, checked against
  the bucket's `checksum.md5_fq`); for HG002 Q6's assembly.
- **Agreement.** Q6's contig rule: the assembly relates two variants when each
  variant's whole reference span lies in one `hap1` and one `hap2` contig
  (the same two for both), inside dipcall's BED, with a phased `0|1` or `1|0`
  call at the same alleles. A pair agrees when the assembly relates it with
  the phaser's relation. A pair that does not agree is never trusted.
- **Calibration pool.** Agreeing pairs of HG002 (`hprc_HG002`, truth v5.0q)
  and of the five Platinum Pedigree members NA12879, NA12881, NA12882,
  NA12885 and NA12886 (`phifi_<sample>`, truth the Platinum Pedigree truthset
  inside its high-quality regions), all 21 chromosomes, 30x. A unit is
  (genome, chromosome, 100 kb index).
- **Test genomes.** NA12877 (`phifi_NA12877`, Platinum Pedigree HiFi, truth
  the Pedigree truthset as above, heterozygous calls of the truthset as phaser
  input) and HG005 (`k1_hg005`, HPRC HiFi aligned by winnowmap, phaser input
  GIAB v4.2.1's heterozygous calls inside its benchmark BED, truth GIAB's trio
  transmission as in the second test).
- **Shift stated now.** NA12877 is the father of all five Pedigree pool
  genomes: he shares one haplotype with each at every locus, and the same
  project, library pipeline and truth. His test is within the family and the
  easiest the design allows. HG005 is unrelated (Han Chinese), sequenced by
  another centre in 2020, aligned by another pipeline, and judged by a truth
  that comes from short-read trio inheritance; that truth phases SNVs only, so
  HG005 has no `indel` cell, and it phases fewer hard-group pairs.

## Hypotheses

A cell is (test genome, phaser, beta, group): `n` test units hold a trusted
pair, `k` hold a wrong trusted pair, and the cell violates when the exact
binomial `P(X >= k | n, beta)` is below the level. A group that certifies
nothing, or has no truth-phased pair, has `n = 0` and cannot violate. There
are 2 x 2 x 3 x 6 = 72 registered cells; the level is 0.05/72 = 6.94e-4
whatever number end up with `n > 0`.

- **K1b-1, validity (primary).** No cell of either test genome violates.
- **K1b-2, the floor (primary).** At `beta = 1e-3` the pooled agreement
  calibration certifies `segdup` for both phasers. Decided on the calibration
  alone.
- **K1b-3, usefulness (primary).** On NA12877 at `beta = 1e-3`, for each of
  `segdup`, `tandem` and `indel`, at least one phaser trusts at least half of
  the group's pairs (all pairs the phaser relates and the truth phases,
  agreeing or not). On HG005 the same for `segdup` and `tandem`.
- **K1b-4 (descriptive).** The same pool's frozen-score certificate without
  agreement (W1b's procedure on this pool) on the same cells: violating cells
  and trusted shares, beside the agreement certificate's.
- **K1b-5 (descriptive).** Per test genome and group: the share of pairs the
  assembly relates, the share that agree, wrong pairs among agreeing and
  other pairs, and pairs wrong in both the phaser and the assembly.

## Risk stated now

Pooled Learn then Test treats units of different genomes at one locus as
independent. Inside a family that is optimistic: a parent and child share a
haplotype, so an assembly and a phaser can fail together at the same locus in
both. A K1b-1 failure confined to NA12877 would point here; one on HG005 would
point to batch or ancestry. The pool's Pedigree genomes come from the
project that also made NA12877's data, so K1b-1 on NA12877 does not test the
cross-project shift that broke W1b; HG005 does.

## Scorer and order

- `benchmarks/k1b_score.py` decides K1b-1 to K1b-3 and reports K1b-4 and
  K1b-5, reading stored fragments, phaser VCFs, dipcall outputs and truths,
  with Q6's contig rule (`q6_assembly_phase.assembly_values`) and the second
  test's cell function (`stratified_score.cells`) unchanged. It is committed
  as a dated amendment with a dry run that reads no test genome: HG002's odd
  chromosomes and the five Pedigree genomes form the pool, HG002's even
  chromosomes stand in for a test genome.
- Runs: the assemblies (`tmp/k1b_asm.sh`, started 19 October before this
  registration; they read no truth), `confidence_run.py` for `phifi_NA12877`
  and `k1_hg005`, then the scorer once.
- If a pool genome's assembly or run fails, the pool is the genomes that
  finished, and the count is reported. A test genome is scored when all 21
  chromosomes finished.
- Every threshold is compared at full precision from the run files. Every
  deviation is a dated amendment below, stating whether results had been
  seen.

## Software and machine

the KCDH compute server, the third test's Python environment, WhatsHap 2.8,
LongPhase 2.0.2, hifiasm 0.25.0-r726 (`~/.local/opt/hifiasm`), dipcall 0.3
(`~/.local/opt/dipcall.kit`), aria2 for the read downloads.

## Amendments

### Amendment 1, 19 October 2026: the scorer and its first dry run

Written before any test genome's run finished; no test genome was read.
`benchmarks/k1b_score.py` (committed at 33d53a1, with `q6_assembly_phase.dipcall_values`
factored out of Q6's table reader, behaviour unchanged) ran with `--dry-run`
(`benchmarks/results/runs/2026-10-19_k1b_dry_run.txt`, 1,171 s). No Pedigree
assembly had finished, so the pool was HG002's odd chromosomes alone, split in
two pseudo-genomes, and HG002's even chromosomes stood in for a test genome:
no agreement cell violates at 6.94e-4; `segdup` abstains at 1e-3 for both
phasers with 1,630 and 1,637 pooled units, below the 2,995 the floor needs,
which is what K1b-2 exists to test with the full pool. A second dry run with
every pool genome runs before the scorer is run on the test genomes.

### Amendment 2, 26 October 2026: the full-pool dry run, NA12885, and the score

Written after the score below; no rule changed. Before the test genomes were
read, the second dry run used every pool genome
(`benchmarks/results/runs/2026-10-26_k1b_dry_run_pool.txt`): `segdup`
certifies at 1e-3 for both phasers, and on HG002's even chromosomes no
agreement cell violates (smallest p 7.6e-3). The frozen-score certificate
of K1b-4 violates there in `tandem` at 1e-3 and 3e-3 for both phasers. NA12885's hifiasm log reports `peak_het: -1`; its haplotypes
are 3.01 and 3.02 Gb and its dipcall heterozygous to homozygous ratio is 1.72
(other genomes 1.57 to 1.92), so it was kept in the pool. The registered
scorer then ran once (`2026-10-26_k1b_score.txt`): K1b-1, K1b-2 and K1b-3 are
supported. The development notebook, 26 October, has the numbers.
