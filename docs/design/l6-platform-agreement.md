# L6: HiFi and Nanopore agreement as a second opinion (development, 22 Nov 2026)

## Question

K1b showed that a pair a phaser relates is almost never wrong when a
single-sample hifiasm assembly gives the same relation, and that pooled
certificates built on agreeing pairs are valid on two new genomes. An
assembly costs about 10 hours and 100 GB of memory per genome. Many labs
that have one person sequenced on two platforms have no assembly. L6 asks
whether a second phasing of the same person from the other platform can
take the assembly's place.

It can only do so if the two platforms make different mistakes. Both are
aligned to the same reference, so a duplication collapsed in GRCh38 can
fool both. Q5 (17 Oct) found that error loci transfer across platforms
within one person (HG002 ONT wrong-pair rate 4.6 times higher for whatshap
and 21.3 times for LongPhase near HiFi error bins), which is evidence
against L6. The plan put the odds at 45 to 55 per cent.

## Status of the data

Every input below has been read and scored before (W1b, X6, X3, K1b). This
is a development analysis. The rules are written here before any
cross-platform number is computed, and the commit of this file comes before
the first run of `benchmarks/l6_platform_agreement.py`.

## Data

People with both platforms phased by the same two phasers at 30x, all 21
confirmatory chromosomes (chr1 to chr22 without chr20):

| person | HiFi platform | ONT platform | truth (both platforms) |
|---|---|---|---|
| HG002 | `hprc_HG002` | `cal_v5` | GIAB v5.0q |
| NA12879, NA12881, NA12882, NA12885, NA12886 | `phifi_<s>` | `ped_<s>` | Platinum Pedigree truthset |
| HG005 | `k1_hg005` | `hg005_ont` | GIAB HG005 trio transmission |

NA12877 has HiFi only and is left out. Wrong flags come from each
platform's cached rows (`confidence_score.load_rows`), so each platform is
scored against its own primary truth, which is the same truth for both
platforms of one person.

## Computation

For a person, a phaser `ph` (whatshap or LongPhase) and a direction
(primary: HiFi checked by ONT; reported: ONT checked by HiFi):

1. Rows of the checked platform: every adjacent pair `(a, b)` the phaser
   relates and the truth phases, with its score and wrong flag.
2. The checking platform's value at each variant: its own `ph` 30x VCF's
   allele and block, looked up by `(chrom, pos, ref, alt)`. A variant the
   checking platform did not call, or left unphased, has no value. The
   checking relation is `q6_assembly_phase.relation` of the two values: it
   exists only when both variants sit in one block of the checking phasing.
3. `agree` is `k1_assembly_agreement.agrees(checked relation, checking
   relation)`, exactly the K1b rule with the second platform in place of
   the assembly.
4. Group: `indel` if either variant is not an SNV, else the stratum of
   `stratified_gate.strata` (`segdup`, `lowmap`, `tandem`, `difficult`,
   `easy`).

Per (person, phaser, direction, group) the output prints: pairs; matched
share (both variants called by the checking platform); relates share
(checking relation exists); agree share; wrong pairs; wrong pairs the
checking platform relates (`W`); shared wrong pairs (wrong and agree, so the
checking platform makes the same mistake, `S`); shared wrong share `S / W`;
wrong rate among agreeing pairs; wrong rate among all pairs.

Because the relation is binary, a wrong pair that the checking platform
relates is either shared (same wrong relation) or caught (opposite
relation). `S / W` is therefore the share of wrong pairs that agreement
fails to remove, among those it can judge.

## Gate (decides whether step 2 runs)

The gate passes when, for at least one phaser `ph`, in the HiFi-checked-by-
ONT direction and in group `segdup`:

- (a) HG002 has `S / W < 0.20`, and at least 3 of the 5 Pedigree members
  have `S / W < 0.20`. A cell with `W < 10` does not count as passing.
- (b) In each of those passing cells, the agree share of `segdup` pairs is
  at least 0.5, so a low shared share does not come only from rejecting
  most pairs.

HG005 and the ONT-checked-by-HiFi direction are reported and do not decide.
If the gate fails, L6 stops and is logged as a negative with these numbers.

## Step 2 (only if the gate passes)

K1b's certificate with platform agreement in place of assembly agreement.
`learn_then_test_by_group` on pooled agreeing HiFi pairs, unit (person,
chromosome, 100 kb), betas 3e-4, 1e-3, 3e-3, delta 0.05.

- Tests: HG005 against the pool of all six people; and each of the six pool
  people against the pool of the other five (leave one out).
- Cells: 7 tests x 2 phasers x 3 betas x 6 groups = 252; level 0.05/252.
- L6-2 validity: no agreeing-pair cell with `P(X >= k | n, beta)` below the
  level.
- L6-3 usefulness: `segdup` certifies at 1e-3 on the full pool for at least
  one phaser, and that phaser trusts at least 0.5 of HG005's `segdup` pairs
  at 1e-3.

These two rules are development rules on data already seen. A positive
result would need a registered test on a new person with both platforms
(for example a further Platinum Pedigree member or an HPRC genome with ONT
data) before it could be claimed.
