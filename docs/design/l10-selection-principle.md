# L10: certificates hold on the labels' population (a synthesis, post hoc)

Round 8, 22 November 2026. Writing and re-analysis of archived outputs only;
no new run. Every number below is read from the development notebook or from the named
file in `benchmarks/results/runs/`. A cell marked "not found" could not be
read from either.

## The principle

Let a Learn then Test threshold be calibrated on labelled pairs drawn from a
population `C`, and scored on test pairs from a population `T`. Write `T`
as the part `T_on` that looks like `C` (the same label rule would have
labelled it) and the part `T_off` that does not.

1. On `T_on`, from the same read batch, the certificate holds: violations
   occur at no more than the nominal rate.
2. On all of `T`, the realized unit error rate at the threshold is
   approximately `(1 - s) r_on + s r_off`, where `s` is the off-population
   share of trusted units and `r_off`, `r_on` are the unit error rates of the
   two parts. A violation needs `s (r_off - r_on)` large enough to push that
   mixture above beta; when `r_off <= r_on`, there is no violation however
   large `s` is.

Part 2 is arithmetic once `s`, `r_on` and `r_off` are known for the test
genome. The claim with content is that, across this project, every
violation sits where `s > 0` and `r_off > r_on`, and that no matched test
violates. The read-batch condition in part 1 is stated because the project
has failures caused by batch alone (see the counterexamples section).

## Known statistics

This is the exchangeability assumption of split calibration seen from the
side of its failure. Learn then Test (Angelopoulos, Bates, Candès, Jordan
and Lei 2025, prior-art matrix row I2, status `paper`) guarantees risk
control for test points exchangeable with the calibration points; a test
population that the label rule did not sample is not exchangeable with it.
Nothing here is a new theorem. The rows on covariate shift (Tibshirani,
Barber, Candès and Ramdas; Gibbs, Cherian and Candès) are `abstract` status
and are not cited for support. What this note adds is an account of where
the assumption broke in phasing, and why: in each case the label rule
needs a second method to work, and that method fails at the same loci as
the phaser.

## Every certificate test in the project

"Matched" means the test pairs come from the label rule that chose the
calibration pairs. "Easier" means the test population is a subset that the
calibration's label rule would also label, selected towards easy links.

| Experiment | Calibration labels from | Test pairs from | Relation | Violating / cells | Source, the development notebook date |
|---|---|---|---|---|---|
| E1, HG002 hac | HG002 ONT sup, v5.0q BED | HG002 ONT hac, v5.0q BED | matched (other basecaller) | 0 / 18 | `2026-10-03_confidence_test.txt`, 4 Oct |
| E1, HG001 | HG002 ONT sup, v5.0q (GIAB v4.2.1 het SNVs in its BED) | HG001 ONT, Platinum 2017 confident regions | mismatched | 18 / 18 | same file, 4 Oct |
| E1 post hoc, HG001 in GIAB | as above | HG001 pairs whose SNVs are GIAB HG001 v4.2.1 het SNVs in its BED | matched | 0 / 18 | `2026-10-04_confidence_hg001_giab.txt`, 4 Oct |
| Track C2 | as above | HG001 against the Pedigree truth, adjudicable pairs (including outside GIAB) | mismatched | 12 / 18 | `2026-10-04_truth_adjudication.txt`, 4 Oct |
| X6 S1 | HG002 sup v5.0q and NA12878 ONT Pedigree truth | `easy` pairs, Pedigree members and HG005 | matched | 0 / 126 | `2026-10-13_stratified_score.txt`, 13 Oct |
| X6 S2 | NA12878 ONT, Pedigree truth | Pedigree members ONT, Pedigree truth, `difficult` and `indel` | matched | 0 / 540 (240 abstain) | same file, 13 Oct |
| X6 S3 | HG002 R10, v5.0q | HG005 R10, GIAB trio transmission (SNVs) | easier | 0 / 72 | same file, 13 Oct |
| X8 V1 | HG002 ONT sup, v5.0q | HG001 ONT sup, Pedigree truth, `easy` | mismatched regions, `easy` | 0 / 18 | `2026-10-14_third_score.txt`, 15 Oct |
| X8 V2 | as above | HG001, Pedigree truth, `difficult` groups and `indel` | mismatched (Pedigree regions reach further) | 2 / 90 (49 certify) | same file, 15 Oct |
| X8 V7 | HG002 R10 (older basecaller), v5.0q | as V1 and V2 | mismatched | 0 (smallest p 0.00224; cells with a unit not found) | same file, 15 Oct |
| W5 part 1 | HG002 ONT, pairs the Pedigree-style parental genotypes label (truthset genotypes) | all pairs of even chromosomes | mismatched, but labelled pairs not easier | 0 / 1,080 | The development notebook 14 Oct (output file not found in `runs/`) |
| W5 part 2, A and B | HG002 ONT, pairs labelled by one Clair3 replicate of the parents | all even pairs, v5.0q | mismatched, noisy labels | 0 / 108 each | `2026-10-15_w5_called.txt`, 15 Oct |
| W5 part 3, `agree` | pairs labelled by both replicates | all even pairs, v5.0q | mismatched | 1 / 108 (`lowmap`, 3e-3) | `2026-10-15_w5_noise.txt`, 15 Oct |
| W5 part 3, sparse filters | `positive`, `agree_positive`, `double` | all even pairs | mismatched | 10, 22, 12 / 108 | same file, 15 Oct |
| W5 part 3, oracle | all odd pairs, v5.0q | all even pairs | matched | 0 / 108 | same file, 15 Oct |
| R3 | label-eligible odd pairs (`agree`) | label-eligible even pairs | matched | 0 / 108 | `2026-10-16_r3_label_eligible.txt`, 16 Oct |
| R3b | HG002 label-eligible pairs | HG005 ONT label-eligible pairs (HG006, HG007 labels) | matched | 0 / 90 | `2026-10-27_r3b_score.txt`, 27 Oct |
| W1 A1, A2 | HG002 HiFi, v5.0q | 39 HPRC genomes, consensus key (dipcall equal to 1000 Genomes transmission), testable cells | easier | 0 / 234 and 0 / 468 | `2026-10-18_w1_score.txt`, 18 Oct |
| W1b B2 | HG002 and 20 HPRC genomes, consensus key | 19 held-out HPRC genomes, consensus key, testable cells | matched | 0 / 342 | `2026-10-18_w1b_score.txt`, 19 Oct |
| W1b B1 | as B2 | five Pedigree children HiFi, full Pedigree truth | mismatched | 163 / 540 (475 certify) | same file, 19 Oct |
| L1 step B, all links | as B2 | same children, all links | mismatched (reproduction) | 163 / 540 | `2026-11-22_l1_key_selection.txt`, 22 Nov |
| L1 step B, key links | as B2 | children, links their assembly-and-transmission key labels | matched | 0 / 540 | same file, 22 Nov |
| L1 step B, random | as B2 | children, same share of links kept at random | mismatched | 120 / 540 | same file, 22 Nov |
| Q3 | HG002 HiFi, inside v5.0q BED, Q100 truth | same genome, outside the BED | mismatched | 11 / 25 (whatshap and LongPhase); chr20 3 / 25 | `2026-10-17_q3_outside_benchmark.txt`, 18 Oct |
| K1a | HG002 HiFi odd chromosomes, pairs the assembly agrees with | even chromosomes, agreeing pairs | matched | 0 / 26 | `2026-10-19_k1_assembly_agreement.txt`, 19 Oct |
| K1b, agreement | HG002 and five Pedigree members, agreeing pairs | NA12877 and HG005, agreeing pairs | matched | 0 / 62 with a trusted unit (72 registered) | `2026-10-26_k1b_score.txt`, 26 Oct |
| K1b-4, frozen score | same pool, all pairs | NA12877 and HG005, all pairs | matched | 0 / 50 with a trusted unit | same file, 26 Oct |
| K4b | K1b pool, within-gene calls with end-to-end chain agreement | NA12877, HG005, same rule | matched | 0 / 8 | `2026-10-27_k4b_score.txt`, 27 Oct |
| L2 | K1b pool, end-variant agreement | NA12877, HG005, same rule | matched | 0 / 8 | `2026-11-22_l2_end_calls.txt`, 22 Nov |
| L7 P1 (registered) | W1b pool (HG002 and 20 HPRC genomes, consensus key) | NA12878 HiFi, Pedigree truth, all links | mismatched | 31 / 108 | `2026-11-24_l7_score.txt`, 24 Nov |
| L7 P2 (registered) | as P1 | NA12878, links her assembly-and-transmission key labels | matched | 0 / 108 | same file, 24 Nov |
| L7 P2, random | as P1 | NA12878, key's share of links kept at random | mismatched | 20 / 108 | same file, 24 Nov |
| L7 P3 (registered) | K1b pool, agreeing pairs | NA12878 30x, agreeing pairs | matched | 0 / 36 | same file, 24 Nov |
| L7 P4 (registered) | L2 pool, end-variant agreement | NA12878, same rule | matched | 0 / 4 | same file, 24 Nov |
| L8, all links (registered) | W1b pool | 19 held-out HPRC genomes, dipcall truth alone, all links | mismatched | 1,221 / 2,052 | `2026-11-24_l8_test.txt`, 24 Nov |
| L8, key links (registered) | W1b pool | same genomes, consensus-key links | matched | 0 / 2,052 | same file, 24 Nov |
| L8, random | W1b pool | same genomes, key's share kept at random | mismatched | 1,159 / 2,052 | same file, 24 Nov |

### Counterexamples and cases outside the principle

- **None among label-population cases.** Every violating row above is
  mismatched, and every matched row has 0 violations.
- **Mismatched without violation.** E1's restriction, V1, V7, W5 parts 1
  and 2, and the non-violating cells of the sparse filters: the principle
  allows these. W5
  part 1 is the informative one: its labelled pairs were not easier
  (NA12879 30x `segdup` 1.04e-2 labelled against 6.03e-3 unlabelled,
  whatshap), so `r_off <= r_on` and no violation is predicted. W5 part 2's
  noisy labels make calibration stricter, which also cannot violate.
- **Failures from read batch, not labels.** W3 pooled (15 Oct): 23 of 1,512
  cells violate, every one with HG002 held out, whose 2025 ONT reads differ
  from the six older genomes in the pool. R2 (16 Oct): the six-genome pool
  alone violates on HG002 (`easy` 2 cells, `tandem` 4 cells at 1e-3). These
  are not label selection, and they are why part 1 of the principle names
  the read batch. They would be counterexamples to a broader reading
  ("matched labels suffice").
- **A mixed case.** K1b's second dry run (26 Oct): the frozen-score
  certificate, calibrated on HG002 odd chromosomes (v5.0q) and five
  Pedigree genomes (Pedigree truth), violates on HG002's even chromosomes
  in `tandem` at 1e-3 and 3e-3 for both phasers. The pool mixes two truths
  and two sequencing batches, so it does not separate label population from
  batch. Counted as unresolved.
- **X8 V2 is consistent but not decomposed.** The violating `lowmap` cells
  sit in the Pedigree truth's wider regions, but no run split HG001's
  `lowmap` pairs into on- and off-population parts, so `s` and `r_off` are
  unknown there. The V2 post hoc (paralogous reads) is another explanation
  and is not excluded.

## Quantitative check

Within one genome, part 2 is an identity: with L1 step A's HG002 counts
(30x LongPhase `segdup`: key share 0.4692, key 5 / 27,367 = 1.83e-4, rest
163 / 30,955 = 5.27e-3) the mixture gives the all-links rate 168 / 58,322 =
2.88e-3 exactly. It cannot fail, so it is no test.

The test with content asks whether the off-to-on ratio measured on one
genome predicts the all-links rate of another. Prediction for the five
Pedigree children at 30x: their key-link pair rate (L1 step C) times
`s_c + (1 - s_c) R`, where `R` is HG002's rest-over-key ratio (step A) and
`s_c` the children's key share (step B pair counts, beta 3e-4 cells).
Observed: the children's all-pair rate against the Pedigree truth, from
L6's HiFi rows. Script: re-analysis of the two output files, no new rows.

| phaser, group | children key share | HG002 `R` | predicted rate | observed rate | predicted / observed |
|---|---|---|---|---|---|
| whatshap `segdup` | 0.464 | undefined (0 key errors) | none | 8.60e-3 (4,889 / 568,627) | none |
| whatshap `lowmap` | 0.587 | 53.8 | 1.92e-3 | 1.60e-3 (969 / 606,256) | 1.20 |
| whatshap `tandem` | 0.651 | 142.8 | 2.19e-2 | 1.57e-3 (1,739 / 1,107,165) | 13.9 |
| LongPhase `segdup` | 0.463 | 28.8 | 2.42e-3 | 5.49e-3 (3,123 / 568,497) | 0.44 |
| LongPhase `lowmap` | 0.591 | 9.0 | 1.80e-4 | 7.57e-4 (455 / 601,454) | 0.24 |
| LongPhase `tandem` | 0.651 | 25.9 | 1.04e-3 | 3.41e-4 (377 / 1,106,193) | 3.04 |

The direction transfers (off-key links are worse in every row), but the
size does not: predictions miss by factors 0.24 to 13.9. The two keys
differ (HG002: dipcall equal to GIAB trio transmission; children:
hifiasm with transmission from 1000 Genomes parents), and whatshap's
`tandem` errors on HG002 are 4 times the children's. So the principle
predicts where violations occur, not how large they are on a new genome
from another genome's ratio. Measuring `r_off` needs a full truth on a
sample of the test genome's own off-population links.

## What would falsify it

- A matched test (same label rule at calibration and test, same read
  batch) with a violating cell at its registered level.
- A mismatched test with a violation where the off-population links are no
  more error-prone than the on-population ones, measured on the same
  genome.
- A key-restricted rescoring (as L1 step B) that leaves as many violations
  as a random restriction of the same size.

L7 (registered 22 Nov, `docs/design/l7-na12878-plan.md`) tests the first
and third on NA12878 HiFi: P1 predicts W1b's thresholds violate on all
links, P2 that they hold on her key links, P3 and P4 that matched
agreement certificates hold.

## Status

Post hoc. The rows were collected after the results, and the principle was
written to fit them. The rows that were registered tests and agree with it
are E1, X6 S1 to S3, X8 V1, V2 and V7, W1 A1 and A2, W1b B1 and B2, K1b,
K4b and R3b; in none of them was the principle stated before scoring. L1,
L2, Q3, K1a, R3 and W5 are development results. The first test registered
with this principle stated in advance is L7.
