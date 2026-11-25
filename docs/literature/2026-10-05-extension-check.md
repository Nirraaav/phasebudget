# Literature check for the extension, 5 October 2026

Section 6 of the extension plan (development repository) listed four items to check
before registration. This note records what was read and what it changes.
Rows are in section I2 of `01-prior-art-matrix.md`.

"Full text" means the arXiv PDF was converted with `pdftotext` and the
passages were read directly. Abstracts were read on the arXiv abstract pages.

## 1. Mondrian (category-conditional) conformal prediction

Vovk 2012, "Conditional validity of inductive conformal predictors"
(arXiv 1209.2673v2), full text, sections 1, 3, 4 and 7.

- Section 4 defines the conditional inductive conformal predictor: a fixed
  inductive taxonomy `K` assigns each example a category, and the p-value is
  computed only from calibration examples of the test example's category
  (equation 10). Proposition 3: under exchangeability, the error probability
  given the category is at most `epsilon`.
- Section 3, Propositions 2a and 2b: the guarantee also holds conditionally on
  the training set, "of the PAC type", with probability `1 - delta` over the
  calibration set. Proposition 2b's condition is a binomial tail.
- Near the end of section 7 the paper says the label-conditional predictor
  "also enjoys the training conditional version" of its guarantee.

**Bearing.** Track A calibrates Learn then Test separately within each stratum,
a category fixed before calibration, with an exact binomial test. That is
Vovk's conditional construction applied to Learn then Test's risk. The
procedure is not new, and the report must not present it as new. What track A
can contribute is the measurement: whether per-stratum calibration fixes the
cross-genome failure, and how much each stratum must abstain.

## 2. Group-conditional risk control

- Gibbs, Cherian and Candès (arXiv 2305.12616, abstract): finite-sample
  coverage over each of a finite collection of subgroups, which may overlap,
  obtained by treating conditional coverage as coverage over covariate
  shifts.
- Blot, Angelopoulos, Jordan and Brunel (arXiv 2406.17819, abstract):
  approximate conditional control of expected risk, with the conditioning
  class chosen from the data.

**Bearing.** Both target overlapping or data-chosen groups and, for Blot et
al., risk in expectation. Track A's strata are disjoint by construction (the
first stratum that holds), so the simpler per-category construction of
section 1 suffices. If the strata were made to overlap, Gibbs et al. is the
method to use.

## 3. Covariate shift

Tibshirani, Barber, Candès and Ramdas (arXiv 1904.06019, abstract): weighted
conformal prediction when the likelihood ratio between test and training
covariate distributions is known or can be estimated.

**Bearing.** A test genome with a different mix of strata shifts only the
stratum label. A per-stratum statement does not depend on that mix, so it
needs no reweighting. A statement about all pairs pooled would need it. Shift
within a stratum (another platform, more paralogous copies) is not covariate
shift with a known ratio, and none of these papers covers it.

## 4. Published per-region phasing error guarantees

Europe PMC, 5 October 2026, `resultType=lite`, titles read:
- `(phasing OR haplotype) AND (conformal OR "risk control" OR "learn then
  test")`: 278 hits. The only sequencing result is SAVANA (already row I4),
  which uses Mondrian conformal prediction for somatic structural variants,
  not phasing.
- `(phasing OR "phase set") AND ("difficult regions" OR stratification OR
  "segmental duplication") AND ("error rate" OR confidence) AND (guarantee OR
  calibrat*)`: 388 hits. The first 25 by relevance contain no phasing error
  guarantee. One is the 2026 Cell paper "A complete diploid human genome
  benchmark for personalized genomics", the HG002 Q100 truth, read on
  1 October.
- `mondrian AND (variant OR genotype OR haplotype OR phasing)`: 130 hits; the
  only genomic method is SAVANA again.

These are searches, not proof of absence. No deliverable may rest on them.

## What changes

- The extension plan's track A is an application of a known construction.
  The stratified run plan should cite Vovk 2012 for the procedure and claim
  only the measurement.
- The 2026 Cell paper (Hansen et al., DOI 10.1016/j.cell.2026.06.016) was
  already read on 1 October (preprint full text in part, PMC12458380; Cell
  abstract; `2026-10-01-followup.md` section 5). GIAB's v5.0q README says
  v5.0q is built from the Q100 v1.1 assembly. An earlier version of this
  note called that "recalled"; it is checked. Still open before
  registration: whether v5.0q's BED leaves out difficult regions the
  assembly covers, which bounds calibration in the segmental duplication
  stratum.
