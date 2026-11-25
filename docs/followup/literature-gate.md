# Literature gate G-E: error-controlled phasing

Run on 2 October 2026 for section 6 of `error-controlled-phasing.md`. The
rows, with bibliographic fields and statuses, are in
`docs/literature/01-prior-art-matrix.md`, section I. Rows from sections F and
H (ProbHap, PhaseQuality, HiPhase, LongPhase, LongPhase 2, pHapCompass) are
cited from there. Only `paper`, `abstract` and `docs` rows support a
statement below; `recalled` rows are listed as unread.

## 1. What was searched and how it was read

- **Tool documentation**, fetched as raw files from each repository's default
  branch: HapCUT2 `outputformat.md` and `README.md`; WhatsHap `guide.rst`;
  HiPhase `user_guide.md` and `methods.md`; LongPhase `README.md` (v2.0.2);
  the GATK 3.7-0 ReadBackedPhasing page in `broadinstitute/gatk-docs`; the
  VCF v4.3 specification; the LongPhase 2 manuscript repository README.
- **Papers**, read from Europe PMC JATS XML, the PMC BioC service, or arXiv
  PDFs converted to text, and searched for the passages that bear on the
  question. Section and passage coverage is stated per row.
- **Europe PMC queries**, including preprints: conformal prediction with
  phasing or haplotype; conformal with variant calling or genotyping;
  conformal with genomics and variants; risk control with phasing; selective
  classification with genomics; phase or phasing quality or confidence with
  calibration; phase confidence with long-read phasers; switch error with
  abstention, reject option or selective; haplotype assembly or read-based
  phasing with abstention, conformal, risk control or selective prediction;
  phasing with false discovery rate or error control and reads;
  compound heterozygosity with long-read phasing in diagnosis; PM3 and in
  trans.
- **arXiv site search**: conformal with haplotype, phasing, or genotype;
  haplotype assembly uncertainty calibration. The arXiv API was
  rate-limited and was not used for search.
- **Web search**: conformal prediction with haplotype phasing and switch
  error; per-variant phasing confidence calibration; GIAB and GA4GH phasing
  benchmark guidance; matched-completeness comparison of phasers; conformal
  methods in genotype imputation and variant calling; the ClinGen PM3
  recommendation.

No personal information was sent to any service. Crossref was queried
without a contact address.

## 2. What each source provides

### 2.1 Phase confidence in phasers

- **HapCUT2** (Edge, Bafna and Bansal 2017, `paper`; documentation `docs`)
  has the richest output. Its block file carries, per SNV, a pruning flag, a
  phred-scaled switch quality and mismatch quality, and the number of
  informative fragments. The switch quality is computed only with
  `--error_analysis_mode 1`, and automatic block splitting is documented as
  broken. The default prunes at a fixed posterior of 0.8. The paper plots
  error against the share pruned, for HapCUT2's and ProbHap's scores on
  chromosome 1 of 11x PacBio, and compares tools on the variants all of them
  phased. It does not report whether the posteriors are calibrated. We did
  not find any later measurement of HapCUT2's calibration.
- **ProbHap** (Kuleshov 2014, section F, `paper`) gives a per-SNP transition
  score, the probability that the phase relative to the previous SNP is
  right, and uses it to split blocks.
- **WhatsHap** (`docs`) leaves a variant unphased when "the phasing algorithm
  could not make a phasing decision", and "currently does not add" the `PQ`
  tag. A phased WhatsHap VCF has no per-variant quality.
- **HiPhase** (`docs`) writes `GT` and `PS`. A variant is left unphased when
  treating it as homozygous gives a better solution, and `PF=TR_OVERLAP`
  marks variants excluded inside tandem-repeat calls. Its stats file has
  block-level solver diagnostics, not phase quality.
- **LongPhase** (`docs`) writes `GT` and `PS`. It gives up a SNP pair when the
  read counts for the two configurations are similar, under
  `--edgeThreshold` (default 0.7).
- **GATK ReadBackedPhasing** (`docs`) outputs phase only when `PQ` reaches
  `--phaseQualityThresh`, default 20. The **VCF v4.3 specification**
  (`docs`) defines `PQ` as the phred-scaled probability that alleles are
  ordered incorrectly within the phase set, and says the measure itself is
  not yet specified.
- **PhaseQuality** (Damaraju et al. 2026, section F, `paper`) is a random
  forest that puts variant pairs into three tiers from sequencing features.
  The tiers are not probabilities.
- **LongPhase 2** (section H, `abstract`; README re-read 2 Oct; an unpublished manuscript
  repository) scores phased variants with a graph neural network and
  unphases those it predicts misphased. Its README records a deployed fixed
  threshold of 0.30 and a calibration supplementary figure marked as a
  placeholder.
- **PhaseME** (Majidian and Sedlazeck 2020, `abstract`) checks any phaser's
  VCF against population linkage and splits blocks at conflicts.

### 2.2 Phase confidence in statistical phasing

- **SHAPEIT5** (Hofmeister et al. 2023, `paper`) gives each rare
  heterozygous genotype "the probability of the reported phase". It reports
  switch error by threshold with the share retained, and says a 0.99
  threshold "controls the SER to a maximum of around 2%" on WGS and 1 per
  cent on WES. That is a measurement against truth, not a bound.
- **Lassen et al. 2024** (`paper`) validates SHAPEIT5's posterior against 96
  trios and against WhatsHap read-backed phasing of 453,514 variant pairs,
  calls it well calibrated, and chooses 0.9 as an empirical cutoff for
  compound heterozygous calls.

### 2.3 Conformal prediction and risk control

- **Split conformal** (Lei et al. 2018, `abstract`) gives finite-sample
  marginal coverage for any model under exchangeability. Vovk et al. 2005 and
  Papadopoulos et al. 2002 are the origins; both are `recalled`.
- **Conformal risk control** (Angelopoulos et al., ICLR 2024, `paper`)
  bounds `E[L_{n+1}(λ̂)] ≤ α` when each loss is non-increasing in `λ`,
  bounded, and at most `α` at `λ_max`. Proposition 1 shows the bound can
  fail by an arbitrary amount for non-monotone losses.
- **Learn then Test** (Angelopoulos et al. 2025, `paper`) controls any risk
  in high probability by multiple testing. Its selective classification
  example controls "the error conditionally on predicting" with binomial tail
  p-values on i.i.d. calibration data. **Risk-controlling prediction sets**
  (Bates et al. 2021, `abstract`) are the earlier high-probability form.
- **Conformal selection** (Jin and Candès 2023, `abstract`) selects units
  while controlling the false discovery rate among those selected.
- **Dependent data:** block permutations for serial dependence (Chernozhukov
  et al. 2018, `abstract`); a coverage penalty for split conformal under
  mixing (Oliveira et al. 2024, `abstract`); weighted quantiles under drift
  (Barber et al. 2023, `abstract`); local spatial exchangeability (Mao et al.
  2024, `abstract`); groups drawn from distinct distributions (Dunn et al.,
  `abstract`).

### 2.4 Selective classification

- **El-Yaniv and Wiener 2010** (`abstract`) name the risk-coverage trade-off
  and characterize it for noise-free models.
- **Geifman and El-Yaniv 2017** (`paper`) define selective risk as
  `E[ℓ g] / E[g]` and give SGR, which picks a confidence threshold so that
  `Pr{R > r*} < δ` over an i.i.d. sample while maximizing coverage.
- **Hallberg Szabadváry et al. 2025** (`abstract`) obtain a reject option
  with a proved error rate by accepting only singleton conformal sets.

### 2.5 Genomics uses of conformal methods

- **SAVANA** (Elrick et al. 2025, `paper` for the passages read) offers
  Mondrian conformal prediction for somatic SV calls from long reads, under
  an i.i.d. assumption.
- **Papangelou et al. 2025** (`abstract`) review conformal prediction in
  genomic medicine.
- We did not find a conformal or risk-controlled method for phasing,
  haplotype assembly or germline small-variant calling.

### 2.6 Benchmarking guidance and clinical use

- **Krusche et al. 2019** (`paper`, phasing passages), from the GA4GH
  Benchmarking Team with GIAB, say phasing assessment methods "have not been
  fully developed or standardized" and that phasing must be judged on
  "completeness of phasing coverage" as well as accuracy. They propose no
  metric that joins the two.
- **Richards et al. 2015** (`paper`, in-trans passages) make in-trans
  detection moderate evidence (PM3) and describe parental testing as the
  way to establish phase.
- **Del Gobbo and Boycott 2025** (`paper`, read-based phasing section) say
  laboratories "primarily rely on familial genotyping or direct read-based
  phasing if possible", and describe long-read phasing by block N50 and genes
  covered.
- The ClinGen PM3 recommendation was not obtained (`recalled`).

## 3. Section 2 of the plan against prior work

Each part of the question, with the prior work that already covers it.

1. **A trust flag per pair with "error among trusted at most `α`",
   keeping as many pairs as possible.** As an abstract problem this is
   selective classification with guaranteed risk. Geifman and El-Yaniv 2017
   state it and solve it for i.i.d. data, and Learn then Test gives the same
   guarantee by multiple testing. Conformal selection (Jin and Candès 2023)
   controls the matching false discovery rate. The plan's claim that it
   "states a new formal problem" therefore does not stand. What is new, as
   far as these searches reach, is the phasing instance: a structured,
   dependent output, scored against a truth with its own errors.
2. **Risk-coverage benchmarking of phasers.** Partly covered. HapCUT2's
   paper compares its pruning and splitting against ProbHap's scores as
   error against the share pruned, and compares tools on the variants all of
   them phased. SHAPEIT5 reports switch error by confidence threshold with the
   share retained. Krusche et al. 2019 state that phasing accuracy and
   completeness must be read together. The risk-coverage curve is El-Yaniv
   and Wiener's.
3. **A tool-agnostic per-boundary score.** Partly covered. PhaseME is a
   tool-agnostic post-processor, but it uses population linkage and not
   reads. PhaseQuality scores variant pairs from sequencing data, in tiers.
   HapCUT2, ProbHap and GATK compute model posteriors inside one phaser.
4. **A conformal step with a finite-sample bound.** The methods exist
   (conformal risk control, Learn then Test, SGR, conformal selection). Their
   use in long-read genomics is limited, as far as we found, to SAVANA's
   Mondrian conformal prediction for SV calls.
5. **The guarantee under dependence, with windows or chromosomes as
   units.** General tools exist: block permutations, mixing penalties,
   spatial and hierarchical conformal prediction. Dunn et al.'s two-layer
   model matches pairs grouped in windows.
6. **The least abstention a target `α` needs, from the vote model.** No prior
   work found addresses it.
7. **Clinical use.** In-trans calls feed PM3, and laboratories use
   read-based phase when they can. Biobank compound-heterozygosity studies
   threshold statistical phase confidence at 0.9 or 0.99, chosen
   empirically.

## 4. What remains open

These are reports of the searches in section 1, not properties of the
literature.

- We did not find a published finite-sample error guarantee, in expectation
  or in high probability, for read-based phasing or haplotype assembly. The
  closest results are fixed posterior thresholds (HapCUT2, ProbHap, GATK),
  tiers (PhaseQuality), an unpublished classifier with a fixed threshold
  (LongPhase 2), and empirical error at a threshold for statistical phasing
  (SHAPEIT5, Lassen et al.).
- We did not find a calibration measurement of HapCUT2's switch quality, of
  GATK's `PQ`, or of ProbHap's transition score against a phase truth.
- We did not find a risk-coverage comparison of current long-read phasers
  (WhatsHap, HiPhase, LongPhase, HapCUT2) under one tool-agnostic score, or
  a published proposal for comparing phasers at matched completeness.
- We did not find conformal, risk-control or selective-classification
  methods applied to phasing, haplotype assembly or germline small-variant
  calling.
- We did not find a prediction of the risk-coverage curve of phasing from a
  read-vote model.
- We did not find clinical guidance that sets an error rate for a read-based
  in-trans call. The ClinGen PM3 document, which may bear on this, was not
  read.

## 5. Changes the plan needs

1. **Drop "states a new formal problem".** Cite Geifman and El-Yaniv 2017,
   Learn then Test and Jin and Candès 2023 for the problem, and claim the
   phasing instance, the dependence and the evidence score.
2. **Match the procedure to the loss.** Section 3.2 names conformal risk
   control for "switch errors among trusted pairs, over trusted pairs". Per
   window, that ratio is not monotone in the threshold, so Theorem 1 of
   Angelopoulos et al. does not cover it, and their Proposition 1 shows the
   miss can be large. There are three options:
   - Learn then Test or SGR, for a high-probability bound on the selective
     risk;
   - conformal selection, for the false discovery rate among trusted pairs;
   - a monotone loss such as trusted errors over all compared pairs in the
     window, which bounds a different quantity.

   The design must name the quantity bounded: an expected per-window ratio,
   a pooled ratio, or an FDR. These differ when windows differ in size.
3. **State the unit honestly.** A finite-sample guarantee "on a new genome"
   needs genomes as units. Section 7 calibrates on three genome-platform
   pairs, which makes the slack `B/(n+1)` of conformal risk control at least
   a quarter of `B`. Make the guarantee over windows, and keep the
   cross-genome claim as the empirical test that section 7 already
   registers.
4. **Add the prior scores as baselines:**
   - HapCUT2 switch quality (`--error_analysis_mode 1`);
   - ProbHap's transition score, where it runs;
   - PhaseQuality tiers;
   - LongPhase's `--edgeThreshold` behaviour and WhatsHap's own unphased
     sites;
   - LongPhase 2, if it is released before registration.
5. **Cite the benchmark precedents.** Cite HapCUT2's pruning trade-off and
   intersection comparison, SHAPEIT5's threshold table and Krusche et al.
   2019 where the benchmark is introduced. Its claim is then a multi-tool,
   tool-agnostic comparison at matched completeness against GIAB v5.0q.
6. **Read before registering:** the HapCUT2 supplement (the trade-off curves
   and the posterior derivation), the ClinGen PM3 recommendation, and the
   full texts of Jin and Candès 2023 and Dunn et al., whose exchangeability
   conditions decide which unit can be used.

## 6. Recommendation

**The gate passes with changes.**

The pass condition in section 6 of the plan holds: these searches found no
published method that gives a finite-sample error guarantee for read-based
phasing. The plan stays a method project and does not become a benchmark and
extension of an existing guarantee.

The changes in section 5 are needed because two of the plan's claims are
covered by prior work. The formal problem is selective classification with
guaranteed risk, published for i.i.d. data in 2017 and 2025. Risk-coverage
comparison of phase confidence scores appears in the HapCUT2 paper. The
procedure named in section 3.2 also does not, as stated, cover the loss the
plan defines.

The remaining contribution has four parts:
- the guarantee for dependent, window-structured phasing output;
- a read-based, tool-agnostic evidence score;
- the vote-model prediction of the risk-coverage curve;
- a multi-tool benchmark of long-read phasers at matched
  completeness.

Each is supported only as "we did not find". No deliverable may depend on
those absences, so the paper should state each one as a result it shows,
not as a gap it fills.
