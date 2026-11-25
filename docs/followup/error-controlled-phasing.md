# Follow-up project: error-controlled haplotype phasing

Written on 1 October 2026, after every registered test of the PhaseBudget
paper was scored, except the ONT tool-level test. It is a plan. It claims
nothing, and its literature statements are `recalled` until the literature
gate (section 6) has run.

## Amendments

1. **2 October 2026, after the literature gate**
   (`docs/followup/literature-gate.md`). The gate passes with changes. We
   found no published finite-sample error guarantee for read-based
   phasing. Three things change:
   - *The formal problem is not new.* Keeping the error among trusted
     pairs at or below `alpha` is selective classification with guaranteed
     risk (Geifman and El-Yaniv 2017). That work assumes independent,
     identically distributed units. The contribution is its use on a
     dependent genomic output: a tool-agnostic evidence score, a guarantee
     over windows, the benchmark at matched completeness, and the
     vote-model risk-coverage theory.
   - *The calibration method changes.* Conformal risk control (Angelopoulos
     et al. 2024) needs a loss monotone in the threshold. The error among
     trusted pairs is a ratio and is not monotone. Section 3.2 now uses
     Learn then Test (Angelopoulos et al. 2025), with fixed-sequence testing
     over a grid of thresholds and 1 Mb windows as units.
   - *The guarantee is over windows, not genomes.* Three calibration
     genomes are too few for a per-genome guarantee. Each test genome is a
     separate check that the window-level guarantee carries over to a new
     individual or platform, and every check is reported.

   Development gates G-A to G-C were run before this amendment (development notebook,
   2 October). The evidence score is not frozen yet; it will be fixed in
   the registration.

2. **4 October 2026, after the registered test was scored.** E1 failed on
   HG001, E2 was partly supported and E5 supported (development notebook, 4 October). The
   next stage is the extension plan (development repository): per-stratum guarantees,
   more test genomes, a check of the HG001 truth, and relations between any
   two variants in a gene.

## 1. What the PhaseBudget results say about where to go

Read selection is not the lever. At WhatsHap's default cap, selection makes
no detectable difference, inside the tool or in the library. With tied sites
left unphased, most of the effect measured at tight caps disappears. What
selection still changes is completeness: how many sites the reads can
decide at all.

The errors that remain at default settings are few, concentrated, and
mostly invisible to the user:
- **They are rare in easy sequence.** On HiFi at full depth, `whatshap
  phase` has 2 switch errors in 1,438,833 easy-region pairs. In tandem
  repeats and homopolymers it has 9 in 17,127, which is 5.25 per 10,000
  (`benchmarks/results/runs/2026-10-01_strata_hifi_full_tools.txt`).
- **They grow at low depth and on ONT.** At 10-fold ONT, `whatshap phase`
  has 1.24 per 10,000 overall and 12.22 in tandem repeats. LongPhase has
  1.63 overall and 1.89 in tandem repeats
  (`2026-10-01_strata_ont_10x_tools.txt`).
- **Tools differ in what they leave unphased.** On ONT at full depth,
  LongPhase has 0 tandem-repeat errors against 13 for `whatshap phase`, and
  it compares 970 fewer tandem-repeat pairs. WhatsHap leaves tied sites
  unphased; the library guessed at them. On ONT chromosome 20 at 10-fold,
  21 tied sites held 17 of the library's 31 errors
  (`2026-10-01_ties_ont_chr20_10x.txt`).
- **Wrong calls can be ranked, but not calibrated.** The library's
  per-boundary posterior margin separated wrong calls from right ones with
  AUC 0.997 on HiFi (G9). Its probabilities were off by a factor of 2 to 5,
  and a frozen correction failed its registered test on ONT (C).
- **Read quality varies, and the alleles do not show it.** A read's gap
  fraction predicts its allele error 7 to 40 times over
  (`2026-10-01_gap_error_ont_chr20.txt`), while every realigned allele
  carries phred 30.

So today a phased VCF says which variants are in a block and nothing about
how far to trust any one relation inside it. Each tool also decides on its
own, by unstated rules, which sites to leave out. Clinical uses need
exactly this information. Calling compound heterozygosity, for example,
asks whether two variants in one gene are in trans.

## 2. The question

Given reads and any phaser's output, can we report, for every pair of
adjacent phased variants, whether its relation is trusted? We want a
guarantee that the switch error among trusted pairs is at most a target
`alpha`, and we want to keep as many pairs trusted as possible.

This is selective prediction with a risk guarantee, applied to haplotype
phasing. The project would:
1. **Benchmark** existing phasers on risk-coverage curves, so that tools
   are compared at matched completeness instead of at whatever completeness
   each one picks.
2. **Build a method:** a tool-agnostic per-boundary evidence score, computed
   from the reads and the phaser's output, and a conformal risk-control step
   that turns the score into a threshold with a finite-sample bound on the
   error among trusted pairs.
3. **Give theory:** the guarantee under dependence along the genome, with
   windows or chromosomes as the exchangeable units; and, from the vote
   model with abstention, the least abstention a target `alpha` needs given
   a depth and an allele error. The PhaseBudget vote model becomes a
   predicted risk-coverage curve.

## 3. Components

### 3.1 Evidence score per boundary

For each adjacent pair `(j, j+1)` that the phaser puts in one block,
computed from the BAM and the phased VCF. No truth is read. The candidate
features are:
- **Support margin:** for each read calling both sites, whether it agrees
  with the called relation, weighted. This is the tie count of
  `benchmarks/tie_rescore.py`, made per boundary.
- **Indirect support:** reads with calls on both sides of the boundary that
  do not call both sites (the `c` count of `benchmarks/explore_levels.py`).
- **Read reliability:** the gap fraction of each supporting read, and its
  allele disagreement with the phaser's haplotypes elsewhere.
- **Sequence context:** whether either site is in a tandem repeat, a
  homopolymer, a segmental duplication or a low-mappability region (GIAB
  v3.6 strata), and the distance to the nearest indel.
- **The library posterior margin** where the library can compute it (G9).

The first score is the support margin alone. A learned combination comes
second, and only if a development gate shows the features add to it.

### 3.2 Risk control

Split calibration by chromosome or by 1 Mb window. Choose the threshold `t`
by Learn then Test (amendment 1): test, for each threshold on a fixed grid
from strict to loose, the null that the error among trusted pairs exceeds
`alpha`, using window-level statistics, and stop at the first acceptance. The units must be
exchangeable, which is tested and not assumed:
- within a genome, across windows;
- across genomes of the same platform and depth.

Where exchangeability across genomes fails, use calibration conditional on
platform and depth, and report it.

### 3.3 Theory

- State the guarantee for windows as units, with the loss "switch errors
  among trusted pairs, over trusted pairs", following the conformal risk
  control literature (to be verified, section 6).
- Derive, from the vote model with abstention at ties and below a margin,
  the risk-coverage curve as a function of the informative-count
  distribution and `q`. Test it the way the ONT prediction was tested, with
  a better regime rule (section 7).

## 4. Data

| data | truth | status in this repository |
|---|---|---|
| HG002 HiFi (Sequel II 11 kb) | GIAB v5.0q | run, fragments kept, BAM streamed |
| HG002 ONT R10.4.1 sup 60x | v5.0q | run, BAM streamed from ONT open data |
| HG005 HiFi | GIAB transmission | run |
| HG002 ONT R10.4.1 hac 60x | v5.0q | available (`ont-open-data`, `giab_2023.05/.../hg002_hac_60x`, 215 GB), unused |
| HG001 ONT R10.4.1 sup 60x | Platinum Genomes NA12878, pedigree-phased: 49,908 of 49,908 chromosome-20 het SNVs phased | available (197 GB), unused |
| HG001 HiFi | Platinum | to locate |
| HG003, HG004 ONT sup | to check: GIAB parental truths may be unphased | available |

HG001 and HG002 hac are new data. Keep them for confirmatory tests only.

## 5. Development gates, in order

Each gate uses only chromosome 20 of data already seen, or documents. A
gate that fails stops the project or changes it, and the failure goes into
The development notebook.

- **G-A, the score separates errors from correct calls.** On chromosome 20,
  for `whatshap phase` on HG002 HiFi full depth and ONT full and 10-fold:
  - compute the support margin for every compared pair;
  - report the AUC for wrong pairs and the risk-coverage curve.
  - *Go* if trusting the top 99 per cent of pairs cuts errors at least
    5-fold on ONT at 10-fold, and at least 3-fold on HiFi.
  - *No-go* if errors are spread evenly over the score, which would mean the
    remaining errors are systematic: truth errors, paralogs, or allele
    miscalls that every read shares.
- **G-B, the tandem-repeat errors are caught.** The share of
  tandem-repeat and homopolymer errors among the lowest-scoring 1 per cent
  of pairs. If it is low, add sequence-context features and re-run G-A once.
- **G-C, the truth is not the error.** Inspect the 32 HiFi `whatshap
  phase` errors genome-wide:
  - read support for the called relation against the truth's;
  - agreement with the HG002 transmission truth;
  - agreement with HiPhase and LongPhase.

  If most "errors" are supported by every tool and every read, they are
  truth problems, and HiFi cannot calibrate anything. Then the project is
  ONT and low depth only.
- **G-D, enough errors to calibrate.** At `alpha` near 1 per 100,000, a
  HiFi genome has about 30 errors. Simulate the conformal bound's looseness
  at that count, by resampling windows on HG002. If it is too loose, the
  targets are ONT and low depth, where errors are 10 to 100 times more
  common.
- **G-E, literature** (section 6).

## 6. Literature gate (must run before any design is registered)

Every entry is `recalled` and may not support a claim until it is read:
- **HapCUT2:** its phase-quality score and pruning of low-confidence
  variants (Edge, Bafna and Bansal 2017). What it outputs, and whether
  anyone has measured its calibration.
- **WhatsHap, HiPhase and LongPhase:** what each writes as a per-variant or
  per-block quality, and their rules for leaving a site unphased.
- **GATK ReadBackedPhasing and its `PQ` field.**
- **Conformal prediction and risk control:** split conformal; conformal
  risk control (Angelopoulos and coauthors); Learn then Test; and conformal
  methods for dependent or spatial data.
- **Selective classification and risk-coverage curves:** El-Yaniv and
  Wiener; Geifman and El-Yaniv.
- **Any genomics use of conformal methods,** especially for variant calling,
  phasing or haplotype assembly.
- **GIAB guidance on phasing benchmarks:** how switch error and phased
  fraction are reported together, and any matched-completeness proposal.
- **Clinical phasing:** how compound-heterozygosity calls use read-based
  phase today.

Rows go into `docs/literature/01-prior-art-matrix.md`, section I, with
their status. If a published method already gives a finite-sample error
guarantee for read-based phasing, the project changes to a benchmark and
extension of that method, and says so.

## 7. Registration, learned from PhaseBudget

- Register each confirmatory test before its data are phased, with GitHub
  push times as the third-party record.
- Score every threshold rule at full precision from committed per-unit
  sums, never from printed text (the R2 rounding error).
- Use decision rules that test the claim, not a fixed line in effect
  units. A regime call should ask whether an interval excludes 0, not
  whether an estimate passes 0.5 (both R2 failures).
- Leave tied sites unphased in every pipeline, as the tools do. Report
  completeness beside accuracy in every table.
- Keep development to chromosome 20 and confirmation to new genomes.

**Confirmatory design, to be written after the gates:**
- **Calibrate:** on HG002 (HiFi, ONT sup) and HG005.
- **Test the guarantee:**
  - on HG001 ONT sup with the Platinum truth;
  - on HG002 ONT hac, the same genome with a new error profile;
  - at depths 10, 20 and full.
- **Primary claim:** realized error among trusted pairs at most `alpha` in
  at least the stated share of test genome, platform and depth cells, with
  phased-pair retention reported.
- **Secondary claim:** at matched completeness, the ranking of existing
  phasers changes relative to their default-output ranking.

## 8. Schedule

Machine estimates come from the PhaseBudget runs: an ONT genome streamed
per chromosome takes about 1 to 2 days, and disk is the limit.

| weeks | work |
|---|---|
| 1 | gates G-A to G-E on chromosome 20 and on documents |
| 2 | evidence-score tool from BAM and VCF, unit tests, a differential check against `tie_rescore.py` |
| 3 | risk-coverage benchmark of the four phasers on the seen data |
| 4 | theory: the guarantee and the vote-model risk-coverage prediction, with simulation checks |
| 5 | registration; fetch HG001 ONT and HG002 hac per chromosome |
| 6 to 7 | confirmatory runs and scoring |
| 8 to 10 | paper, supplement, release |

## 9. Risks

- **The remaining errors are systematic** (G-A, G-C). Abstention cannot
  catch errors that every read supports. Then the project narrows to ONT and
  low depth, or changes to sequence-context allele calling.
- **HiFi errors are too rare to calibrate** (G-D). Then the claims are
  limited to ONT and low depth.
- **Prior art already gives a guarantee** (section 6). Then the project
  becomes a benchmark and extension.
- **Exchangeability fails across genomes.** Then the guarantee is stated
  per platform and depth, and the failure is reported as a result.
- **Disk and network:** 197 to 215 GB per new genome, streamed one
  chromosome at a time with the existing disk guards.

## 10. Alternatives kept in reserve

- **Sequence-context allele calling for phasing in tandem repeats and
  homopolymers.** Most remaining HiFi errors sit there. Fixing allele calls
  at those sites is a method project, but its absolute gain is small (9 of
  32 HiFi errors).
- **Paralog-aware phasing in segmental duplications.** The clinical value
  is high, and specialized tools exist (Paraphase and others, recalled).
- **Read-based plus population-panel phasing for low-coverage long-read
  cohorts.** The upside is large and the field is crowded. Months of work,
  with compute this machine may not have.
