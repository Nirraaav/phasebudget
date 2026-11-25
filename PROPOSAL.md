# The PhaseBudget proposal

## Exact and adaptive optimization for coverage-constrained haplotype assembly

**Course:** DH607, Autumn 2026 / Spring 2027, Advanced Research Project
**Working package name:** `phasebudget`
**Duration:** 24 weeks
**Author:** Nirav Bhattad (23B3307), working alone
**Deliverables:** a Python library with exact solvers, a reproducible benchmark
suite, accuracy-completeness-runtime phase diagrams, and a research report

**Outcome, 25 November 2026:** Part A of `docs/submission/paper.tex` answers
the selection question; the project then moved to certified phase relations
(`README.md`).

## Hypothesis

Haplotype assembly trades computational tractability against retention of
informative sequencing reads. This project characterizes coverage-constrained
read selection using exact optimization, quantifies the gap between deployed
heuristics and optimal surrogate objectives, and develops a solver-guided
adaptive selection method whose effectiveness is studied jointly with the
fidelity of the weighted MEC objective.

### The hypothesis as the measurements left it, 26 September 2026

The form submitted to DH607 states the hypothesis as "choosing the kept reads
optimally, rather than by rule of thumb, makes phasing more accurate". On
simulation that holds at low caps and at high error, and not at the deployed
cap: at cap 15 the selector WhatsHap ships reaches 0.922 to 0.968 of the
additive optimum on the 60-variant grid and 0.915 to 0.916 at chromosome
scale, and its switch count, summed over 25 replicates, is within one of the
optimum's at error rates up to 0.10. (Until 26 September these figures came
from the paper's version of the selector, which WhatsHap's code does not
follow; OQ-14.) The testable form is therefore conditional: **the accuracy
value of optimal read selection is concentrated where the cap is small
relative to depth or the reads are noisy, and vanishes at the deployed cap for
accurate reads.** The real-data run decides whether real reads, with
correlated errors and uneven depth, fall inside that region. Part A of the paper has
the real-data tables.

### The real-data run, 26 September 2026

The registered HG002 run (`docs/design/real-data-plan.md`)
confirms the deployed gap (0.937 of the optimum at cap 15) and cannot
decide the accuracy hypotheses: 35-fold HiFi reads let every sensible
pipeline phase chromosome 20 with no real switch error from cap 5 up, and
the few differences at cap 3 are too small for the registered bootstrap.
The conditional hypothesis is therefore untested on real data rather than
refuted; the reads where it could be tested are noisier ones. The adaptive
loop gave no gain and at cap 3 exposed a defect in stitched acceptance.

### Component 3 as the measurements left it, 26 September 2026

The adaptive loop lowers switch error, by a factor of about 3 at chromosome
scale once rounds are accepted window by window. A control that chooses the
fragile boundaries at random keeps almost all of that, so the proposal's
mechanism, MEC margins steering reselection, is not what produces the gain;
reselecting under perturbed weights and rescoring against every read is. The
central hypothesis's first half, that phase-confidence signals decay as the
objective stops tracking the truth, holds on simulation. Its second half, that
the adaptive advantage follows the signal, does not arise: the advantage does
not depend on the signal. OQ-13 (`docs/literature/open-questions.md`) has the numbers.

## 1. Problem

Read-based haplotype assembly reconstructs the two chromosome copies of a
diploid individual from sequencing reads spanning heterozygous variants. The
data reduce to a fragment matrix `R in {0,1,-}^(n x m)` with confidence weights
on observed alleles, and the standard objective is weighted minimum error
correction: find the haplotype pair and read assignment minimizing the total
weight of corrected alleles.

Exact wMEC is available through a dynamic program whose state is a bipartition
of the reads active at a variant, so its cost is exponential in active coverage
and linear in the number of variants. Sequencing depth is set by the
experiment, not by the assembler, so reads are discarded before solving. Which
reads to keep is an optimization problem in its own right, and it is solved in
deployed tools by a heuristic.

Three questions follow, and they are the three components of the project.

## 2. Component 1: when is wMEC a faithful surrogate?

The objective is not the accuracy. wMEC counts corrected alleles; a user cares
about switch errors against the truth. The two coincide when the data are
informative enough, and the first component measures where they stop.

The instrument is the **objective-truth gap**: the wMEC cost of the true
haplotype minus the optimal cost. It is zero when the truth is an optimal
solution and positive when the objective prefers something else. Simulated
instances with known truth vary coverage, read span, heterozygosity, error
rate, and base quality independently.

**Evidence:** an exact 0-1 ILP for wMEC using binary read-assignment and
haplotype variables with standard XOR linearization and a symmetry-breaking
constraint, cross-checked against the coverage-parameterized dynamic program
and against exhaustive enumeration on small instances. All three are
implemented and agree; `tests/test_ilp.py` and `tests/test_dp.py` carry the
comparisons.

**Expected outcome:** a phase diagram of the parameter range in which the wMEC
optimum coincides with the truth, and the boundary at which it stops doing so.

## 3. Component 2: exact coverage-constrained read selection

Each read induces an interval `[l_r, u_r]` over ordered heterozygous variants.
For a maximum permitted active coverage `c`,

```text
sum over r with j in [l_r, u_r] of x_r <= c
```

must hold at every variant `j`. For an additive score `w_r` the problem is to
maximize `sum_r w_r x_r` subject to that. Because the interval-incidence matrix
has the consecutive-ones property, the linear relaxation is integral, and the
problem has an exact min-cost-flow formulation: a path of `c`-capacity arcs
along the variants, a unit arc of cost `-w_r` per read, and `c` units pushed
end to end. The argument is written out in
`docs/design/selection-integrality.md`.

This is not claimed as new. Maximum-weight interval selection on `k` machines
is a textbook flow problem (Arkin and Silverberg 1987), and Mäkinen, Staneva,
Tomescu and Valenzuela 2015 already bring interval-scheduling machinery to read
selection, for a maximum-minimum-coverage objective that this project
implements alongside the additive one. The observation that phasing read
selection is an instance of the additive problem is an observation. The
contribution is the measurement: how far the deployed heuristic sits from the
optimum of its own objective, and whether closing that gap improves phasing
accuracy. Measured on simulation as of 24 September 2026, that gap is 7 to 13
per cent.

**Evidence:** the outcome difference between published read selectors and the
additive exact selector,
across coverage, read span, and cap, with worst cases reported alongside means.

**Theoretical extension, now partly settled:** the tractability depends on
counting a read's whole span against the cap. Counting only the variants it
calls makes the linear relaxation fractional and the selection problem NP-hard,
by reduction from independent set, and already APX-hard with one hole per read
by Bar-Yehuda et al. 2006. `docs/design/coverage-convention.md` gives both and
the tests check both embeddings. This part is settled and is not a
contribution. A separate question is whether the tractable class extends to
redundancy-sensitive or non-additive objectives, for instance one that pays
once for covering a boundary rather than once per read. This is fenced off as
OQ-3 and nothing in the deliverable depends on it.

## 4. Component 3: adaptive solver-guided read selection

After solving wMEC on a coverage-limited subset, compute the fragility of every
adjacent phase decision using its **MEC margin**: the difference between the
optimum and the best solution with that local phase relation reversed. Rather
than re-solving once per boundary, the exact margins are extracted block-wise
from the forward and backward tables of the dynamic program, and verified
against brute-force constrained re-solving on small instances.
`docs/design/margin-extraction.md` derives the extraction and
`tests/test_dp.py` runs the verification.

Reads spanning low-margin boundaries then receive increased weight, the exact
selector chooses a revised subset, and wMEC is solved again. The outer loop
uses a fixed iteration cap and retains the best solution encountered, scored
against the full read set so that subsets are comparable, which prevents
oscillatory reweighting from degrading the answer.

**The central hypothesis** ties this component to the first: the adaptive
signal is most useful precisely where wMEC itself tracks biological truth. As
the objective-truth gap widens in high-error or low-information regimes,
MEC-derived margins should become less informative, and the adaptive advantage
should shrink.

## 5. Evaluation

Compared methods: random downsampling, the greedy heuristic, exact static
selection, adaptive selection, and at least one external phasing tool.
Reported quantities: switch and flip error, phased fraction, block N50 and
QAN50, wMEC score under the full read set, selection-objective gap, runtime,
and peak memory. Metric definitions are in `docs/design/metrics.md`.

Simulation provides ground truth and parameter control. The principal
conclusions are then validated on **GIAB HG002 chromosome 20**, which is the
conventional target and makes the numbers comparable to published ones.

## 6. How the project can be shown wrong

**Solver correctness.** The three solvers must agree wherever two of them can
run. One disagreement invalidates every measurement made with the faster one.

**Margin correctness.** Block-wise margins must equal constrained re-solving on
every small instance. The check is exhaustive over a bounded domain, not a spot
check.

**Selection exactness.** The flow selector must match subset enumeration on
small instances, and its output must be feasible at every variant.

**Accuracy claim.** Adaptive selection must reduce switch error against exact
static selection at the same cap on preregistered simulated workloads, and the
reduction must survive on real data. A null result is reportable and would say
that margins do not carry usable information at realistic error rates.

**Surrogate claim.** If closing the selection-objective gap does not improve
accuracy, the honest conclusion is a negative result about the objective. The
first benchmark run already points that way, and the report will say so rather
than quietly changing the objective.

**Reproducibility.** Every benchmark records its commit, command line, machine,
and seeds. Simulated results are labelled simulated wherever they appear.

## 7. Six-month timeline

| Weeks | Work | Milestone |
|---|---|---|
| 1-3 | prior-art matrix, formal specification, exact oracles, three-solver agreement | M0: novelty gate passed or the thesis narrowed |
| 4-6 | objective-truth gap over the full parameter grid | M1: phase diagram for component 1 |
| 7-9 | exact selection at scale, heuristic gap, worst cases | M2: component 2 measured |
| 10-13 | margin extraction verified, margin-versus-switch-error analysis | M3: the central hypothesis tested on simulation |
| 14-17 | adaptive loop, sensitivity to boost and iteration cap | M4: component 3 measured, or a negative result recorded |
| 18-21 | GIAB HG002 chromosome 20, external tool comparison | M5: real-data validation |
| 22-24 | phase diagrams, report, artifact audit | M6: reproducible release candidate |

## 8. Risks and mitigations

| Risk | Consequence | Mitigation |
|---|---|---|
| The adaptive method is already published | the novelty claim fails | milestone 0 precedes the engineering; the measurement contribution stands alone |
| Margins turn out to be uninformative | component 3 has no signal | the decay curve is itself the result; OQ-5 states what would settle it |
| The surrogate objective and accuracy diverge | optimizing selection exactly does not help | already anticipated; the two are reported in separate columns everywhere |
| Real data does not fit the state limit after selection | no real-data validation | OQ-9 fixes the cap from the read length distribution before the run |
| The complementary-haplotype model misfits real genotypes | accuracy loss attributed to the wrong cause | the general bipartition model is implemented as a small-scale comparison |
| Real data acquisition consumes the project | no result | the simulation track is complete and self-contained; real data is validation, not the deliverable |

## 9. Division of work

The project was done by one person, Nirav Bhattad, who owned every part: the
solvers and their oracles, selection, the adaptive loop, the metrics, the
real-data pipeline, the external tool comparison and the evaluation.

## 10. Initial references

This list was written before the review. The status of every entry, and the
entries added since, are in `docs/literature/01-prior-art-matrix.md`; cite from
there, not from here.

- Lancia, G. et al. *SNPs Problems, Complexity, and Algorithms*. ESA, 2001.
- Cilibrasi, R., van Iersel, L., Kelk, S., Tromp, J. *On the Complexity of
  Several Haplotyping Problems*. 2005/2007.
- Bonizzoni, P. et al. *On the Minimum Error Correction Problem for Haplotype
  Assembly in Diploid and Polyploid Genomes*. 2016.
- Patterson, M. et al. *WhatsHap: Weighted Haplotype Assembly for
  Future-Generation Sequencing Reads*. J. Comput. Biol., 2015.
- Fischer, S. M., Marschall, T. *Selecting Reads for Haplotype Assembly*. 2016.
- Edge, P., Bafna, V., Bansal, V. *HapCUT2*. Genome Research, 2017.
- Pirola, Y. et al. *HapCol*. Bioinformatics, 2016.
- Mäkinen, V., Staneva, V., Tomescu, A. I., Valenzuela, D. *Interval Scheduling
  Maximizing Minimum Coverage*. 2015. Interval-scheduling machinery applied to
  read selection, for a different objective.
- Arkin, E. M., Silverberg, E. B. *Scheduling Jobs with Fixed Start and End
  Times*. Discrete Applied Mathematics, 1987.
- Hoffman, A. J., Kruskal, J. B. *Integral Boundary Points of Convex
  Polyhedra*. 1956. Total unimodularity of interval matrices.
- Zook, J. M. et al. Genome in a Bottle benchmark sets. HG002 chromosome 20 is
  the evaluation target.
