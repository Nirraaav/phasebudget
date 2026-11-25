# The novelty gate

## The rule

No deliverable in this repository may depend on a negative claim about the
literature. "Nobody has done X" is not a thesis, it is a bet against a search
that was never exhaustive.

The three components are therefore stated so that each produces a usable result
whichever way the review goes.

## Component 1: when is wMEC a faithful surrogate?

**Not a novelty claim.** This is a measurement. The objective-truth gap as a
function of coverage, span, error rate, and heterozygosity is a phase diagram,
and its value does not depend on whether someone has drawn it before. If a
published version exists, this one is a replication with the margin analysis
attached, and the report says so.

**Failure mode:** the gap turns out to be negligible everywhere in the
plausible parameter range, in which case component 3's central hypothesis has
nothing to attach to and the honest conclusion is that MEC-derived margins are
uniformly usable. That is a reportable result.

## Component 2: exact coverage-constrained selection

**Not a novelty claim, by construction.** The integrality argument is Hoffman
and Kruskal and Arkin and Silverberg, the flow network is Carlisle and Lloyd, and the fact that the phasing
instance is an interval instance is an observation, not a theorem. The
contribution is the measurement: how far the deployed heuristic sits from the
optimum of its own surrogate objective, and whether closing that gap changes
phasing accuracy.

**Failure mode, already partly observed:** closing the gap on the surrogate
does not reliably improve accuracy, because the surrogate is not the accuracy.
The first benchmark run shows adaptive selection beating exact selection on
switch error while scoring *lower* on the surrogate objective. If that holds
up, the result is a negative one about the objective, and it is more useful
than a win would have been.

The optional extension, classifying redundancy-sensitive objectives, is fenced
off in OQ-3 and nothing depends on it.

## Component 3: adaptive solver-guided selection

**This is where a novelty claim would live, and it is not yet made.** The parts
are all standard: constrained re-optimization for sensitivity, iterative
reweighting, exact selection. The question the review must answer is whether
the composition, applied to read selection for haplotype assembly, exists
already, in phasing or in a neighbouring problem. OQ-1 carries it.

**Exit criterion for milestone 0.** Either

1. name one precise question in components 1 to 3 that the literature does not
   answer, and state what evidence would settle it; or
2. narrow the project to the measurement contribution, which stands on its own,
   and say so in the report before any further engineering.

Rebranding an existing method is not an acceptable outcome of this gate.

## What has already been given up

**The phrase "WhatsHap's read selection" as a single object**, given up
2026-09-10 to Patterson et al. 2015. The paper prunes at random and lists
proper selection as an open direction in its discussion. The Fischer and
Marschall heuristic postdates it. Any comparison has to name which procedure
and which year it means.

**The claim that `phasebudget.dp` is a reimplementation of WhatsHap**, narrowed
2026-09-10 to the same paper. WhatsHap does not need the all-heterozygous
assumption and says so explicitly. The dynamic program here is the WhatsHap
recurrence specialized to complementary haplotypes, which is a restriction of
it, and the general model is available only as the brute-force oracle in
`exact.solve_bipartition_bruteforce`. OQ-7 carries the cost of the restriction
and it is no longer a small question.

**The framing of component 3 as adaptivity without qualification**, narrowed
2026-09-10 to Beretta et al. 2018. HapCHAT already holds the name "adaptive
haplotype assembly" for adapting a per-column error bound. Component 3 has to
be stated as adaptive *read selection* against a fragility signal, and the
report has to distinguish itself from HapCHAT explicitly rather than by
omission.

**The stated mechanism of the adaptive loop**, narrowed 2026-09-10 to this
repository's own measurement. `benchmarks/benchmark_margin_confound.py` shows
the subset margin has rank correlation 0.26 to 0.56 with the full-matrix
margin, so the loop does not measure the fragility of the phase decision. It
measures the fragility of the phase decision in the subset it just chose. The
accuracy gain survives this, because a full-matrix oracle does no better at a
matched round budget, but the mechanism sentence has to be the narrow one.

**Component 1 as a new question**, given up 2026-09-24 to Majidian, Kahaei and
de Ridder 2020. They ask when MEC prefers a wrong haplotype, across coverage,
error rate and platform, with a theorem and simulations. Component 1 is now an
extension: weighted MEC, exact optima, the gap as a magnitude, and its link to
margin informativeness.

**The margin as a new quantity**, given up 2026-09-24 to Kuleshov 2014.
ProbHap runs forward-backward over the same read-bipartition DP and reports a
per-adjacent-pair transition score, and says the MEC case of its DP gives
confidence associated with the MEC. The max-product margin here is that
construction in the other semiring. HapCUT2's switch quality is a third
per-boundary confidence. What component 3 can still claim is only the use of
such a confidence to reselect reads under a coverage cap.

**The gapped-coverage variant as owned theory**, given up 2026-09-24 to
Bar-Yehuda et al. 2006. One hole per read is already APX-hard through an
embedding of (2,2)-union graphs. OQ-11 is closed.

**The margin as a phasing-specific object**, given up 2026-09-25 to Kohli and
Torr 2008. A margin is a pairwise min-marginal energy minus the optimum, their
equations 5 and 6, computed exactly on a chain by forward and backward passes.
The report names it as such.

**Uncertainty-driven data selection as an idea**, given up 2026-09-25 to the
active-learning literature (margin sampling, as surveyed by Settles 2010) and
to Weilguny et al. 2023, whose BOSS-RUNS steers nanopore sequencing toward the
sites of highest genotype uncertainty. Component 3 is margin sampling applied
to boundaries, choosing among reads already sequenced under a coverage cap.

**Component 1 as the first real-data observation**, given up 2026-09-25 to
Duitama et al. 2012. On NA12878 fosmid data the trio truth has higher MEC than
HapCUT in every block, and HapCUT's MEC barely correlates with switch error.

**Margin-guided reselection as the mechanism**, given up 2026-09-26 to this
repository's own control. Moving the fragile boundaries to random places keeps
almost all of the adaptive gain at both scales (OQ-13). Component 3 is
restated as reselection under perturbed weights with rescoring against every
read, accepted window by window at chromosome scale; the margin is one
perturbation among several that work as well. That is closer to ensemble and
local-search methods than to margin sampling, and it is stated as a measured
effect, not as a new method.

## Where the gate stands on 25 September 2026

The review is complete as far as open access allows; section G of the prior-art
matrix lists the six rows that stay `recalled` and what stands in for each. No
component rests on novelty. Component 2 is a measurement nothing read reports.
Component 1 is an extension of Majidian et al. and Duitama et al., and it has
one result of its own that does not depend on the literature: the closed form
of Majidian et al. is an upper bound, shown by a tested counterexample and
measured in `benchmarks/benchmark_majidian_bound.py`. Component 3 is stated as
margin sampling applied to read reselection under a cap, and it stands or falls
on its measured accuracy change. Exit criterion 2 is met; criterion 1 is met by
the question section G states, and it is narrower than on 24 September.

## Where the gate stood on 24 September 2026

Exit criterion 2 is met and criterion 1 is met narrowly. The project stands on
the measurement of component 2, which no source found reports and which does
not depend on that absence. The one precise question the literature was not
found to answer is whether reselecting reads against a per-boundary phase
confidence improves accuracy under a fixed coverage cap; the evidence that
would settle it is `benchmarks/benchmark_sensitivity.py` and its real-data
successor. Section F of the prior-art matrix records the search.
