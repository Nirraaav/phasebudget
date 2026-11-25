# Milestone-0 prior-art matrix

## How to read the status column

| Status | Meaning |
|---|---|
| `recalled` | written down from background knowledge, not yet checked against the paper |
| `abstract` | the abstract has been read |
| `paper` | the full text has been read and the row reflects it |

Every row starts at `recalled`. A `recalled` row may not be cited in the report
and may not support a claim of novelty. Upgrading these rows is the first task
of milestone 0, and the deadline is in `open-questions.md` as OQ-1 and OQ-2.

## A. Complexity of the objective

| Work | Claim as recalled | Status |
|---|---|---|
| Lancia, Bafna, Istrail, Lippert, Schwartz 2001, [SNPs problems, complexity, and algorithms](https://doi.org/10.1007/3-540-44676-1_15), ESA, LNCS 2161, 182-193 | defines the SNP matrix, the fragment and SNP conflict graphs, and three problems: MFR (remove fewest fragments), MSR (remove fewest SNPs) and LHR. It does not define MEC; Cilibrasi et al. say MEC appears there only "implicitly". MFR is NP-hard with one gap per fragment and polynomial for gapless fragments, by network flow | `paper` |
| Lippert, Schwartz, Lancia, Istrail 2002, [Algorithmic strategies for the SNP haplotype assembly problem](https://doi.org/10.1093/bib/3.1.23), Brief. Bioinform. 3(1), 23-31 | defines MEC, the fewest entries to change so the matrix is feasible; Theorem 2: MSR and MEC reduce from MAXCUT, so both are NP-hard on general matrices | `paper`; this is the citation for MEC |
| Greenberg, Hart, Lancia 2004, [Opportunities for combinatorial optimization in computational biology](https://doi.org/10.1287/ijoc.1040.0073), INFORMS J. Comput. 16(3), 211-231; read as the author copy | section 3.2 names MEC "minimum letter flips" and introduces "its weighted version in which, for each entry, there is a weight for its flipping", proportional to the confidence of the base call | `paper`, section 3.2; this is the citation for wMEC |
| Cilibrasi, van Iersel, Kelk, Tromp 2007, Algorithmica; read as [arXiv q-bio/0508012v2](https://arxiv.org/abs/q-bio/0508012) | MEC is APX-hard when each fragment has at most one gap, and NP-hard in the gapless case; the paper questions earlier claims that MEC is NP-hard on binary (hole-free) matrices | `paper`, arXiv v2 |
| Bonizzoni, Dondi, Klau, Pirola, Pisanti, Zaccaria 2016, [On the MEC Problem for Haplotype Assembly in Diploid and Polyploid Genomes](https://doi.org/10.1089/cmb.2015.0220); read as the [HAL author copy](https://inria.hal.science/hal-01388448) | Theorem 5: MEC is not in APX under the Unique Games Conjecture. Theorem 7: `O(log nm)`-approximable, and FPT in the total number of corrections. Theorem 11: gapless MEC, without the all-heterozygous assumption, is FPT in fragment length, `O(3^l l n)`. Section 5: a 2-approximation for binary MEC. Theorems 13, 15, 17: `k`-ploid MEC is APX-hard, and FPT in (coverage, `k`) and in (fragment length, `k`) | `paper`, author copy |

If MEC is APX-hard, exact solving is only available at bounded coverage, which
is the premise of the whole project. Two corrections from the full texts. The
unconditional APX-hardness of diploid MEC is Cilibrasi et al., for fragments
with one gap; Bonizzoni et al. add inapproximability within any constant, and
that result is conditional on the Unique Games Conjecture. A sentence in the
report that says "MEC is inapproximable" has to carry that condition. The
parameterization this project relies on, (coverage) for the diploid case, is
Theorem 15 of Bonizzoni et al. at `k = 2`, and it is also the WhatsHap result.

## B. Exact and heuristic solvers

| Work | Model | Exact | Coverage handling | Status |
|---|---|---|---|---|
| WhatsHap (Patterson, Marschall, Pisanti, van Iersel, Stougie, Klau, Schönhuth 2015), [J. Comput. Biol. 22(6)](https://doi.org/10.1089/cmb.2014.0157) | wMEC, general bipartition | yes, fixed-parameter in coverage | random pruning in the 2015 paper | `paper` |
| HapCol (Pirola, Zaccaria, Dondi, Klau, Pisanti, Bonizzoni 2016), [Bioinformatics 32(11)](https://doi.org/10.1093/bioinformatics/btv495) | k-cMEC: wMEC with at most `k` corrections per column | yes, FPT in (`k`, coverage) | handles coverage 25 on 256 GB; says WhatsHap "cannot handle coverages higher than 20"; no read selection | `paper`, sections 1, 2 and the results summary |
| HapCHAT (Beretta, Patterson, Zaccaria, Della Vedova, Bonizzoni 2018), [BMC Bioinformatics 19:252](https://doi.org/10.1186/s12859-018-2253-8) | MEC with a per-column bound that the solver raises until feasible | yes, under the current bound | merges near-identical reads, then greedy selection on phred score with random tie-breaking | `paper` |
| HapCUT2 (Edge, Bafna, Bansal 2017), [Genome Research 27(5)](https://doi.org/10.1101/gr.213462.116) | maximum likelihood, greedy likelihood cuts; "provides no performance guarantees on the optimality" | no | none needed, linear in reads; prunes variants with posterior below 0.8 and splits blocks at poor linkages after solving | `paper`, main text and methods; supplement not read |
| ReFHap (Duitama, Huebsch, McEwen, Suk, Hoehe 2010, ACM BCB), also described in [Duitama et al. 2012, NAR 40(5)](https://doi.org/10.1093/nar/gkr1042) | Maximum Fragments Cut, NP-hard by reduction from max-cut (Theorem 1); solved by a max-cut heuristic | no | none needed | `paper` |
| ProbHap (Kuleshov 2014) | likelihood generalizing MEC, exact DP over read bipartitions | yes | DP "practical for coverages of up to 10-12x"; above that, merges reads likely to share a haplotype; no selection | `paper`, sections 1, 3.4, 4.4, 4.5, 5 |
| SDhaP (Das and Vikalo 2015), [BMC Genomics 16:260](https://doi.org/10.1186/s12864-015-1408-5) | MEC as correlation clustering, SDP relaxation with randomized rounding; diploid and polyploid | no | none needed | `paper`, methods |
| H-PoP and H-PoPG (Xie, Wu, Wang, Jiang 2016), [Bioinformatics 32(24)](https://doi.org/10.1093/bioinformatics/btw537) | polyploid balanced optimal partition, NP-hard; heuristic DP | no | "can deal with long reads and deep read coverage" | `paper`, abstract and introduction |
| MixSIH (Matsumoto and Kiryu 2013), [BMC Genomics 14(S2):S5](https://doi.org/10.1186/1471-2164-14-S2-S5) | two-component mixture model; a "minimum connectivity" score per segment used to keep reliable regions | no | none needed | `paper`, background and results |
| He, Choi, Pipatsrisawat, Darwiche, Eskin 2010, [Bioinformatics 26(12)](https://doi.org/10.1093/bioinformatics/btq215) | MEC by DP in `O(m 2^k n)`, `k` the longest read in variants; partial MaxSAT for long reads | yes | exponential in read length, not coverage; on HuRef it excluded 3725 reads, per Chen et al. | `paper`, abstract and sections 3.2, 3.3 |
| Chen, Deng, Wang 2013, [Bioinformatics 29(16)](https://doi.org/10.1093/bioinformatics/btt349) | MEC by integer linear programming per block, with or without the all-heterozygous assumption | yes | reports the first optimal MEC solutions for HuRef | `paper`, abstract, introduction and methods |

The dynamic program in `phasebudget.dp` is the WhatsHap algorithm, specialized.
It is reimplemented here because the project needs its internal tables, not
because the algorithm is new. Nothing in this repository claims it. Reading the
paper settled three things that the `recalled` row had wrong or vague.

**WhatsHap does not assume all-heterozygous sites.** Section 2 states the
assumption and then declines it: "Under the all-heterozygous assumption, where
all columns correspond to heterozygous sites, h1 must be the complement of h2.
Some methods need to make this assumption for computational reasons, and this
is why we mention it, although our tool does not need it." The project's
complementary-haplotype decision is therefore a restriction *relative to
WhatsHap*, not a shared premise, and `phasebudget.dp` is the WhatsHap DP
specialized to the all-heterozygous case rather than a faithful port. That
raises the stakes on OQ-7: the restriction is a divergence from the reference
algorithm and its cost has to be measured, not assumed small.

**The 2015 paper's state is an unlabelled bipartition and it exploits the
symmetry this project deliberately gives up.** The paper stores `2^(|F(j)|-1)`
entries per column "because C(j, (R, S)) = C(j, (S, R))". `phasebudget.dp`
carries a labelled bipartition plus the allele of `h` at `j`, which is
`2^(k+1)` states, twice the paper's count. That is a deliberate cost:
`docs/design/margin-extraction.md` needs the symmetry present to make the phase
relation well defined at a boundary.

**The 2015 paper downsamples at random.** Section 4: "we pruned the data sets
to target coverages of 5x, 10x and 15x by removing randomly selected reads that
violated the coverage constraints until no more such reads exist." Proper
selection is listed as unfinished work in the discussion: "Another possibility
is to just properly select reads up to the manageable coverage and to discard
the rest." So the deployed heuristic postdates the paper, and the paper itself
names coverage-constrained read selection as an open direction.

**No margin, confidence, or sensitivity quantity is computed from the DP
tables.** The word confidence appears twice in the paper and both times refers
to the per-entry sequencing weight, not to the firmness of a phase decision.
There is no backward pass. Block-wise margin extraction is not in WhatsHap.

HapCHAT is the nearest published thing to this project's third component and it
has to be distinguished carefully, because a reader will ask. HapCHAT is
adaptive in its *error bound*: it computes a per-column correction budget from
coverage and an error-rate parameter and raises it by
`z + floor(log2(z) + 1)` until the instance becomes feasible. It is not
adaptive in its read selection. Its coverage reduction is a one-pass pipeline,
read merging by a likelihood ratio followed by "a greedy selection approach
(based on the Phred score), with ties broken at random", and it does not
re-select after solving. It computes no margin from solver internals. So
HapCHAT occupies the name "adaptive haplotype assembly" while adapting a
different quantity, and it is evidence that the *threshold*-adaptive idea is
taken and the *selection*-adaptive one is, so far, not visibly.

## C. Read selection

| Work | Objective | Guarantee | Status |
|---|---|---|---|
| Fischer and Marschall 2016, [Selecting Reads for Haplotype Assembly](https://doi.org/10.1101/046771) | lexicographic dynamic coverage, static coverage, and worst-quality score; optional component bridging | heuristic | `paper` |
| WhatsHap as deployed, described in Martin et al. 2016, section 2.2 | the Fischer and Marschall heuristic, greedy on a per-read score with bridging reads; default cap 15x, compared only with random downsampling; "exceeding 15x coverage mostly leads to adding redundant information" | heuristic | `paper` for section 2.2 and the discussion; source at `1f92cf6` read |
| Mäkinen, Staneva, Tomescu, Valenzuela 2015, [Interval scheduling maximizing minimum coverage](https://arxiv.org/abs/1508.07820); journal version Discrete Appl. Math. 225 (2017) 130-135, with Wilzbach added as author | maximize the least per-position coverage under a cap | exact, max-flow plus floor search | `paper`, arXiv v2 of 30 Oct 2015, which is also the text the University of Helsinki repository deposits as the journal article (hdl 10138/307575); cite the journal version |
| Arkin and Silverberg 1987, [Scheduling jobs with fixed start and end times](https://doi.org/10.1016/0166-218X(87)90037-0), Discrete Appl. Math. 18, 1-8 | maximum-value subset of fixed intervals on `k` identical machines, which is the span-coverage selection problem | exact. The Lemma: the clique-by-job matrix has consecutive ones in each column, so it is totally unimodular. Theorem 1: min-cost flow in `O(n^2 log n)` that deletes the cheapest intervals so every clique has at most `k`, and one run solves every `k`. This is a different network from the one `select_optimal` builds | `paper` |
| Hoffman and Kruskal 1956 | the integrality theorem for totally unimodular matrices | theorem | `recalled`; book chapter, not obtained, and not needed: Schrijver's notes state it as Corollary 8.2a and Arkin and Silverberg prove the case used here |
| Schrijver, *A Course in Combinatorial Optimization*, lecture notes of 23 March 2017 | Corollary 8.2a (Hoffman-Kruskal) and Corollary 8.4c: every interval matrix, rows of consecutive ones, is totally unimodular | theorem | `paper`, sections 8.2 to 8.4. The book section 19.3 cited in `docs/design/selection-integrality.md` is not read |
| Karp 1972, Reducibility among combinatorial problems, in *Complexity of Computer Computations*, 85-103; read as the reprint in *50 Years of Integer Programming* | the Main Theorem lists CLIQUE (3) and NODE COVER (5) as complete; independent set follows because a set is independent exactly when its complement is a node cover | theorem | `paper`, sections 3 and 4 and the Main Theorem |
| Garey, Johnson and Stockmeyer 1976, Theor. Comput. Sci. 1, 237-267 | maximum independent set is NP-hard on cubic graphs | theorem | `recalled`; Elsevier and academia.edu only, not obtained, and nothing depends on it |
| Bar-Yehuda, Halldórsson, Naor, Shachnai, Shapira 2006, [Scheduling split intervals](https://csaws.cs.technion.ac.il/~reuven/PDF/BarHal04.pdf), SIAM J. Comput. 36(1) | maximum-weight independent set in `t`-interval graphs | APX-hard for `t >= 2` even on unweighted (2,2)-union graphs (Corollary 2.3); `2t`-approximation | `paper`, sections 1 and 2 read |

Reading the WhatsHap paper split this section's first two rows apart. The 2015
paper prunes at random and lists proper selection as an open direction, so
"WhatsHap's selection" is not one object: the published algorithm downsampled
randomly, and the Fischer and Marschall heuristic became the deployed procedure
afterwards. Any sentence in the report comparing against "WhatsHap's read
selection" has to say which of the two it means. The quotations are in
section B.

HapCHAT is a third published selection procedure and it is neither of these:
greedy on phred score with random tie-breaking, preceded by a read-merging step
that has no analogue in this repository. It is a fourth baseline this project
does not yet implement.

The last three rows arrived with `docs/design/coverage-convention.md`, which
classifies the called-allele coverage convention. The first two support the
hardness of that convention and are textbook results rather than prior art for
this project; the reduction itself is self-contained and tested, so nothing
depends on the rows being read. The split-interval row is different: OQ-11
asks whether one hole per read is already hard, that is the `d = 2` case, and
the question cannot be called open until the row is at `paper`. No claim about
the one-hole case may be made before then.

The Mäkinen et al. objective is a maximum-minimum coverage objective. It is
not the additive total-weight objective solved by `select_optimal`, and its
existence does not establish the flow reduction for that objective. Rows four
and five are the reason the exact selector in this repository is not
a contribution. Maximum-weight `k`-colourable interval selection is a solved
problem with a textbook flow algorithm, and the phasing instance is an instance
of it. What is open is whether anyone has connected the two, and what the
outcome difference between the deployed heuristic and the additive selector is.
Both are questions the review has to answer, not premises.

## D. Sensitivity, margins, and adaptivity

| Work | Relation | Status |
|---|---|---|
| Kohli and Torr 2008, [Measuring uncertainty in graph cut solutions](https://doi.org/10.1016/j.cviu.2008.07.002), CVIU 112(1); read as the Brookes preprint | min-marginal energies, equations 5 and 6: the least energy with one or two variables fixed. The margin here is a pairwise min-marginal minus the optimum. Exact max-marginals come from max-product message passing on tree-structured graphs, and a chain is a tree | `paper`, sections 1 and 2 |
| Settles 2010, Active Learning Literature Survey, Wisconsin CS TR 1648 | margin sampling (Scheffer et al. 2001): query the instance whose two most probable labels are closest; applied to HMMs and CRFs | `paper`, sections 3.1 and 4 |
| Weilguny et al. 2023, [BOSS-RUNS](https://doi.org/10.1038/s41587-022-01580-z), Nat. Biotechnol. 41 | per-site genotype posteriors decide, during nanopore sequencing, which molecules to read in full; data go to the sites of highest uncertainty | `paper`, results and discussion |
| Confidence scores in phasing tools | HapCUT2 prunes variants with posterior below 0.8 and splits blocks where the posterior of the adjacent configuration is low; ProbHap computes a per-adjacent-pair transition score by forward-backward; both use it after solving, to prune or split; see section F | `paper` for both |
| HapCHAT 2018 | adaptive, but in the per-column error bound rather than in read selection | `paper` |
| WhatsHap 2015 | computes no phase-decision confidence from its DP tables, and has no backward pass | `paper` |

This is the row that decides whether the third component of the project is new.
The specific object here is the MEC margin of an adjacent phase decision,
extracted block-wise from the DP tables, used to reweight reads for a second
exact selection. Each of the three parts is standard on its own. The review has
to establish whether the composition is, and the honest outcome may be that it
is not.

Two of the five rows are now read, and both are negative for the nearest
alternatives: WhatsHap derives no phase confidence from the tables that would
support it, and HapCHAT adapts a threshold rather than a selection. That does
not answer OQ-1. It removes the two candidates most likely to have answered it
already, which is a different and much weaker statement.

What was searched on 2026-09-10, so that the gap in the search is on the record
rather than implied: full-text reads of WhatsHap 2015 and HapCHAT 2018; keyword
searches combining read selection, haplotype assembly, solver-guided, iterative
reweighting, sensitivity, margin, active learning, and coverage-constrained
optimization. Those searches surfaced HapCHAT and RFhap and nothing that
composes constrained re-optimization with read selection. Three things that
would change the picture were **not** searched: the LP-sensitivity and
active-learning literatures on their own terms rather than through phasing
keywords, the polyploid and metagenomic assembly literatures where the same
composition could appear under other names, and the WhatsHap issue tracker and
release notes, where a deployed heuristic can change without a paper. Until
those are done, OQ-1 stays open and no novelty claim rests on this section.

All three were searched on 25 September 2026; section G has the results. The
first three rows of the table replace the two `recalled` rows that stood here.

## E. Evaluation data

| Resource | Use | Status |
|---|---|---|
| GIAB HG002 v4.2.1 benchmark, [README](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/README_v4.2.1.txt) and files read 25 Sep | the main benchmark VCF is unphased: every heterozygous call in the first 17,054 checked is `0/1` with no phase set. Phased versions are in `GRCh38/SupplementaryFiles/` | source read |
| Wagner et al. 2022, [Benchmarking challenging small variants with linked and long reads](https://doi.org/10.1016/j.xgen.2022.100128), Cell Genomics 2(5) | the v4.2.1 paper. Two phased HG002 files: `..._phased_MHCassembly_StrandSeqANDTrio.vcf.gz`, conservative, keeping calls where trio phasing and Strand-seq plus HiFi phasing with WhatsHap agree; and `..._hifiasm_v11_phasetransfer.vcf.gz`, phase transferred from a trio-hifiasm assembly, "draft phased vcfs, which have not been evaluated for accuracy" | `paper`, phasing methods |
| Zook et al. 2016, [Extensive sequencing of seven human genomes](https://doi.org/10.1038/sdata.2016.25), Sci. Data 3:160025; read as the 2015 bioRxiv preprint | data from 11 technologies for seven genomes, including the Ashkenazi trio of which HG002 is the son; cell lines, DNA and data are public | `abstract`; the data index at `giab_data_indexes` is read and lists the HG002 datasets |
| GIAB HG002 v5.0q, [README](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/NIST_HG002_v5.0q_variant-benchmarksets_README.md), read 1 Oct 2026 | assembly-based small and structural variant draft benchmarks from GIAB (PI Justin M. Zook), made with the DeFrABB pipeline and "based on the V1.1 release of the T2T HG002 Q100 assembly"; the README asks users to cite the Q100 preprint until the benchmark itself is published. The run uses `HG002_GRCh38_v5.0q_smvar.vcf.gz` inside `HG002_GRCh38_v5.0q_smvar.benchmark.bed` | source read |
| ONT open data, release `giab_2023.05`, [bucket](https://ont-open-data.s3.amazonaws.com/giab_2023.05/), `analysis/variant_calling/hg002_sup_60x/hg002.haplotagged.bam` | HG002 R10.4.1 reads basecalled in super-accuracy mode, aligned to GRCh38 and haplotagged by the release's workflow (`docs/design/ont-run-plan.md`) | path read; no paper |
| WhatsHap `compare`, [user guide](https://github.com/whatshap/whatshap/blob/main/doc/guide.rst) | switch errors between adjacent phased variants in a block; two consecutive switches reported as one flip; Hamming distance; either haplotype may be matched | documentation read |

Chromosome 20 is chosen because it is the conventional target, which makes the
numbers comparable to published ones; Damaraju et al. 2026 use HG002 and HG005
in the same way. That the target is conventional is still an impression and not
a counted fact, so the report says only that HG002 is a GIAB reference sample.
The dataset is not yet fixed and OQ-9 depends on it; OQ-12 carries the choice of
truth file. In a sample of 16,779 heterozygous calls from the start of chromosome
1, the conservative file phases 74.6 per cent and the hifiasm file 99.3 per cent.
Simulated truth is complete, so `phasebudget.metrics` has never had to skip an
unphased truth pair. The real-data evaluation has to compare only adjacent
pairs that the truth phases in one phase set, as `whatshap compare` does; with
that rule the conservative file costs coverage of the evaluation and not
correctness.

## F. Re-check of 24 September 2026

A second search, run claim by claim after the scope was sharpened to the
measurement. Keyword searches, citation lists of the three nearest papers
(Semantic Scholar for Fischer and Marschall 2016 and HapCHAT 2018; OpenAlex
and Crossref for Mäkinen et al.), full-text reads where a PDF was reachable,
and the current WhatsHap source. What it found changes three components out of
four.

| Work | What it does | Bearing | Status |
|---|---|---|---|
| Kuleshov 2014, [Probabilistic single-individual haplotyping](https://www.cs.cornell.edu/~kuleshov/papers/eccb2014.pdf), Bioinformatics 30(17) (ProbHap) | forward-backward over the read-bipartition DP; a per-SNP "transition score", the probability that the phase relative to the previous SNP is right; breaks blocks where it is low; states that with MEC factors the same DP gives "exact solutions to the MEC objective" and "confidence probabilities associated with the MEC"; merges reads when coverage is too deep | the sum-product form of this project's margin extraction. Margins from forward and backward DP tables are not new | `paper`, sections 3.1, 4.4, 5 read |
| HapCUT2 (Edge, Bafna, Bansal 2017), [output format](https://github.com/vibansal/HapCUT2/blob/master/outputformat.md) | phred-scaled per-variant switch quality; blocks split at poor linkages | a per-boundary phase confidence, used to prune and split, not to reselect | `paper` and documentation; see section B |
| MixSIH (Matsumoto and Kiryu 2013) | named by Kuleshov as the only earlier method with confidence scores; its minimum connectivity score selects reliable regions, not reads | same family | `paper`, see section B |
| Damaraju et al. 2026, [Determinants of haplotype phasing accuracy](https://www.biorxiv.org/content/10.64898/2026.05.04.722832v2), bioRxiv | coverage, read length and reference effects on WhatsHap-style phasing on HG002 and HG005; PhaseQuality, per-variant-pair confidence tiers "based solely on sequencing data"; error plateaus near 20x | real-data neighbour for component 1 and for margins; no read selection | `paper`, results, PhaseQuality and methods |
| Hosseini et al. 2025, [pHapCompass](https://arxiv.org/abs/2512.04393), arXiv | distribution over phasings with read-assignment uncertainty | uncertainty quantification, no selection | `paper`, abstract and sections 1 to 2; pHapCompass-short chooses SNP pairs greedily, not reads |
| Majidian, Kahaei, de Ridder 2020, [MEC-based haplotype assembly: considerations for long read data](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0234470), PLOS One | when a wrong haplotype has lower MEC than the true one, by theorem (errors in a column above half its depth) and by simulation across coverage 2 to 100, error rate, and Illumina, PacBio and ONT; coverage 25 needed for PacBio RS | component 1's question, already asked and answered in unweighted form with a heuristic solver. Their Theorem 2 is an upper bound, not the probability it is stated as; see section G | `paper` |
| Duitama et al. 2012 (online 2011), [NAR 40(5)](https://doi.org/10.1093/nar/gkr1042) | on real NA12878 fosmid data the trio gold standard has higher MEC than HapCUT in every block; HapCUT's MEC correlates 0.11 with switch error, and within bins of allele error rate the correlation is negative, -0.5 and -0.4 | the objective-truth gap on real data, as a correlation | `paper`, the MEC section |
| Fischer and Marschall 2016, discussion | "We are currently comparing our heuristic approach to the flow-based algorithm proposed by Mäkinen et al. (2015)." | the comparison component 2 makes was planned by the heuristic's authors | `paper` |
| Martin et al. 2016, [WhatsHap: fast and accurate read-based phasing](https://www.biorxiv.org/content/10.1101/085050v2.full), bioRxiv | read selection compared only with random downsampling; default cap 15x | no optimal or flow comparison | `paper`, section 2.2 and discussion |
| Garg 2018, Computational Haplotyping: Theory and Practice, PhD thesis, Saarland (Marschall group) | "we prune the input data sets to a target maximum coverage using the read-selection method introduced by Fischer and Marschall (2016)", 5x for PedMEC-G and 5x or 15x for wMEC; full text searched for Mäkinen, flow, and read selection | no comparison with a flow or exact selector | `paper`, searched; the pruning passage read |
| Citations of Mäkinen et al. | Crossref counts 4; OpenAlex lists Fischer and Marschall 2016 and one strain-composition paper | no later paper runs the planned comparison | searched |
| Citations of Fischer and Marschall | six on Semantic Scholar: WhatsHap 2016, pedigree phasing, HapCHAT, two diplotyping papers, Garg 2018 | none compares selectors against an optimum | searched |
| WhatsHap source, commit `1f92cf6`, 22 Sep 2026 | `readselect.pyx` is still the slice-and-bridge heuristic; coverage counted over `range(begin, end)` | the deployed selector is the one benchmarked, and it uses span coverage | source read |
| Carlisle and Lloyd 1995, [On the k-coloring of intervals](https://doi.org/10.1016/0166-218X(95)80003-M), Discrete Appl. Math. 59, 225-235 | maximum-weight `k`-colourable interval set in `O(k S(n))`, Theorem 3. The network of Lemma 4 is a chain of capacity-`k` clique edges of cost zero and one capacity-1 edge of cost minus the weight per interval, with a flow of size `k` | the network `select_optimal` builds, with the same bound. It is not a faster algorithm than the one here; it is the one here, and the design document should attribute it | `paper`, sections 1 and 3 |

### What each component can still claim

**Component 1, the objective-truth gap.** The question is published. Majidian
et al. 2020 measure when MEC prefers a wrong haplotype, across coverage and
error rate, and give the column-majority theorem that explains it. What this
repository adds is narrower: weighted MEC, an exact optimum rather than
HapCUT, the gap as a magnitude rather than a yes-or-no, and the link to margin
informativeness. The report must cite Majidian et al. as the prior
measurement and present component 1 as an extension of it.

**Component 2, the certified optimality gap.** Nothing found reports it. The
closest thing is a sentence of intent in Fischer and Marschall 2016 with no
visible follow-up in the citing literature, the thesis from the same group, or
the current source. That is strong negative evidence and it is still negative
evidence; the novelty gate forbids resting a deliverable on it, and the
measurement does not need to rest on it. The one way to turn it positive is to
ask the WhatsHap authors whether the comparison was run.

**The gapped-coverage theory.** Settled in the literature. The one-hole case
is APX-hard by Bar-Yehuda et al. 2006 through an embedding that
`tests/test_gapped.py` checks. There is no theorem here for the project to
own; the contribution is recognizing that the convention decides tractability.

**Component 3, margins and adaptive reselection.** The margin as a quantity is
taken: ProbHap computes the sum-product version from the same DP tables, and
HapCUT2 and MixSIH report per-boundary confidence. Using such a confidence to
reweight reads and reselect under a coverage cap was not found in any of the
sources above, all of which use the confidence only to prune or split. That is
the one composition still unclaimed, and it is exactly as strong as the search
behind it.

## G. Completion of 25 September 2026

Every row that could be obtained was read against its full text, and the three
searches that section D listed as missing were run. PDFs of everything read were
kept in the development repository and are not redistributed here. Text was extracted with
`pdftotext` or from the Europe PMC full-text XML, and a row marked `paper`
names the sections read where the whole text was not.

### Still `recalled`, and why

| Row | Reason | What stands in for it |
|---|---|---|
| Hoffman and Kruskal 1956 | book chapter, not obtained | Schrijver's lecture notes, Corollary 8.2a; the Lemma of Arkin and Silverberg for the case used here |
| Garey, Johnson and Stockmeyer 1976 | Elsevier and academia.edu only | nothing depends on it; the reduction in `coverage-convention.md` is tested |
| Scheffer, Decomain, Wrobel 2001, Active hidden Markov models for information extraction, IDA 2001, LNCS 2189 | Springer, not obtained | the definition of margin sampling as quoted by Settles 2010, section 3.1 |

Every other row is `paper` or `abstract`. The second pass, later on 25
September with PDFs supplied by the author, added Lancia et al. 2001, Lippert
et al. 2002, Greenberg et al. 2004, ReFHap 2010, Carlisle and Lloyd 1995,
Karp 1972, Garg 2018, and the Zook et al. abstract.

### Two attributions the second pass corrected

**MEC is Lippert et al. 2002, and wMEC is Greenberg et al. 2004.** The matrix
credited Lancia et al. 2001 with formulating error correction. Lancia et al.
define MFR, MSR and LHR and not MEC, and Cilibrasi et al. say MEC is in that
paper only "implicitly". Lippert et al. define MEC and prove it NP-hard from
MAXCUT; Greenberg et al. add the per-entry weight proportional to base-call
confidence, which is this project's objective. HapCol attributes both the same
way. The report cites Lippert et al. for MEC and Greenberg et al. for wMEC.

**The selection network is Carlisle and Lloyd's.** `select_optimal` builds a
chain of capacity-`c` slack arcs and one unit arc of cost `-w_r` per read, and
pushes `c` units; that is the network of Lemma 4 of Carlisle and Lloyd 1995,
whose Theorem 3 gives the `O(c S(n))` bound the implementation also meets with
a heap. Arkin and Silverberg solve the same problem with a different network
that deletes intervals, in `O(n^2 log n)`. An earlier explainer called Carlisle and
Lloyd a faster algorithm than the one here; it is the same algorithm. So the
problem's tractability is Arkin and Silverberg's, the network and the bound are
Carlisle and Lloyd's, and the repository claims neither.

One integrality detail for `selection-integrality.md`. Schrijver's interval
matrix has consecutive ones in each row. The selection constraint matrix has a
row per variant and a column per read, and the ones are consecutive down each
column, which is the transpose. Total unimodularity is closed under
transposition, and Arkin and Silverberg state the column form directly, so the
argument holds; the design document should cite the Lemma of Arkin and
Silverberg rather than the unread book section.

### Majidian et al. Theorem 2 is an upper bound

Theorem 1 of Majidian, Kahaei and de Ridder 2020 is a sufficient condition for
the true haplotype to have higher MEC than some other haplotype: a column with
more than half of its entries wrong. Theorem 2 computes the probability that
no column meets it and calls that the probability of obtaining a minimum MEC
for the true haplotype. The step between them is the sentence "if the
antecedent is not satisfied, the MEC approach works properly", which is not
proved and is false. `benchmarks/benchmark_majidian_bound.py` carries a
five-read instance in which every column is less than half wrong and `00010`
costs 2 against the truth's 3; `tests/test_benchmark_majidian_bound.py` checks
it against both exact solvers. Their equation 7 is therefore an upper bound on
the probability that the truth is an MEC optimum.

The same script measures how loose the bound is. Unweighted MEC, exact dynamic
program, 30 variants, span 5, 100 replicates per cell:

| cov | err | Theorem 2 | no bad column | truth optimal | no bad column, truth not optimal |
|---:|---:|---:|---:|---:|---:|
| 4 | 0.05 | 0.832 | 0.860 | 0.760 | 0.100 |
| 4 | 0.10 | 0.588 | 0.660 | 0.490 | 0.170 |
| 4 | 0.20 | 0.149 | 0.120 | 0.020 | 0.100 |
| 6 | 0.05 | 0.924 | 0.920 | 0.860 | 0.060 |
| 6 | 0.10 | 0.767 | 0.750 | 0.610 | 0.150 |
| 6 | 0.20 | 0.284 | 0.310 | 0.070 | 0.240 |
| 8 | 0.05 | 0.960 | 0.950 | 0.940 | 0.010 |
| 8 | 0.10 | 0.867 | 0.850 | 0.800 | 0.050 |
| 8 | 0.20 | 0.432 | 0.420 | 0.140 | 0.280 |

The formula agrees with the event it describes to within sampling error, which
checks the implementation. It overstates the event that matters by a factor of
more than 7 at error 0.20 and coverage 4, and by 3 to 4 at coverages 6 and 8. The
limits: simulated reads of fixed span with uniform independent errors, and
coverage below the 10 to 100 range of their figures, because the exact solver
is exponential in coverage. The direction of the error is a theorem and holds
at any coverage; its size is measured only here.

### What was searched, and what it found

- **The WhatsHap issue tracker and changelog.** Searches for read selection,
  max-coverage, flow, Mäkinen, downsampling, confidence and switch quality.
  No issue or release discusses an optimal or flow-based selector. Issue 170
  has a maintainer writing that raising `--max-coverage` "does not effect the
  result very much"; issue 225 notes the cap is shared across a pedigree.
  Release 2.2 of 26 January 2024 adds `--algorithm=heuristic`, which trades
  MEC optimality for higher coverage. Issue 160 (2018) proposed an
  `--algorithm` switch with "autoselect" as a possible value.
- **Sensitivity on its own terms.** Kohli and Torr define min-marginal
  energies and use them as the confidence of a MAP labelling. The margin here
  is a pairwise min-marginal gap on a chain model, computed by forward and
  backward passes. It is the max-product counterpart of ProbHap's transition
  score and an instance of Kohli and Torr's quantity.
- **Active learning and adaptive acquisition.** Margin sampling, as surveyed by
  Settles, picks the instances whose two best labels are closest, including
  for HMMs. BOSS-RUNS steers nanopore sequencing toward the sites of highest
  genotype uncertainty. Both use an uncertainty signal to decide which data to
  use. Neither concerns phase, neither chooses among reads already sequenced,
  and neither works under a coverage cap.
- **Polyploid and metagenomic assembly.** H-PoP, SDhaP, pHapCompass (which
  greedily chooses SNP pairs to add to its current phasing, not reads) and
  Floria (beam search over reads in windows, then a flow over blocks) were
  read. snpphaser (Zhang, Nie and Wang, bioRxiv 2025) phases high-quality
  variants first and filters the rest against them; the abstract is read and
  it selects variants, not reads. The full text is now read: it phases
  high-quality SNPs, then "selectively retains only high-confidence low-quality
  SNPs". That is the nearest structural analogue found to component 3, a
  second pass steered by a first solution, and it chooses variants rather than
  reads and has no coverage cap. nTChap (BMC Genomics, 2026; PMC12888297) is
  read: iterative two-round read clustering and consensus, with no selection
  and no confidence feedback. RFhap (Gonzalez et al., bioRxiv 2026), which the
  10 Sep search surfaced, bins reads by parent with trio k-mers before de novo
  assembly and is outside read-based phasing.
- **The real-data neighbour.** Damaraju et al. 2026 run Clair3 and WhatsHap on
  HG002 and HG005, find error plateaus near 20x, and assign confidence tiers
  to variant pairs with a random forest on alignment features. Their GRCh38
  truth is built from the parents' calls by Mendelian inheritance, which is a
  third choice of truth for OQ-12.

### What each component can claim after the review

**Component 1.** An extension of Majidian et al. 2020, with Duitama et al.
2012 as the real-data precedent: the truth has higher MEC than a heuristic
optimum on real reads, and the correlation between MEC and switch error is
weak or negative. The exact measurement adds one result neither has: in 1
to 28 per cent of all simulated instances no column is more than half wrong and
MEC still prefers a wrong haplotype, so the closed form is not a substitute for
solving.

**Component 2.** Unchanged. Nothing read reports the optimality gap of a
deployed selector, and the measurement does not depend on that absence.

**Component 3.** Every ingredient is published. The margin is a min-marginal
gap (Kohli and Torr), the per-boundary confidence is in ProbHap and HapCUT2,
and steering data toward low-margin decisions is margin sampling. The
combination found nowhere in the sources read is margin-guided reselection
from reads already sequenced, under a fixed coverage cap, to improve phasing.
The report states component 3 as margin sampling applied to that problem, and
its value rests on the measured accuracy change, not on novelty.

## H. Re-check of 27 September 2026, after the genome-wide run

Run before any claim from the genome-wide run goes into a preprint. The
searches: read selection and coverage in long-read phasing, 2025 and 2026;
HiPhase, LongPhase and their benchmarks; exact MEC work since 2025; phase
confidence and calibration; methylation-aware phasing; ONT R10 phasing on
HG002; RECOMB 2026 accepted papers; information-theoretic phasing; polyploid
and pedigree coverage caps; convex-cost flow for interval selection; the
HG002 Q100 benchmark.

Status in this section is stricter than the table key. Pages were read
through a fetch tool that returns a summary of the page, not the page. A
row read that way is `abstract` even where the summary drew on the full
text, and every number quoted here must be checked against the paper before
it is cited. A row known only from a search-result snippet is `recalled`.

### New rows

| Work | What it does | Bearing on this project | Status |
|---|---|---|---|
| Holt et al. 2024, [HiPhase](https://doi.org/10.1093/bioinformatics/btae042), Bioinformatics 40(2) btae042; [PMC10868326](https://pmc.ncbi.nlm.nih.gov/articles/PMC10868326/) | wMEC variant solved by A* search per putative block, "without down-sampling" at 30x HiFi; on HG002 replicates NG50 480 kb and 929 switchflips against WhatsHap's 407 kb and 1,615; says WhatsHap "down-samples the data in order to run efficiently, potentially biasing the data" | The strongest objection to the cap framing: a production HiFi phaser runs uncapped. Main text read 28 Sep from the Europe PMC full-text XML: "no down-sampling of the data"; A* "to search for the optimal diplotype" with a heuristic "generated by chaining sub-problem solutions". Whether that heuristic is admissible, and so whether the search is exact, is in the supplement, which is not read | `paper`, main text; supplement not read |
| Lin et al. 2022, [LongPhase](https://doi.org/10.1093/bioinformatics/btac058), Bioinformatics 38(7); [v2.0 release notes](https://github.com/twolinin/longphase/releases) | graph-based SNV and SV co-phasing, no coverage cap; v2.0 adds 5mC co-phasing and reports switch error about 0.055 per cent on ONT | Uncapped, fast, widely used; the baseline a reader expects beside WhatsHap. Preprint (bioRxiv 2021.09.09.459623 v1) searched on 28 Sep: phasing is a greedy construction of the longest pair of disjoint paths through a graph over adjacent SNPs and SVs; reads are downsampled only to build the evaluation sets (60x, 50x, 25x); no read selection. The journal version is not open access | `paper` for the preprint, searched: the algorithm and data passages read; release notes `abstract` |
| LongPhase 2 authors, [manuscript repository](https://github.com/ythuang0522/longphase-2-paper) (GitHub user ythuang0522), in preparation, read 27 Sep 2026 | co-phases SNVs, indels, SVs and 5mC; a graph neural network (GATv2 plus transformer) scores phased variants and unphases those predicted misphased; trained on HG002 ONT at 10 to 60x | Direct competition for any learned or calibrated per-boundary confidence (T4, T9). Not a publication; no preprint found | `abstract` (repository README) |
| Pfennig and Akey 2026, [LongHap](https://www.biorxiv.org/content/10.64898/2026.03.11.710820v2), bioRxiv | DAG of allele co-occurrence, 5mC-informed read assignment, Viterbi-like decoding; no downsampling mentioned. HG002 Revio HiFi at 38x: switch error WhatsHap and HapCUT2 about 0.246 per cent, LongPhase 0.192, LongHap 0.196 | Genome-wide switch error of current tools on HG002 is two orders of magnitude above this project's pooled rates, because the denominators differ (their variant sets and truth against SNVs inside the v4.2.1 BED here). Any comparison must use one protocol | `abstract` |
| HapBridge, [bioRxiv 2025.11.07.687303](https://www.biorxiv.org/content/10.1101/2025.11.07.687303v4) | methylation-guided correction of switch errors and bridging of blocks after SNV phasing; ONT R9, R10 and HiFi | A post-phasing repair, like two-level orientation, with a different signal | `recalled` (search snippet) |
| Zhang, Nie and Wang 2025, snpphaser, [bioRxiv 2025.12.21.695866](https://www.biorxiv.org/content/10.64898/2025.12.21.695866v1) | already in section G. New numbers: HG002 switch error ONT ultra-long 0.22 per cent (WhatsHap 0.32, LongPhase 0.23), HiFi 0.28 (WhatsHap 0.43, LongPhase 0.30); drops links spanned by fewer reads than a threshold | Uses spanning-read count per boundary as a reliability signal, which is the quantity the direction in the research plan (development repository) section 13 optimizes | `abstract` for the new numbers |
| Damaraju et al. 2026, bioRxiv preprint indexed in PMC ([PMC13317610](https://pmc.ncbi.nlm.nih.gov/articles/PMC13317610/), PubMed 42388792); corrected 30 Sep: not peer reviewed, and no journal version was found | as in section F; depth downsampling by `samtools --subsample` only, never WhatsHap's own cap; tools Clair3 with WhatsHap, LongPhase, PMDV at 10 to 30x; PhaseQuality tiers capture 82 to 99 per cent of errors | Confirms the gap: nobody in this line varies the selector or its cap. The PMC record is the preprint itself (verification report of 30 Sep) | `paper` for the bioRxiv text (section F) |
| Hansen et al. 2026, [A complete diploid human genome benchmark for personalized genomics](https://doi.org/10.1016/j.cell.2026.06.016), Cell 189(16):4857-4875.e31, 65 authors, first Nancy F. Hansen, last Adam M. Phillippy; preprint [bioRxiv 2025.09.21.677443](https://doi.org/10.1101/2025.09.21.677443), PMC12458380; [Q100 page](https://genomeinformatics.github.io/Q100/) | the T2T-HG002 Q100 v1.1 diploid assembly, telomere to telomere for all 46 chromosomes, with no detectable errors across 99.4 per cent of the diploid genome; adds 701.4 Mb of autosomal sequence and both sex chromosomes absent from v4.2.1; describes using the assembly to build phased GIAB variant benchmarks against GRCh37, GRCh38 and CHM13. The preprint does not use the name v5.0q; GIAB's v5.0q README says v5.0q is "based on the V1.1 release of the T2T HG002 Q100 assembly" and names this preprint as the citation | The citation for the assembly behind GIAB v5.0q, the primary HG002 truth. Cite it for the assembly and the README (section E) for what v5.0q is | `paper` for the preprint (Europe PMC JATS XML of PMC12458380, read 1 Oct 2026: abstract, introduction, the comparison with v4.2.1, the variant-benchmark and phasing-error passages, the discussion; methods and supplement not read); `abstract` for the Cell version (Europe PMC, PMID 42561913) |
| Zhang et al. 2026, [QHap](https://arxiv.org/abs/2603.25762), arXiv | phasing as Max-Cut, solved by a GPU simulated-bifurcation heuristic | Heuristic; no bearing on exactness or selection | `abstract` |
| Barak, Gibney and Jain, WABI 2025, [LIPIcs.WABI.2025.4](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.WABI.2025.4); journal version Algorithms Mol. Biol. 2026 | haplotype-aware long-read error correction built on MEC; NP-hardness of the de novo formulation | Shows MEC-based formulations still reach WABI and AMB in 2025 and 2026 | `recalled` (search snippet) |
| Kamath, Şaşoğlu and Tse 2015, [Optimal haplotype assembly from high-throughput mate-pair reads](https://doi.org/10.1109/ISIT.2015.7282588), ISIT, pp. 914-918; [arXiv 1502.01975](https://arxiv.org/abs/1502.01975) | "a simple formula for the coverage needed for haplotype assembly, under a generative model", via connections with decoding convolutional codes (arXiv abstract) | The information-theoretic prior art for T7 | `abstract` (Crossref fields and arXiv abstract read 30 Sep, re-read 1 Oct) |
| Chen, Kamath, Suh and Tse 2016, [Community recovery in graphs with locality](http://proceedings.mlr.press/v48/chena16.pdf), ICML, PMLR 48, 689-698 | pairwise noisy same-or-different measurements between nearby nodes; nearly linear algorithms that reach the information limit for exact recovery; phasing is the motivating example | T7's natural frame: a boundary is a local pairwise measurement. Any error model in section 13 must be stated against this | `abstract` |
| Vikalo group, ISIT 2014 (authors not retrieved), [Haplotype assembly: an information theoretic view](https://users.ece.utexas.edu/~hvikalo/pubs/isit14.pdf), ISIT | necessary and sufficient read counts for haplotype recovery with and without errors | Same | `recalled` (search snippet) |
| Mahajan, Jain and Kashyap 2024, [On the coverage required for diploid genome assembly](https://arxiv.org/abs/2405.05734), ISIT | information-theoretic coverage and read-length conditions for diploid assembly; practical assemblers need far more than the bound | Assembly, not phasing; a model of how such a result is stated | `abstract` |
| pHapCompass (Hosseini et al.) | already in section F; now in the RECOMB 2026 proceedings | Probabilistic phasing with uncertainty reached RECOMB | `paper` (arXiv, section F) |
| WhatsHap polyphase 2020 (first author from memory: Schrinner), [Genome Biology](https://doi.org/10.1186/s13059-020-02158-1) | polyploid phasing by read clustering and haplotype threading | A regime where exact wMEC is out of reach and coverage caps bind hardest | `recalled` (search snippet) |
| Separable convex-cost min-cost flow, e.g. [Végh 2012, arXiv 1110.4882](https://arxiv.org/abs/1110.4882) | integer min-cost flow with separable convex arc costs is polynomial; piecewise-linear convex costs reduce to parallel linear arcs | The algorithmic basis of the objective proposed in the research plan (development repository) section 13; textbook material, claimed by nobody here | `abstract` |

### What changed

- **The cap framing is weaker than on 25 September.** HiPhase, LongPhase,
  LongHap and snpphaser all phase HiFi or ONT at full depth without
  selecting reads. For a single diploid sample the question "which reads
  under a cap" matters only to WhatsHap and to exact solvers. It still binds
  hard for pedigrees (PedMEC is exponential in the summed coverage), for
  polyploids, and wherever an exact certificate is wanted.
- **Component 2 still stands.** No work found varies WhatsHap's selector
  or its cap and measures the result. Damaraju et al. downsample sequencing
  depth only. The G1 measurement (0.933 of the certified optimum,
  genome-wide) remains unreported in anything read. As before, no
  deliverable depends on that absence.
- **The truth is now behind the field.** The genome-wide run scored
  against GIAB v4.2.1 regions with a transmission truth. GIAB v5.0q and the
  Q100 assembly give a near-complete phase truth including difficult
  regions, where Damaraju et al. find most errors. A preprint scored only on
  v4.2.1 would be asked why.
- **Denominators differ from the tool papers.** Published HG002 switch
  errors are 0.1 to 0.4 per cent; the pooled rates here are about 0.003
  per cent. Both can be right, because the variant sets, regions and truth
  differ, but no sentence may compare them directly.
- **Learned confidence is crowded.** PhaseQuality (random forest tiers) is
  published and LongPhase 2 (graph neural network) is being written up. T9
  and a calibration-only paper would enter a crowded field.

### Verification of 30 September 2026

`2026-09-30-verification-and-sweep.md` records a second pass: bibliographic
fields checked against Crossref and Europe PMC, and a sweep for 2023 to 2026
work. Changes to this matrix:

- **Damaraju et al. 2026** is a bioRxiv preprint. The row above is
  corrected.
- **Q100.** The Q100 paper is Hansen et al., *Cell* 189(16):4857-4875.e31,
  2026, doi 10.1016/j.cell.2026.06.016. Its content status stays
  `recalled`: the benchmark's construction must be read before it is
  described.
- **The Vikalo-group row** is Si, Vikalo and Vishwanath, ITW 2014, not
  ISIT, with a journal version in *IEEE Trans. Inf. Theory* 63(6), 2017;
  status `abstract`.
- **Kamath, Şaşoğlu and Tse 2015** (ISIT, pp. 914-918, doi
  10.1109/ISIT.2015.7282588) rises to `abstract`.

New rows:

| Work | What it does | Bearing | Status |
|---|---|---|---|
| Battistella et al. 2025, Ralphi, *Genome Res.* 35(12):2617-2625, doi 10.1101/gr.280569.125 | graph-based deep RL for haplotype assembly; compares against WhatsHap at its default 15x cap on HG002 and reuses WhatsHap's read selection for preprocessing | shows WhatsHap's selector is reused downstream; the paper varies neither the selector nor the cap | `paper`, searched |
| Teo, Squair, Courtine and Skinnider 2024, *Nat. Commun.* 15:8805, doi 10.1038/s41467-024-53089-5 | a benchmark of scATAC-seq methods published as a Registered Report, with post hoc analyses reported separately | precedent for a pre-registered methods benchmark with the registered and exploratory split this project uses | `paper`, searched |
| Mallona et al. 2025, *PLoS Comput. Biol.* 21(11):e1013658 | perspective on continuous benchmarking; recommends pre-registration when developers benchmark their own methods | the rationale for this project's registrations | `paper`, searched |
| Ben-Yashar 2023, *Theory and Decision* 94(1):83-95, doi 10.1007/s11238-022-09872-1 | simple majority rule with an even number of homogeneous voters | nearest published treatment of the even-count identity; whether it states the identity is unknown from the abstract | `abstract` |

### Verification of 1 October 2026

Mock review 4 asked for the v5.0q paper to be read and for three cited
works to have rows. `2026-10-01-followup.md`, section 5, records what was
read. Changes:

- **Q100** (Hansen et al. 2026) is now `paper` for the bioRxiv preprint and
  `abstract` for the Cell version; the row in section H is rewritten. The
  preprint names no "v5.0q"; the link from the assembly to v5.0q is GIAB's
  README, now a row in section E.
- **Kamath, Şaşoğlu and Tse 2015**: the row in section H now says
  `abstract`, as the 30 September prose already did.

| Work | What it does | Bearing | Status |
|---|---|---|---|
| Berend and Sapir 2005, [Monotonicity in Condorcet jury theorem](https://doi.org/10.1007/s00355-003-0293-z), Soc. Choice Welf. 24(1):83-92; read as the author's preprint (`tzin.bgu.ac.il/~lsapir/Papers/CJT.pdf`, footer "10/08/2003") | Theorem 1: for voters with competence at least 1/2 and a random subcommittee of size m deciding by majority with a coin for ties, the probability of a correct decision is equal for m odd and m + 1, and non-decreasing from m even to m + 1 | The citation for f(2k) = f(2k - 1) in the vote model | `paper`, author's preprint: introduction, definitions, Theorem 1 and equation (4) read; the lemmas not checked |
| Andersson and Angner 2026, [A 250-year-old theorem walks into a courtroom](https://doi.org/10.1007/s00283-026-10556-9), Math. Intelligencer 48(3):270-273, CC BY 4.0 | a column on Condorcet's jury theorem; Remark 1 gives P_2m(p) = P_2m-1(p) under fair tie-breaking by a coupling argument, for homogeneous voters | A second, open statement of the identity, with the coupling proof the supplement uses | `paper` (full-text PDF read 1 Oct, `2026-10-01-followup.md` section 1; Remark 1 re-checked 1 Oct evening through the publisher page) |
| Dwarshuis et al. 2024, [The GIAB genomic stratifications resource for human reference genomes](https://doi.org/10.1038/s41467-024-53260-y), Nat. Commun. 15:9029, 14 authors; PMC11489684 | BED files that define genomic contexts (hard-to-map, GC-rich and others) on GRCh37, GRCh38 and CHM13, for stratifying benchmark performance | The citation for scoring by stratification in the Limits paragraph | `abstract` (Europe PMC, read 1 Oct); the Data availability passage was also read on 1 Oct (`2026-10-01-followup.md` section 2) |

## I. Follow-up: error-controlled phasing

Rows for the literature gate in `docs/followup/error-controlled-phasing.md`,
section 6, run on 2 October 2026. `docs/followup/literature-gate.md` has the
searches, what each source provides and the verdict. Bibliographic fields
were checked against Crossref, Europe PMC, the arXiv abstract page or the
publisher's proceedings page. Text was fetched as raw files (PDF, JATS XML,
BioC or Markdown) and read directly, not through a summarizing fetch tool.
A fourth status, `docs`, marks a tool's documentation or a specification
read at the source given in the row. Rows already in sections F and H are not
repeated; the gate cites ProbHap, PhaseQuality (Damaraju et al. 2026),
HiPhase, LongPhase, LongPhase 2 and pHapCompass from there.

### I1. Phase confidence and unphasing rules in phasers

| Work | What it does | Bearing | Status |
|---|---|---|---|
| Edge, Bafna, Bansal 2017, [HapCUT2](https://doi.org/10.1101/gr.213462.116), *Genome Res.* 27(5):801-812, PMC5411775 | Prunes a variant when the largest posterior of its four configurations is below 0.8 (default); for long reads, splits a block where the posterior of the adjacent configuration is below a threshold. On chromosome 1 of 11x PacBio it compares these scores with ProbHap's as trade-off curves: pruning 3.1 per cent of variants cut mismatch error 55.1 per cent (ProbHap under 15 per cent), and block splitting gave lower switch error at fixed AN50. On fosmid data it compares tools on the subset of variants that all tools phased | The nearest published risk-coverage comparison of phase confidence scores. The thresholds are fixed posteriors, not calibrated to an error rate; the main text reports no calibration measurement | `paper`, main text (post-processing results and methods); supplement, which holds the curves (Supplemental Fig. S4), not read |
| HapCUT2 [output format](https://github.com/vibansal/HapCUT2/blob/master/outputformat.md) and [README](https://github.com/vibansal/HapCUT2/blob/master/README.md), master branch, read 2 Oct 2026 | Block file per SNV: discrete pruning flag (RefHap rule: as many reads consistent as inconsistent), phred-scaled switch quality (0 to 100) and mismatch quality, and the count of haplotype-informative fragments. Switch quality is empty unless `--error_analysis_mode 1`; `--split_blocks 1` is documented as broken. Default prunes SNVs with mismatch quality below 6.98 (error probability 0.2); `--threshold` adjusts pruning | A per-boundary score exists in a shipped tool but is off by default and not written by default | `docs` |
| WhatsHap [user guide](https://github.com/whatshap/whatshap/blob/main/doc/guide.rst), `main`, read 2 Oct 2026 | "Variants can be left unphased for two reasons: Either the variant type is not supported or the phasing algorithm could not make a phasing decision." Phase sets are the connected components of the read-variant graph. On `PQ`: "WhatsHap currently does not add this tag." With `--distrust-genotypes` a variant switched to homozygous appears unphased. Polyphase cuts blocks by `--block-cut-sensitivity` | WhatsHap writes no per-variant or per-boundary quality; its unphasing rule is "no decision", not an error target | `docs` |
| HiPhase [user guide](https://github.com/PacificBiosciences/HiPhase/blob/main/docs/user_guide.md) and [methods](https://github.com/PacificBiosciences/HiPhase/blob/main/docs/methods.md), `main`, read 2 Oct 2026 | Writes `GT` and `PS`; a heterozygous call treated as homozygous in the optimal solution is left unphased with no `PS`; `PF=TR_OVERLAP` marks variants excluded because they lie inside a tandem-repeat call; `--min-vcf-qual` and `--min-mapq` filter inputs. The stats file reports per block `pruned_solutions` ("if 0, then this block has a guaranteed optimal solution") and `cost_ratio` | No per-variant phase quality. The only block-level signals are solver diagnostics | `docs` |
| LongPhase [README](https://github.com/twolinin/longphase/blob/main/README.md), `main`, v2.0.2, read 2 Oct 2026 | Output VCF carries `GT` and `PS` only. `--edgeThreshold` (default 0.7): "give up SNP-SNP phasing pair if the number of reads of the two combinations are similar"; `--baseQuality` and `--edgeWeight` down-weight low-quality edges | The unphasing rule is a fixed read-ratio threshold, with no stated error target and no output quality | `docs` |
| GATK 3.7-0 ReadBackedPhasing [tool documentation](https://github.com/broadinstitute/gatk-docs/blob/8fcf44bb0686f2f7d442aade181ff6ed508a97de/gatk3-tooldocs/3.7-0/org_broadinstitute_gatk_tools_walkers_phasing_ReadBackedPhasing.html), gatk-docs commit `8fcf44b` | Builds up to `2^n` local haplotypes over at most `--maxPhaseSites` (10) successive heterozygous SNPs and outputs phase only when the `PQ` score reaches `--phaseQualityThresh` (default 20.0); diploid SNPs only; phase recorded in the `HP` tag | A phred-scaled phase quality with a fixed threshold. The documentation does not define `PQ` further. The WhatsHap guide notes that RBP omits `PQ` on the first variant of a set, so it is relative to an earlier variant | `docs` |
| [VCF v4.3 specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf), section 1.6.2, hts-specs | "PQ (Integer): Phasing quality, the phred-scaled probability that alleles are ordered incorrectly in a heterozygote (against all other members in the phase set). We note that we have not yet included the specific measure for precisely defining 'phasing quality'; our intention for now is simply to reserve the PQ tag for future use" | The format reserves a phase-quality field and leaves its definition open | `docs` |
| Majidian and Sedlazeck 2020, [PhaseME](https://doi.org/10.1093/gigascience/giaa078), *GigaScience* 9(7):giaa078, PMC7379178 | Assesses and corrects any single-individual phasing VCF against population linkage from the 1000 Genomes Project; reports a per-block quality ratio (share of SNVs without conflict); splits blocks at conflicting SNVs; Hamming error down 22.4 per cent on average across five technologies | A tool-agnostic post-processor of phased VCFs, using population data rather than reads, with no error guarantee | `abstract`, plus the overview and quality-assessment paragraphs |
| Hofmeister, Ribeiro, Rubinacci, Delaneau 2023, [SHAPEIT5](https://doi.org/10.1038/s41588-023-01415-w), *Nat. Genet.* 55(7):1243-1249, PMC10335929 | Statistical phasing. A phasing confidence per rare heterozygous genotype, "the probability of the reported phase", from 0.5 to 1; filtering at 0.99 "controls the SER to a maximum of around 2% for WGS data and around 1% for WES data" while keeping over 75 and over 40 per cent of variants with MAC 2 to 5; SER is reported by threshold with the count filtered (Extended Data Fig. 5); compound heterozygous LoF calls are reported with and without the 0.99 filter | A published per-variant phase confidence with a risk-retention table, for statistical phasing. "Controls" is an empirical measurement against truth, not a finite-sample bound | `paper`, the phasing-confidence and compound-heterozygosity sections; methods not read |
| Lassen et al. 2024, [Exome-wide evidence of compound heterozygous effects](https://doi.org/10.1016/j.xgen.2024.100602), *Cell Genomics* 4(7):100602, PMC11293579 | Uses SHAPEIT5 posterior PP on UK Biobank exomes; measures SER against 96 trios and agreement with WhatsHap read-backed phasing of 453,514 variant pairs; reports "well-calibrated PP" and picks "the empirical cutoff of PP ≥0.9" for compound heterozygous calls | Read-based phase is used as the validation truth for statistical phase; the cutoff is chosen empirically | `paper`, the phasing-validation passages |
| Heinz, Marin, Meyerson, Li 2026, HipHap: Haplotype Assignment and Confidence Scoring for Diploid Reference Genomes, bioRxiv, posted 29 Sep 2026, DOI 10.64898/2026.09.24.754153 | Aligns each read to both haplotypes of a diploid assembly, assigns it to one, and reports a per-read assignment quality (HapQ) | A per-read confidence, not a per-relation one; it phases no variants and states no error bound. Found in the 10 October re-check (`2026-10-10-paper-recheck.md`), which also records what its code shows | `abstract` (read 10 Oct) |

### I2. Conformal prediction and risk control

| Work | What it does | Bearing | Status |
|---|---|---|---|
| Vovk, Gammerman, Shafer 2005, *Algorithmic Learning in a Random World*, Springer, doi 10.1007/b106715 | The conformal prediction monograph | Background citation only | `recalled`; Crossref fields verified, book not read |
| Papadopoulos, Proedrou, Vovk, Gammerman 2002, Inductive confidence machines for regression, ECML, LNCS 2430, 345-356, doi 10.1007/3-540-36755-1_29 | Inductive (split) conformal prediction | Origin of the split method | `recalled`; Crossref fields verified, not read |
| Lei, G'Sell, Rinaldo, Tibshirani, Wasserman 2018, [Distribution-free predictive inference for regression](https://doi.org/10.1080/01621459.2017.1307116), *JASA* 113(523):1094-1111; arXiv 1604.04173 | Full and split conformal inference with finite-sample marginal coverage for any regression estimator | The standard citation for split conformal | `abstract` |
| Angelopoulos and Bates 2023, [Conformal prediction: a gentle introduction](https://doi.org/10.1561/2200000101), *Found. Trends Mach. Learn.* 16(4):494-591; arXiv 2107.07511 | Tutorial; covers structured outputs, distribution shift, time series and models that abstain | A route into the extensions | `abstract` |
| Bates, Angelopoulos, Lei, Malik, Jordan 2021, [Distribution-free, risk-controlling prediction sets](https://doi.org/10.1145/3478535), *J. ACM* 68(6):1-34; arXiv 2101.02703 | Calibrates set-valued predictors on a holdout set so the expected loss is controlled, with finite-sample guarantees in high probability | Risk control in high probability, the alternative to control in expectation | `abstract` |
| Angelopoulos, Bates, Fisch, Lei, Schuster, [Conformal risk control](https://proceedings.iclr.cc/paper_files/paper/2024/hash/f3549ef9b5ff520a7e41ff3cc306ab2b-Abstract-Conference.html), ICLR 2024; arXiv 2208.02814 | Theorem 1: if each loss `L_i(λ)` is non-increasing and right-continuous in `λ`, bounded by `B`, and `L_i(λ_max) ≤ α`, the threshold `λ̂ = inf{λ : n/(n+1) R̂_n(λ) + B/(n+1) ≤ α}` gives `E[L_{n+1}(λ̂)] ≤ α`; the bound is tight up to `2B/(n+1)`. Proposition 1: for non-monotone losses the procedure can miss by an arbitrary amount | The procedure named in section 3.2 of the plan. A per-window ratio "errors among trusted over trusted" is not monotone in the threshold, so Theorem 1 does not apply to it as stated | `paper`, sections 1 and 2 |
| Angelopoulos, Bates, Candès, Jordan, Lei 2025, [Learn then Test](https://doi.org/10.1214/24-AOAS1998), *Ann. Appl. Stat.* 19(2); arXiv 2110.01052 | Risk control as multiple hypothesis testing, with finite-sample high-probability guarantees for any model, monotone or not. Section "Selective classification with accuracy control": predict when the top score exceeds `λ`, control "the error conditionally on predicting" with exact binomial tail p-values; calibration data i.i.d. | The published form of "error among trusted at most alpha", for i.i.d. points | `paper`, abstract and the selective classification section |
| Barber, Candès, Ramdas, Tibshirani 2023, [Conformal prediction beyond exchangeability](https://doi.org/10.1214/23-AOS2276), *Ann. Statist.* 51(2); arXiv 2202.13415 | Weighted quantiles and a randomization so coverage degrades gracefully when exchangeability fails | One tool for genomic dependence and drift | `abstract` |
| Chernozhukov, Wüthrich, Zhu 2018, [Exact and robust conformal inference methods for predictive machine learning with dependent data](https://arxiv.org/abs/1802.06300), COLT, PMLR 75:732-749 | Block structures in the permutation scheme; exact under exchangeability, approximately valid under serial dependence | Blocks of adjacent pairs are the natural analogue along a chromosome | `abstract` |
| Oliveira, Orenstein, Ramos, Romano 2024, [Split conformal prediction and non-exchangeable data](https://arxiv.org/abs/2203.15885), *JMLR* 25(225) | Split conformal stays valid for many non-exchangeable processes with a small coverage penalty, by concentration and decoupling | The penalty a mixing assumption along the genome would cost | `abstract` |
| Mao, Martin, Reich 2024, [Valid model-free spatial prediction](https://doi.org/10.1080/01621459.2022.2147531), *JASA* 119(546):904-914 | Local spatial conformal prediction; spatial data locally approximately exchangeable under infill asymptotics | Spatial conformal; the genome is a one-dimensional index | `abstract` |
| Dunn, Wasserman, Ramdas, [Distribution-free prediction sets for two-layer hierarchical models](https://doi.org/10.1080/01621459.2022.2060112), *JASA* 118(544):2491-2502, online 9 May 2022; arXiv 1809.07441 | Conformal sets when observations come in groups from distinct distributions; repeated subsampling gives coverage guarantees | Windows or genomes as groups and pairs as observations is this two-layer setting | `abstract` |
| Jin and Candès 2023, [Selection by prediction with conformal p-values](https://www.jmlr.org/papers/v24/22-1176.html), *JMLR* 24(244):1-41; arXiv 2210.01408 | Selects test units whose unobserved outcome passes a threshold while controlling the proportion of falsely selected units (FDR), from conformal p-values and Benjamini-Hochberg; holds under exchangeability conditions | Trusting pairs is selecting units; the error among trusted pairs is a false discovery proportion | `abstract` |
| Vovk 2012, [Conditional validity of inductive conformal predictors](https://arxiv.org/abs/1209.2673), ACML 2012 (per arXiv); arXiv 1209.2673v2 | Section 3, Propositions 2a and 2b: inductive conformal predictors are training-conditional, "of the PAC type", and Proposition 2b's condition is a binomial tail, `delta >= bin_{n,E}(floor(epsilon (n+1)) - 1)`. Section 4: the conditional ICP takes an inductive taxonomy `K` and computes p-values within the test example's category; Proposition 3: under exchangeability the error probability given the category is at most `epsilon`. Section 7 notes the training-conditional version holds for label-conditional ICPs | Track A's per-stratum Learn then Test is this construction (calibrate within each category of a fixed taxonomy, with a binomial PAC bound) applied to Learn then Test's risk. The procedure is not new | `paper`, sections 1, 3, 4 and 7, from the arXiv PDF |
| Gibbs, Cherian, Candès, [Conformal prediction with conditional guarantees](https://arxiv.org/abs/2305.12616), arXiv 2305.12616 (v3, 16 Sep 2024) | Reformulates conditional coverage as coverage over a class of covariate shifts; for a finite collection of subgroups, exact finite-sample coverage over each group | The route if strata overlapped; track A's strata are disjoint by construction | `abstract`; journal version not checked |
| Tibshirani, Barber, Candès, Ramdas, [Conformal prediction under covariate shift](https://arxiv.org/abs/1904.06019), arXiv 1904.06019 (v3, 6 Jul 2020) | Weighted conformal prediction when the test and training covariate distributions differ and their likelihood ratio is known or estimable from unlabelled test points | A test genome with a different mix of strata is a shift in the stratum label only; calibrating within each stratum makes each stratum's guarantee independent of the mix, so no reweighting is needed for per-stratum statements | `abstract`; venue not checked |
| Blot, Angelopoulos, Jordan, Brunel, [Automatically adaptive conformal risk control](https://arxiv.org/abs/2406.17819), arXiv 2406.17819 (v2, 27 Mar 2025) | Approximate conditional control of expected risk, building on Gibbs et al., with a data-driven conditioning class | Risk in expectation and data-chosen groups; track A fixes its groups in advance and controls risk in high probability | `abstract` |
| Browning and Browning 2022, [Genotype error biases trio-based estimates of haplotype phase accuracy](https://doi.org/10.1016/j.ajhg.2022.04.019), *AJHG* 109:1016-1025; bioRxiv 10.1101/2022.04.06.487354 | Estimates the genotype error rate from parent-offspring trios and the bias it puts into trio-based switch error; in UK Biobank White British trios the observed switch error rate is 2.4 times the estimated true rate (1.4e-3 against 5.8e-4) | Label noise from transmission is the first question for Y1 (the October directions memo (development repository)); its gate measures it on the GIAB Ashkenazi trio | `abstract` (Europe PMC) |
| Miller and Piccolo 2021, [trioPhaser](https://doi.org/10.1186/s12859-021-04470-4), *BMC Bioinformatics* 22:559 | Combines Mendelian inheritance logic with SHAPEIT4 to phase trios; reports concordance with WhatsHap, not confidence or error bounds | Family phasing as a phaser, not as calibration labels | `abstract` |

### I3. Selective classification

| Work | What it does | Bearing | Status |
|---|---|---|---|
| El-Yaniv and Wiener 2010, [On the foundations of noise-free selective classification](https://www.jmlr.org/papers/v11/el-yaniv10a.html), *JMLR* 11(53):1605-1641 | Names the risk-coverage (RC) trade-off of classification with a reject option and characterizes it for noise-free models | The origin of the risk-coverage curve the benchmark would use | `abstract` |
| Geifman and El-Yaniv 2017, [Selective classification for deep neural networks](https://papers.nips.cc/paper_files/paper/2017/hash/4a8423d5e91fda00bb7e46540e2b0cf1-Abstract.html), NeurIPS 30; arXiv 1705.08500 | Selective risk `R(f, g) = E[ℓ g] / E[g]`; Selection with Guaranteed Risk (SGR) finds a threshold on a confidence score so that `Pr{R > r*} < δ` over an i.i.d. sample (Theorem 3.2, binomial tail bound with a union over the binary search), maximizing coverage | Exactly the guarantee in section 2 of the plan, for i.i.d. data and a fixed classifier | `paper`, sections 1 to 3 |
| Hallberg Szabadváry, Löfström, Johansson, Sönströd, Ahlberg, Carlsson 2025, [Classification with reject option: distribution-free error guarantees via conformal prediction](https://doi.org/10.1016/j.mlwa.2025.100664), *Mach. Learn. Appl.* 20:100664; arXiv 2506.21802 | Accepts only singleton conformal sets in binary classification and proves the resulting error rate; error-reject curves | A conformal route to a reject option with an error bound | `abstract` |

### I4. Conformal methods in genomics

| Work | What it does | Bearing | Status |
|---|---|---|---|
| Elrick et al. 2025, [SAVANA](https://doi.org/10.1038/s41592-025-02708-0), *Nat. Methods* 22(7):1436-1446, PMC12240814 | Somatic SV and copy-number calling from long reads; optional Mondrian conformal prediction on a random-forest score, with a 70/30 split; states that MCP bounds the per-class error "provided that the training and test examples are independently and identically distributed" and that other SV callers' quality scores "are not formally guaranteed to be well calibrated" | The nearest use of conformal prediction in long-read variant calling that we found. Not phasing; the guarantee is for i.i.d. calls | `paper`, the MCP results paragraph and methods |
| Papangelou, Kyriakidis, Natsiavas, Chouvarda, Malousi 2025, [Reliable machine learning models in genomic medicine using conformal prediction](https://doi.org/10.3389/fbinf.2025.1507448), *Front. Bioinform.* 5:1507448 | Perspective review of conformal prediction in genomic medicine, with drug-response and subtyping examples | No sequencing-level application | `abstract` |

### I6. Truth sets for the test

| Work | What it does | Bearing | Status |
|---|---|---|---|
| Eberle et al. 2017, [A reference data set of 5.4 million phased human variants validated by genetic inheritance from sequencing a three-generation 17-member pedigree](https://doi.org/10.1101/gr.210500.116), *Genome Res.* 27(1):157-164, PMID 27903644 | Platinum Genomes: deep Illumina sequencing of a 17-member, three-generation pedigree; a phased catalog of 4.7 million SNVs and 0.7 million indels consistent with inheritance | The HG001 (NA12878) truth of `docs/design/confidence-run-plan.md` | `abstract` (Europe PMC, read 2 Oct) |
| Kronenberg et al. 2025, [The Platinum Pedigree: a long-read benchmark for genetic variants](https://doi.org/10.1038/s41592-025-02750-y), *Nat. Methods* 22:1669-1676, PMID 40759746 | Variants in CEPH-1463 filtered by Mendelian inheritance across PacBio HiFi, Illumina and ONT: 4.7 million SNVs, 767,795 indels, 537,486 tandem repeats and 24,315 SVs over 2.77 Gb of GRCh38 | A candidate secondary truth for HG001 that reaches difficult regions; not in the registration | `abstract` (Europe PMC, read 2 Oct) |

### I5. Benchmarking guidance and clinical use

| Work | What it does | Bearing | Status |
|---|---|---|---|
| Krusche et al. 2019 (GA4GH Benchmarking Team, with GIAB), [Best practices for benchmarking germline small-variant calls in human genomes](https://doi.org/10.1038/s41587-019-0054-x), *Nat. Biotechnol.* 37:555-560, PMC6699627 | "Methods for assessing accuracy of phasing have not been fully developed or standardized"; "assessing phasing requires us to consider not only phasing variant accuracy, but also completeness of phasing coverage"; no comparison method supports phase sets | The GIAB-associated statement that phasing accuracy and completeness must be read together. It proposes no matched-completeness metric | `paper`, the phasing passages (full text through the PMC BioC service) |
| Richards et al. 2015, [ACMG/AMP standards and guidelines for the interpretation of sequence variants](https://doi.org/10.1038/gim.2015.30), *Genet. Med.* 17(5):405-424, PMC4544753 | PM3: for a recessive disorder, a variant found in trans with a pathogenic variant is moderate evidence; phase is determined by "testing parental samples" | The clinical criterion that a trusted in-trans call feeds. No error rate is set for the phase call | `paper`, the in cis / in trans passages |
| Del Gobbo and Boycott 2025, [The additional diagnostic yield of long-read sequencing in undiagnosed rare diseases](https://doi.org/10.1101/gr.279970.124), *Genome Res.* 35(4):559-571, PMC12047273 | Section "Improved read-based phasing": "Laboratories primarily rely on familial genotyping or direct read-based phasing if possible"; cites long-read cases resolving compound heterozygosity (CC2D2A, Werner syndrome) | Read-based phase is in clinical use for in-trans calls, reported by block N50 and genes covered, not by error rate | `paper`, the read-based phasing section |
| ClinGen Sequence Variant Interpretation, recommendation for the in trans criterion PM3 | Points for PM3 by whether phase is confirmed | Would say how clinical curation weights unconfirmed phase | `recalled`; search snippets only, the document was not obtained |

### I7. Round 6 search, 17 October 2026 (cited in the combined paper's Part III)

Bibliographic fields checked against Crossref on 18 October 2026. Status is
what the 17 October search read, as recorded in
the October directions memo (development repository) (round 6, prior work).

| Work | What it does | Bearing | Status |
|---|---|---|---|
| Kong et al. 2008, [Detection of sharing by descent, long-range phasing and haplotype imputation](https://doi.org/10.1038/ng.216), *Nat. Genet.* 40(9):1068-1075 | Uses relatives sharing segments by descent as surrogate parents to phase SNP-array data | Relatives as a phase source predate Q1; Q1's use is calibration of a read-based certificate | `abstract` |
| O'Connell et al. 2016, [Haplotype estimation for biobank-scale data sets](https://doi.org/10.1038/ng.3583), *Nat. Genet.* 48(7):817-820 | SHAPEIT3; scores statistical phasing against IBD1 segments of 384 sibling pairs | Siblings as phasing truth are not new | `paper` |
| Chen et al. 2025, [Genome-wide profiling of highly similar paralogous genes using HiFi sequencing](https://doi.org/10.1038/s41467-025-57505-2), *Nat. Commun.* 16:2340 | Paraphase: phases all copies of 160 segmental-duplication regions | A dedicated caller for paralogs, not a flag for a generic phaser's errors (Q2) | `abstract` |
| Jamshidi et al. 2025, [Long-range PCR and Nanopore sequencing for localisation and phasing variants](https://doi.org/10.1186/s12920-025-02251-z), *BMC Med. Genomics* 18:186 | Clinical cis or trans workflow; 15 of 15 HG002 pairs concordant | Concordance on few pairs, no bounded error rate (Q4) | `paper` |
| Devaney et al. 2026, [Sensitivity of HiFi long-read genome sequencing for difficult-to-detect pathogenic variants](https://doi.org/10.1016/j.ajhg.2026.04.001), *Am. J. Hum. Genet.* 113(5):1036-1048 | Clinical HiFi pipeline with WhatsHap, homologous genes sent to Paraphase | No phase error rate reported (Q4) | `paper` |
| Akbari et al. 2022, [Genome-wide detection of imprinted differentially methylated regions using nanopore sequencing](https://doi.org/10.7554/eLife.77898), *eLife* 11:e77898 | The imprinted DMR list W2 uses (supplementary file 1) | Data source only | `abstract` |
