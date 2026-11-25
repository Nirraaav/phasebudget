# Bibliography verification and novelty sweep, 30 September 2026

This note has two parts. Part A checks the bibliography entries in
the earlier read-selection paper (development repository) that are marked unverified and gives a corrected entry
for each. Part B reports a sweep of 2023 to 2026 work on five questions that
bear on the paper's contributions. Nothing in this note edits the paper or the
prior-art matrix; the corrections are for the owner to apply.

## How status is assigned here

Statuses follow `01-prior-art-matrix.md`, with the stricter reading of
section H.

- Bibliographic fields (authors, venue, volume, pages, DOI) marked "verified"
  come from the publisher's Crossref record, fetched from
  `https://api.crossref.org/works/<DOI>`, or from PubMed and PMC metadata.
  Verifying a citation does not change the status of the row that describes
  the paper's content.
- `paper, searched` means the full text was retrieved as raw JAT XML from the
  Europe PMC REST service (`/rest/<PMCID>/fullTextXML`), and the passages
  quoted or summarised below were read directly. The rest of the paper was
  not read. This follows the convention of the LongPhase row in section H.
- `abstract` means the abstract, or a publisher landing page, was read.
- `recalled` means the item is known only from a search-result snippet.

## A. Bibliography verification

### Summary of corrections

| Key | What was wrong or missing in `paper.tex` | Source checked |
|---|---|---|
| `pirola2016` | pages and DOI missing | Crossref, DOI 10.1093/bioinformatics/btv495 |
| `kuleshov2014` | pages and DOI missing; the pages carry an `i` prefix because the paper is in the ECCB 2014 proceedings issue | Crossref, DOI 10.1093/bioinformatics/btu484 |
| `makinen2015` | volume and pages were correct; DOI missing; the key says 2015 and the journal year is 2017 | Crossref, DOI 10.1016/j.dam.2016.08.021 |
| `damaraju2026` | author list missing; there is no journal. PMC13317610 and PubMed 42388792 index the bioRxiv preprint itself | PubMed esummary for 42388792; PMC XML front matter; Crossref for the bioRxiv DOI; Europe PMC |
| `q100` | the authors are 65 named people, first author Nancy F. Hansen; volume, issue, pages and DOI missing | Crossref, DOI 10.1016/j.cell.2026.06.016, and the preprint DOI 10.1101/2025.09.21.677443 |
| `kamath2015` | ISIT pages and IEEE DOI missing | Crossref, DOI 10.1109/isit.2015.7282588; arXiv abstract page 1502.01975 |
| `chen2016` | correct; editors, arXiv id and PMLR URL can be added | PMLR proceedings page `proceedings.mlr.press/v48/chena16.html`; arXiv API |

### Details

**pirola2016.** Crossref: Yuri Pirola, Simone Zaccaria, Riccardo Dondi,
Gunnar W. Klau, Nadia Pisanti, Paola Bonizzoni. "HapCol: accurate and
memory-efficient haplotype assembly from long reads." *Bioinformatics*
32(11):1610-1617, issue date 1 June 2016. DOI 10.1093/bioinformatics/btv495.

**kuleshov2014.** Crossref: Volodymyr Kuleshov. "Probabilistic
single-individual haplotyping." *Bioinformatics* 30(17):i379-i385, 2014
(online 22 August 2014). DOI 10.1093/bioinformatics/btu484. The `i` pages are
the ECCB 2014 proceedings supplement of that issue; the author's copy is at
`cs.cornell.edu/~kuleshov/papers/eccb2014.pdf`, which is consistent.

**makinen2015.** Crossref: Veli Mäkinen, Valeria Staneva, Alexandru I.
Tomescu, Daniel Valenzuela, Sebastian Wilzbach. "Interval scheduling
maximizing minimum coverage." *Discrete Applied Mathematics* 225:130-135, July
2017. DOI 10.1016/j.dam.2016.08.021. The volume and pages already in
`paper.tex` are correct and can lose the "unverified" note. The arXiv
preprint is 1508.07820. Renaming the key to `makinen2017` is optional; the
printed year should be 2017.

**damaraju2026.** PubMed esummary for PMID 42388792 gives source "bioRxiv",
publication type "Journal Article; Preprint", full journal name "bioRxiv : the
preprint server for biology", and the article ids PMC13317610 and DOI
10.64898/2026.05.04.722832. The PMC XML front matter has
`<journal-title>bioRxiv</journal-title>` and
`article-version-type="status"` equal to `preprint`, PMC version 2, epub 10
May 2026. Crossref gives posted date 8 May 2026 for the bioRxiv DOI. A title
search of Europe PMC returns only the preprint record (source `PPR`). No
journal version was found. A web search turned up a University of
Washington doctoral dissertation by Nikhita Damaraju (2026, under embargo to
August 2027) that covers the same topic; that item is `recalled`.

Authors, as given by PubMed and Crossref: Nikhita Damaraju, F. Graeme Frost,
Jiayu Fu, David A. D'Onofrio, Joy Goffena, Sophie H. Storz, Zach B.
Anderson, Trent Prall, Miranda P. Zalusky, May Christine V. Malicdan, David
R. Adams, Danny E. Miller.

This corrects section H of the prior-art matrix, which describes Damaraju et
al. as "now peer reviewed". The PMC record is a preprint, and the paper
should cite it as a bioRxiv preprint.

**q100.** Crossref for DOI 10.1016/j.cell.2026.06.016: "A complete diploid
human genome benchmark for personalized genomics." *Cell* 189(16):4857-4875.e31,
August 2026. Publisher item identifier S0092867426007038, which matches the
Cell URL in section H. Crossref lists 65 authors and no group author. The
bioRxiv preprint, DOI 10.1101/2025.09.21.677443, posted 21 September 2025,
has the same 65 authors in the same order.

Authors in Crossref order: Nancy F. Hansen, Nathan Dwarshuis, Hyun Joo Ji,
Arang Rhie, Hailey Loucks, Glennis A. Logsdon, Mitchell R. Vollger, Jessica
M. Storer, Juhyun Kim, Eleni Adam, Nicolas Altemose, Dmitry Antipov, Mobin
Asri, Sofia Barreira, Stephanie C. Bohaczuk, Andrey V. Bzikadze, Sara A.
Carioscia, Andrew Carroll, Kuan-Hao Chao, Yanan Chu, Arun Das, Peter Ebert,
Adam English, Mark Fleharty, Laura E. Fleming, Giulio Formenti, Andrea
Guarracino, Gabrielle A. Hartley, Katharine Jenike, Jenna Kalleberg, Yu Kang,
Robert King, Josipa Lipovac, Mira Mastoras, Matthew W. Mitchell, Shloka Negi,
Nathan D. Olson, Keisuke K. Oshima, Luis F. Paulin, Brandon D. Pickett, David
Porubsky, Jane Ranchalis, Desh Ranjan, Mikko Rautiainen, Harold Riethman,
Robert D. Schnabel, Fritz J. Sedlazeck, Kishwar Shafin, Mile Sikic, Steven J.
Solar, Alexander P. Sweeten, Winston Timp, Justin Wagner, DongAhn Yoo, Ying
Zhou, Erik Garrison, Evan E. Eichler, Michael C. Schatz, Andrew B.
Stergachis, Rachel J. O'Neill, Karen H. Miga, Steven L. Salzberg, Sergey
Koren, Justin M. Zook, Adam M. Phillippy.

The bibliographic entry is verified. The content of the paper has not been
read, so any sentence in `paper.tex` about what the benchmark contains still
rests on the `recalled` row in section H. Damaraju et al. (methods, read in
the PMC XML) state that they used "the VCF from the Telomere-to-Telomere
Consortium's Q100 project (https://github.com/marbl/HG002) created using
diploid assemblies" as the HG002 truth for T2T alignments, which is a second
source for the existence of that truth set.

**kamath2015.** Crossref for DOI 10.1109/isit.2015.7282588: Govinda M.
Kamath, Eren Sasoglu (Şaşoğlu on arXiv and bioRxiv), David Tse. "Optimal
haplotype assembly from high-throughput mate-pair reads." In *2015 IEEE
International Symposium on Information Theory (ISIT)*, Hong Kong, 14-19 June
2015, pages 914-918. The arXiv version is 1502.01975v1, 6 February 2015, with
the comment "Submitted to ISIT 2015" and no journal reference. A bioRxiv copy
exists at DOI 10.1101/014993. The arXiv abstract was read, so the content
status rises from `recalled` to `abstract`.

**chen2016.** PMLR page: Yuxin Chen, Govinda Kamath, Changho Suh, David Tse.
"Community recovery in graphs with locality." In *Proceedings of The 33rd
International Conference on Machine Learning*, edited by Maria Florina Balcan
and Kilian Q. Weinberger, PMLR volume 48, pages 689-698, New York, 20-22 June
2016. URL `https://proceedings.mlr.press/v48/chena16.html`. arXiv 1602.03828.
PMLR proceedings have no DOI.

### Replacement entries for `paper.tex`

```latex
\bibitem{pirola2016} Y.~Pirola, S.~Zaccaria, R.~Dondi, G.~W.~Klau, N.~Pisanti, P.~Bonizzoni. HapCol: accurate and memory-efficient haplotype assembly from long reads. \emph{Bioinformatics} 32(11):1610--1617, 2016. \url{https://doi.org/10.1093/bioinformatics/btv495}
\bibitem{kuleshov2014} V.~Kuleshov. Probabilistic single-individual haplotyping. \emph{Bioinformatics} 30(17):i379--i385, 2014. \url{https://doi.org/10.1093/bioinformatics/btu484}
\bibitem{makinen2015} V.~M\"akinen, V.~Staneva, A.~I.~Tomescu, D.~Valenzuela, S.~Wilzbach. Interval scheduling maximizing minimum coverage. \emph{Discrete Appl. Math.} 225:130--135, 2017. \url{https://doi.org/10.1016/j.dam.2016.08.021}
\bibitem{damaraju2026} N.~Damaraju, F.~G.~Frost, J.~Fu, D.~A.~D'Onofrio, J.~Goffena, S.~H.~Storz, Z.~B.~Anderson, T.~Prall, M.~P.~Zalusky, M.~C.~V.~Malicdan, D.~R.~Adams, D.~E.~Miller. Determinants of haplotype phasing accuracy in long-read human genome sequencing. \emph{bioRxiv} preprint, 2026. \url{https://doi.org/10.64898/2026.05.04.722832}
\bibitem{kamath2015} G.~M.~Kamath, E.~\c{S}a\c{s}o\u{g}lu, D.~Tse. Optimal haplotype assembly from high-throughput mate-pair reads. In \emph{Proc. IEEE Int. Symp. Inf. Theory (ISIT)}, 914--918, 2015. \url{https://doi.org/10.1109/ISIT.2015.7282588}
\bibitem{chen2016} Y.~Chen, G.~Kamath, C.~Suh, D.~Tse. Community recovery in graphs with locality. In \emph{Proc. 33rd Int. Conf. Mach. Learn. (ICML)}, PMLR 48:689--698, 2016. \url{https://proceedings.mlr.press/v48/chena16.html}
\bibitem{q100} N.~F.~Hansen, N.~Dwarshuis, H.~J.~Ji, A.~Rhie, et al. (65 authors), S.~Koren, J.~M.~Zook, A.~M.~Phillippy. A complete diploid human genome benchmark for personalized genomics. \emph{Cell} 189(16):4857--4875.e31, 2026. \url{https://doi.org/10.1016/j.cell.2026.06.016}
```

If the paper's style lists every author, the full Q100 list above can replace
"et al.". The Damaraju entry has twelve authors, so it is written out.

## B. Novelty sweep, 2023 to 2026

### What was searched

- Crossref bibliographic queries for each title in part A.
- Europe PMC REST, restricted to publication years 2023 to 2026:
  read selection or coverage cap with WhatsHap; read selection with haplotype
  assembly or phasing; exact or optimal MEC; phasing with calibration and
  switch error or phasing confidence; switch error with spanning or
  informative reads; pedigree or trio read-based phasing with coverage;
  preregistration with benchmarking in bioinformatics or genomics;
  "Registered Report" with benchmark and sequencing.
- arXiv API: phasing with coverage or read selection; haplotype assembly with
  information-theoretic or coverage terms; minimum error correction.
- Web search: WhatsHap read selection against optimal selection; switch error
  as a function of spanning reads; coverage theory of phasing error;
  preregistered benchmarking; per-pair phase confidence calibration; pedigree
  read-selection budgets; majority vote with an even number of voters and
  random tie-breaking.
- Primary sources: the WhatsHap `CHANGES.rst` on the `main` branch and the
  `guide.rst` user guide, both fetched 30 September 2026.

Section H of the prior-art matrix already covers HiPhase, LongPhase and
LongPhase 2, LongHap, snpphaser, HapBridge, QHap, pHapCompass and the
information-theory rows. Those are not repeated below unless something new
was found.

### B1. WhatsHap's read selection or cap against exact or optimal selection

| Work | What it does | Bearing | Status |
|---|---|---|---|
| Battistella, Maheshwari, Ekim, Berger, Popic 2025, "Graph-based deep reinforcement learning for haplotype assembly with Ralphi", *Genome Res.* 35(12):2617-2625, DOI 10.1101/gr.280569.125, PMC12667386 | Learns a read-to-haplotype partition on fragment graphs with a graph network and actor-critic RL, rewarded by the maximum fragment cut. Downsamples HG002 ONT R10, HiFi and Illumina to 5, 10, 15 and 30x. The text says Ralphi phases "the full read set, whereas WhatsHap downsamples input data sets to 15x by default", and that Ralphi optionally uses "WhatsHap's downsampling routine"; the acknowledgements say the authors reuse WhatsHap's read selection module for preprocessing | Contribution 2: a 2025 Genome Research paper runs WhatsHap at its default cap on HG002 and reuses its selector, without comparing the selector to any optimum or varying the cap. It is a citation for "the shipped selector is reused downstream", and it should be cited among the HG002 long-read comparisons | `paper, searched` |
| WhatsHap `CHANGES.rst`, main branch | v2.8 is dated 8 June 2025. The only algorithmic change near selection since 2023 is the `--algorithm=heuristic` diploid phaser, which the changelog says gives "support for higher coverages and/or larger pedigrees at the cost of not solving the underlying MEC model to optimality anymore" and which it does not recommend "especially for pedigree phasing". No entry after "Smart read selection" (an old release) changes the selector | Confirms that the selector certified in the paper is the one shipped in 2.8. The version pin in the paper can cite this file | source read |
| WhatsHap user guide, `guide.rst` | Describes the default as "read selection: a heuristic that extracts a subset of the reads that is most informative for phasing", with optional read merging. It states that running time "grows quadratically with the coverage" for polyploid phasing | The documentation calls the selector a heuristic, which the paper's framing can quote | source read |
| Hartmann, Shaham, Nathan, Blech, Zeevi 2025, "TinkerHap: a novel read-based phasing algorithm with integrated multimethod support for enhanced accuracy", *GigaScience* 14:giaf138, DOI 10.1093/gigascience/giaf138, PMC12723663 | Read-based phaser that extends haplotypes from overlapping reads and merges statistical or pedigree phasing; evaluated on UK Biobank trios and the GIAB Ashkenazi trio | No read selection and no cap. Relevant to B2 only | `paper, searched` |

No item found compares WhatsHap's selector or its cap with an exact or
certified optimum. That is a report of these searches, and the paper should
not state it as a property of the literature.

### B2. Models of switch error from reads spanning adjacent variants, and when the cap matters

| Work | What it does | Bearing | Status |
|---|---|---|---|
| Damaraju et al. 2026, bioRxiv (part A) | Downsampling by `samtools --subsample` with a fixed seed to 10, 20 and 30x; "phasing error rates plateau at ~20x coverage"; at 10x on GRCh38, error beyond 30 to 50 kb reaches 20 to 30 per cent for ONT HG002. WhatsHap is run per chromosome with `--ignore-read-groups`; no `--max-coverage` value is stated, so the default applies | Empirical coverage curves by distance, with no model of error as a function of informative-read count and no variation of the selector or cap. The plateau near 20x sits above WhatsHap's default cap of 15x, which the paper can note, with the caveat that their depth is sequencing depth and the cap acts on phase-informative reads | `paper, searched` |
| Hartmann et al. 2025, TinkerHap (B1) | Methods: a read's phase "is determined by propagating the majority consensus" of overlapping reads | An algorithm that decides by majority among overlapping reads. It is not a model of the error rate, but a reader may ask how the paper's majority-vote model relates to it | `paper, searched` |
| Si, Vikalo, Vishwanath 2014, "Haplotype assembly: an information theoretic view", *2014 IEEE Information Theory Workshop (ITW)*, 182-186, DOI 10.1109/ITW.2014.6970817; arXiv 1404.0097; journal version "Information-theoretic analysis of haplotype assembly", *IEEE Trans. Inf. Theory* 63(6):3468-3479, 2017, DOI 10.1109/TIT.2017.2686884 | Recasts haplotype recovery as joint source-channel coding over a channel with erasures and errors; derives the number of reads needed for reliable recovery | Corrects section H: the venue is ITW 2014, not ISIT 2014, although the author's file is named `isit14.pdf`; the authors are now known. Prior art for any statement about read counts needed for correct phase. Outside 2023 to 2026, but needed to state the majority-vote model against it | `abstract` (arXiv); bibliographic fields verified in Crossref |
| Chen, Suh, Goldsmith 2016, "Information recovery from pairwise measurements", *IEEE Trans. Inf. Theory* 62(10):5881-5905, DOI 10.1109/TIT.2016.2600566; arXiv 1504.01369; ISIT 2014 and 2015 versions | Fundamental recovery criteria for node variables observed through noisy pairwise differences on a graph, for general channels | A boundary vote is a noisy pairwise difference. This is the general theory behind Chen et al. 2016 (ICML), already cited | `abstract` |
| Das and Vikalo 2014, "Optimal haplotype assembly with statistical pruning", *IEEE GlobalSIP*, 1330-1333, DOI 10.1109/GlobalSIP.2014.7032339 | Title only | Possibly relevant to pruning under exact solving; not read | `recalled` (Crossref title) |
| Ben-Yashar 2023, "An application of simple majority rule to a group with an even number of voters", *Theory and Decision* 94(1):83-95, DOI 10.1007/s11238-022-09872-1 | The abstract: the odd-group assumption of Condorcet's jury theorem "is not necessary"; simple majority rule is applied to an even-numbered homogeneous committee | The paper's even-count proposition is a statement about a homogeneous even committee with random tie-breaking. This is the nearest published treatment found. Whether it contains the identity P(correct, 2k) = P(correct, 2k-1) is not known from the abstract; the full text must be read before the paper cites it for that identity | `abstract` |
| Wikipedia, "Condorcet's jury theorem" | States that for even n the same argument works "if ties are broken by adding a single voter" or by fair coin-flips | Secondary; indicates the identity is known in the jury-theorem literature. Not citable | not a source |
| Shteingart et al. 2020, "Majority voting and the Condorcet's jury theorem", arXiv 2002.03153 | A proof of the jury theorem for odd n | Full text searched for the even case; it restricts to odd n and does not state the identity. Outside the date window | `paper` |
| Bar-Ilan University seminar, M. Mann, 13 June 2022, "The mathematical relationship between even and odd homogeneous groups" | Seminar announcement describing theorems relating even and odd committees | Suggests a paper may exist; none found | `recalled` |

Recommendation for the paper: present the even-count proposition as a short
derivation of a known property of majority voting with random tie-breaking,
and cite the jury-theorem literature for it once a primary statement is read.
Ben-Yashar 2023 is the first candidate to read in full. The contribution is
applying it to boundary votes and measuring it on HG002, not the identity.

### B3. Pre-registered benchmarking in bioinformatics methods evaluation

| Work | What it does | Bearing | Status |
|---|---|---|---|
| Teo, Squair, Courtine, Skinnider 2024, "Best practices for differential accessibility analysis in single-cell epigenomics", *Nat. Commun.* 15:8805, DOI 10.1038/s41467-024-53089-5, PMC11470024 | A benchmark of scATAC-seq differential-accessibility methods published as a Registered Report; "The Stage 1 protocol for this Registered Report was accepted in principle on 25th October 2023", protocol at figshare 10.6084/m9.figshare.24541816.v1. Analyses requested in review are reported in a separate section as "post hoc, unregistered analyses" | Direct precedent for a pre-registered benchmark of bioinformatics methods, with the same split between registered and exploratory analyses that the paper uses. Should be cited where the paper describes its pre-registration | `paper, searched` |
| Mallona, Soneson, Carrillo, et al. 2025, "Building a continuous benchmarking ecosystem in bioinformatics", *PLoS Comput. Biol.* 21(11):e1013658, DOI 10.1371/journal.pcbi.1013658, PMC12614554 | Perspective. Says method developers who benchmark their own tools can use "pre-registration and posting the parameters of the study in advance at a registry", cites Teo et al. as an example, and notes "a tension between preregistration and the concept of 'continuous' benchmarking" | Supports the rationale for pre-registration when the authors built one of the methods compared, which is this paper's situation | `paper, searched` |
| Lange, Wilcke, Hoffmann, Herrmann, Boulesteix 2025, "On 'confirmatory' methodological research in statistics and related fields", *Stat. Med.* 44(25-27):e70303, DOI 10.1002/sim.70303, PMC12600059 | Reviews the exploratory and confirmatory distinction, HARKing, and preregistration for methodological research | The methodological-research citation for pre-registration | `paper, searched` |
| Kohnert and Kreutz 2025, "Computational study protocol: leveraging synthetic data to validate a benchmark study for differential abundance tests for 16S microbiome sequencing data", *F1000Research* 13:1180 (version 2), DOI 10.12688/f1000research.155230.2 | A published study protocol for a computational benchmark, following SPIRIT; its history table records submission as a registered report to PLOS Biology and PLOS ONE, not accepted | A second example of a protocol-first computational benchmark | `abstract` |
| Siepe, Bartoš, Morris, Boulesteix, Heck, Pawel 2024, "Simulation studies for methodological research in psychology: a standardized template for planning, preregistration, and reporting", *Psychol. Methods*, DOI 10.1037/met0000695 | The ADEMP-PreReg template | A template the pre-registration document could be checked against | `recalled` (snippet, and cited in Mallona et al.) |
| Siepe, Bartoš, Boulesteix, Heck, Peikert, Sarafoglou, Pawel 2026, "Why, when, and how to (or not to) preregister a simulation study", preprint | Argues the benefit of preregistration depends on the study's aim and its place in the research cycle | Same | `recalled` (snippet) |

### B4. Calibration of per-variant-pair phase confidence

| Work | What it does | Bearing | Status |
|---|---|---|---|
| Damaraju et al. 2026, PhaseQuality | Random forest classifiers per technology and reference, trained on a balanced subset of HG002 (all incorrect pairs and correct pairs at 2:1) with six features, among them the ratio of the maximum gap between variants to the longest aligned read in the phase set, the maximum gap, and difficult-region overlap. Output is three tiers (HighQ, MediumQ, LowQ), not probabilities. LowQ and MediumQ capture 82 to 99 per cent (ONT) and 91 to 99 per cent (HiFi) of errors; HighQ error rates are 0.005 to 0.35 per cent (ONT) and 0.05 to 0.46 per cent (HiFi); applied to HG005 for cross-sample generalisation | The closest competitor to the paper's calibration section. PhaseQuality ranks pairs into tiers from alignment features; it does not report calibration of a probability, and its features do not include the count of reads calling both variants. The paper should cite it and state the difference in these terms | `paper, searched` |
| LongPhase 2 manuscript repository (section H) | GNN that scores phased variants and unphases predicted misphasings | Unchanged from section H | `abstract` |
| pHapCompass (section H), arXiv 2512.04393v2 | Posterior over polyploid phasings; reports pairwise phasing uncertainty as Hamming distance between posterior samples | Uncertainty quantification, not a calibration measurement against truth in the parts read | as in section H |
| Kuleshov 2014, ProbHap | Posterior and transition scores per SNP (the web search surfaced the author's copy, which describes the transition score as the probability of correctly phasing a SNP relative to the previous one) | The older per-boundary probability the calibration section should compare against | citation verified in part A |

No 2023 to 2026 item found reports a reliability diagram or an expected
calibration error for per-pair phase probabilities from a read-based phaser.
This is a report of the searches listed above.

### B5. Pedigree or trio read-selection budgets for PedMEC

Nothing from 2023 to 2026 was found that studies how a coverage budget should
be split across the members of a pedigree for PedMEC. What was found:

- The WhatsHap v2.7 changelog entry above presents the heuristic phaser as a
  way to handle "larger pedigrees", and recommends the exact algorithm for
  pedigree phasing. This shows the maintainers treat pedigree size and
  coverage as the binding cost.
- The WhatsHap guide's pedigree section describes the `.ped` input and
  recombination costs and says nothing about how the cap is applied across
  samples. Issue 225 (section F of the matrix) is still the only source for
  the cap being shared across a pedigree.
- Garg et al. 2016 (Bioinformatics 32(12):i234, the PedMEC paper) was
  surfaced again. A search snippet, not the paper, says that 2x per member
  of a trio can match 15x on one individual; that is `recalled` until the
  abstract is read. The paper is outside the window and is presumably
  already the PedMEC citation.
- TinkerHap uses trio data for evaluation and can merge pedigree phasing, but
  it does not select reads.

A statement in the paper that pedigree budgets are an open direction is
consistent with what was found, and should be worded as future work rather
than as a gap in the literature.

## Items that change the prior-art matrix

These are recorded here for the owner; the matrix is not edited.

1. Damaraju et al. 2026 is a bioRxiv preprint. Section H's "now peer
   reviewed" is wrong, and the "journal name not retrieved" note is resolved:
   there is none.
2. The Q100 row can carry verified bibliographic fields (Hansen et al.,
   *Cell* 189(16):4857-4875.e31, 2026). Its content status stays `recalled`.
3. The Vikalo-group row is Si, Vikalo and Vishwanath, ITW 2014 (not ISIT),
   with a 2017 IEEE Trans. Inf. Theory version; status `abstract`.
4. Kamath, Şaşoğlu and Tse 2015 rises to `abstract`.
5. New rows: Ralphi (B1), Teo et al. 2024 and Mallona et al. 2025 (B3),
   Ben-Yashar 2023 (B2).
