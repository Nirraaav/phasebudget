# PhaseBudget

**Which reads to keep, and when to trust a phased genome.**

PhaseBudget is a Python library and a set of registered experiments on
long-read haplotype phasing. It asks two questions:

1. **Which reads should an exact phaser keep?** WhatsHap solves the phasing
   objective exactly, but only after cutting coverage to at most 15 reads per
   site with a greedy rule. PhaseBudget solves that cut exactly, as a min-cost
   flow, and measures what the greedy rule costs.
2. **When can a phased relation be trusted?** Phasers report a phase with no
   error rate. PhaseBudget marks each relation between neighbouring variants
   as trusted or not, with a threshold chosen so that, with 95 per cent
   confidence, at most a share $\beta$ of 100 kb stretches holds a wrong
   trusted relation, and it tests that promise on genomes it was not tuned on.

Course project for DH607, Koita Centre for Digital Health, IIT Bombay
(Autumn 2026), by Nirav Bhattad (23B3307).

| | |
|---|---|
| **Start here:** report for a general reader, with a primer on DNA and sequencing, figures and a glossary | [`docs/presentation/report.pdf`](docs/presentation/report.pdf) |
| Paper, the full technical account | [`docs/submission/paper.pdf`](docs/submission/paper.pdf) |
| Registrations | [`docs/design/`](docs/design/): each test's data, thresholds and decision rules, written before its test data were read; the commit times are in the private development repository (see [History](#history)) |
| Where each number comes from | [`docs/presentation/numbers.tsv`](docs/presentation/numbers.tsv): 74 of the report's main numbers, each with its run file and line; `python3 scripts/check_report_numbers.py` checks them |

## Results

All data are public. 48 people: the Genome in a Bottle references HG001,
HG002 and HG005, HG001's family of seven (Platinum Pedigree), and 39 Human
Pangenome (HPRC) genomes of African, admixed American and East Asian
ancestry, with PacBio HiFi and Oxford Nanopore reads at 10 to 30-fold.
**Registered** means the plan and decision rule were committed before the
test data were read (by the author's record; see [History](#history)); **development** means the data had been seen.

### Read selection

| Finding | Evidence |
|---|---|
| WhatsHap 2.8's shipped selector is not the procedure in its paper: it also lowers the scores of waiting reads | `select_whatshap`, a port of the shipped code, picks WhatsHap's reads on 1,500 of 1,500 inputs; the published rule on 178 |
| At the default cap of 15 the shipped selector reaches 0.933 (HG002 HiFi), 0.963 (HG002 ONT) and 0.934 (HG005 HiFi) of the exact optimum | min-cost flow, checked against subset enumeration |
| Inside `whatshap phase`, exact selection lowers ONT switch error by 42 per cent at cap 3 (5.204 to 3.029 per 10,000, $p = 8.4\times10^{-28}$) and changes nothing detectable at caps 5 and 15 | registered |
| Most of the larger gaps measured in the library were guesses at tied sites | post hoc: leaving ties unphased removed 77 to 85 per cent of the gap |

### Certified phase relations

| Finding | Evidence |
|---|---|
| One genome-wide threshold does not transfer to a new genome | registered: 18 of 18 cells violated on HG001 |
| One threshold per GIAB stratum held in the easy stratum in all 582 testable cells, across relatives, basecallers and 39 HPRC genomes of three ancestries, while trusting more than 99.5 per cent of easy relations | registered (four tests) |
| Certifying a new sequencing batch needs at least $\ln(1/\delta)/\beta$ labelled units of that batch, about 2,995 at $\beta = 10^{-3}$; one genome has about 2,900 segmental-duplication units | proved (Proposition `prop:floor` in the paper) |
| Trusting a relation only when a single-sample hifiasm assembly agrees, with six genomes pooled, held in all 62 testable cells on two new genomes and certified segmental duplications at $10^{-3}$ with 0.78 and 0.96 of their pairs trusted | registered (K1b) |
| Within-gene cis or trans calls at a per-call error of $10^{-2}$: 0.60 to 0.89 of pairs certified with assembly agreement; checking the assembly only at a call's two end variants certified 0.90 and 0.89 of pairs 10 to 100 kb apart on NA12878 | registered (K4b; L7 P4) |
| Pooled thresholds failed (163 of 360 cells) because the calibration truth labelled a relation only where an assembly and transmission agreed, which keeps easy relations | development (L1: 0 violating cells on key-labelled pairs, 120 when the same share is kept at random); replicated on 19 held-out HPRC genomes (registered L8: 0 of 2,052 on key links, 1,221 on all links) |
| Five predictions for NA12878's HiFi reads, written before her phase existed, all held | registered (L7) |
| Benchmarks scored inside the GIAB region miss most errors: 53 to 96 per cent of wrong relations lie outside it | registered on HG005 and NA12878; development on HG002 |
| Over 37 certificate tests, every violation came from test relations outside the population the calibration labels came from (or a new sequencing batch); every matched test had none | post hoc synthesis (L10), with two registered tests scored after it |

Tests that failed by their own rules are kept in the record: a majority-vote
model of where selection matters, error-optimal selection, a learned score
in segmental duplications, aligning to CHM13, siblings as labels, pooling
across sequencing batches, a pool-size planner, HiFi and ONT agreement, and
three-way agreement at $10^{-3}$. The paper's appendix "Other directions
tested" reports them.

## Install

```bash
git clone https://github.com/Nirraaav/phasebudget
cd phasebudget
python3 -m pip install -e ".[dev]"
make -C cpp            # optional: the native dynamic program, 30 to 100 times faster
python3 -m pytest
```

Python 3.10 or newer. The library core has no dependencies. Optional
extras: `ilp` (PuLP and numpy, for the integer-program solver; its tests
need a MILP solver such as CBC on the path), `whatshap` (WhatsHap 2.8, for
the live check of the selector port). Without `make -C cpp` everything runs
in pure Python and the native tests skip.

## Using the library

### Select reads exactly, then phase

```python
from phasebudget import CoverageDP, select_optimal, select_whatshap
from phasebudget.metrics import phasing_from_solution, summarize
from phasebudget.simulate import simulate_connected

data = simulate_connected(60, coverage=15.0, read_span=6, error_rate=0.05, seed=1)
data.matrix.max_coverage()             # too deep for the exact dynamic program

best = select_optimal(data.matrix, 5)  # min-cost flow: the heaviest selection with <= 5 reads per site
greedy = select_whatshap(data.matrix, 5)
greedy.objective / best.objective      # the shipped greedy's share of the optimum

subset = data.matrix.subset(best.selected)   # variant indices are kept
dp = CoverageDP(subset)
solution = dp.solve()                  # exact weighted minimum error correction
phasing = phasing_from_solution(subset, solution)
summarize(data.haplotype, phasing, data.positions)   # switch and flip errors, N50, ...
```

### Score each phased relation and choose a certified threshold

```python
from phasebudget.confidence import pair_support, score, learn_then_test_by_group

support = pair_support(data.matrix, phasing)     # reads agreeing and disagreeing with each relation
rows = []
for pair in support:
    called = phasing.alleles[pair.first] ^ phasing.alleles[pair.second]
    true = data.haplotype[pair.first] ^ data.haplotype[pair.second]
    unit = data.positions[pair.first] // 100_000
    rows.append(("easy", unit, score(pair), called != true))   # (group, unit, score, wrong)

learn_then_test_by_group(rows, beta=1e-3, delta=0.05)
# {'easy': None}: a toy instance has far fewer than the 2,995 units a promise at 1e-3 needs
```

On real data, `rows` come from calibration genomes with a truth set, one
group per GIAB stratum, and the returned thresholds are applied to new
genomes. `benchmarks/phase_confidence.py` scores every relation of any
phaser's VCF (only `GT` and `PS` are read) from a BAM:

```bash
PYTHONPATH=src:benchmarks python3 benchmarks/phase_confidence.py \
  --bam reads.bam --vcf phased.vcf.gz --reference GRCh38.fa --sample HG002 --out pairs.tsv
```

### What is in `src/phasebudget`

| Module | Contents |
|---|---|
| `fragments`, `_types` | the fragment matrix (reads as allele calls over global variant indices) and immutable result types |
| `exact` | exhaustive oracles for cost, constrained cost and margins; the reference every fast solver is checked against |
| `dp`, `_native`, `ilp`, `general` | the coverage-parameterised dynamic program (Python and C++), a 0-1 integer program, and the general bipartition model |
| `selection` | exact selection by min-cost flow, connectivity-aware selection, maximum-minimum coverage, and the published selectors (Fischer and Marschall as published, WhatsHap 2.8 as shipped) |
| `confidence` | Wilson scores, Learn then Test per group, the label floor, chain counts for within-gene calls |
| `adaptive`, `posterior`, `twolevel`, `hapchat`, `gapped`, `paired` | the adaptive reselection loop, sum-product margins, exact block orientation, HapCHAT-style read merging, the gapped coverage convention, paired selectors |
| `metrics`, `bootstrap`, `truth`, `simulate` | switch, flip and block metrics, the block bootstrap, truth handling, simulators |

`dp`, `ilp` and `exact` agree on the cost of every instance all three can
solve; `tests/` holds fixed examples and randomised differential tests
against the oracles. `CONTRIBUTING.md` states the correctness contract.

## Repository map

```text
src/phasebudget/     the library
tests/               fixed examples, differential tests against the oracles, pipeline checks
cpp/                 the native dynamic program
benchmarks/          one script per measurement; each prints its parameters and a table
benchmarks/results/  run outputs (runs/) behind every number in the paper
scripts/             data download with md5 checks, read subsampling, run drivers
docs/design/         problem statement, proofs, metric definitions, and every registration
docs/literature/     prior-art matrix with verification status, and the novelty gate
docs/presentation/   the report for a general reader (report.pdf) and numbers.tsv
docs/submission/     the paper (paper.tex, partA.tex, paper.pdf)
docs/followup/       the problem statement of the certificates, and the literature gate
PROPOSAL.md          the course proposal
```

## Reproducing the measurements

Data are not in the repository. `scripts/fetch_data.sh`,
`scripts/fetch_hg005.sh`, `scripts/fetch_platinum.sh` and
`scripts/fetch_pedigree_truth.sh` download the GIAB, Platinum Pedigree and
reference files and check them against the `scripts/*.md5` lists; HPRC and
Platinum Pedigree reads are streamed from their public buckets by
`benchmarks/genome_run.py` and `benchmarks/confidence_run.py`.
`scripts/reproduce_paper.sh` runs the read-selection measurements in order;
the certificate tests are run per registration, as each design note says.
`requirements-lock.txt` pins the environment the genome-wide runs used.
Each run file in `benchmarks/results/runs/` starts with a header describing
the run and its parameters; about 40 per cent also name the commit that
produced it. `docs/presentation/numbers.tsv` links the report's main numbers
to the run file and line they come from, and
`python3 scripts/check_report_numbers.py` checks that each one is still
there (70 of 74 point to a run file; the other 4 are sums or ratios of
numbers that do). The large runs took one to four days
each on a 96-core server; one hifiasm assembly took up to 67 hours.

## What is not claimed

- Min-cost flow for interval selection and Learn then Test are not new. The
  contributions are their application to phasing, the measurement of the
  deployed selector, the label floor, the account of why pooled calibration
  failed, and the registered tests.
- The guarantee assumes independent 100 kb units. In the hard strata wrong
  units cluster (181 of 576 tests); merging units to 200 and 300 kb changed
  no decision (development).
- Every lower-depth subsample in the read-selection runs used one random
  seed. Reruns with seeds 1 to 3 (`PHASEBUDGET_SEEDS`,
  `benchmarks/score_seeds.py`) are in progress.
- NA12877 and NA12878 are related to calibration genomes in some tests;
  HG005 is the unrelated test genome and has no indel truth. Truth sets have
  their own errors, and outside benchmark regions some counted errors may be
  the truth's.
- No claim depends on a negative claim about the literature
  (`docs/literature/02-novelty-gate.md`).

## History

This repository is a single-commit snapshot of the project at submission.
The development history is kept in a separate repository: the commits whose
times date each registration (the hashes quoted in `docs/design/`), the lab
notebook and planning memos that design notes refer to, earlier drafts, and
notes that served only the development machine. Copies of third-party papers
are not redistributed, which is why that repository stays private. Readers
therefore cannot check the registration times themselves: GitHub's push
times for it are saved in
`benchmarks/results/runs/2026-10-01_github_push_events.tsv`, but they come
from the author's account.

## Citation and licence

MIT licence. To cite, use `CITATION.cff`, or:

> Nirav Bhattad. *Read selection and certified phase relations in long-read
> haplotype phasing: exact optimization, error control within genomic
> strata, and why calibration labels must match the test population.* DH607
> course project, IIT Bombay, 2026. https://github.com/Nirraaav/phasebudget
