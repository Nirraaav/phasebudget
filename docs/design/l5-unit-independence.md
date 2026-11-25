# L5: are 100 kb units independent? (round 8, 22 Nov 2026)

## Claim

Every Learn then Test certificate in this repository treats the 100 kb units
of one group as independent draws: the binomial bound counts wrong units
among `n` units as if each unit were wrong independently. Neighbouring units
share reads, repeats and reference problems, so wrong units may sit next to
each other. The only earlier check (`benchmarks/gate_bin_dependence.py`) used
10 kb bins of one HG002 HiFi phasing, pooled over groups. This test asks the
question per group on every genome whose truth labels the whole genome.

The data are already seen. Every row cache below was read by an earlier
registered or development scoring, so this is a development check with its
rules written before it runs, not a blind test.

## Genomes, platforms and truths

The cached rows of `benchmarks/confidence_score.load_rows` (consecutive
pairs, both variants labelled by the platform's primary truth), 21
confirmatory chromosomes, depths 30, 20 and 10 where cached:

| platform | genome, reads | truth | depths |
|---|---|---|---|
| `hprc_HG002` | HG002 HiFi | GIAB v5.0q | 30, 20, 10 |
| `w4_grch38` | HG002 HiFi, minimap2 | GIAB v5.0q | 30 |
| `cal_v5` | HG002 ONT sup v5 | GIAB v5.0q | 30, 20, 10 |
| `cal_r10` | HG002 ONT R10 | GIAB v5.0q | 30, 20, 10 |
| `hg001_v5` | NA12878 ONT | Platinum Pedigree | 30, 20, 10 |
| `ped_<s>` | NA12879, NA12881, NA12882, NA12885, NA12886 ONT | Platinum Pedigree | 30, 20, 10 |
| `phifi_<s>` | the same five, HiFi | Platinum Pedigree | 30, 20, 10 |
| `phifi_NA12877` | NA12877 HiFi | Platinum Pedigree | 30 |
| `k1_hg005` | HG005 HiFi | GIAB trio transmission | 30 |
| `hg005_ont` | HG005 ONT R10 | GIAB trio transmission | 30, 20, 10 |

HPRC genomes scored against the consensus truth are left out: that truth
labels a selected subset of pairs, which is idea 1's subject. Phasers:
whatshap and LongPhase. Groups: `indel`, `segdup`, `lowmap`, `tandem`,
`difficult`, `easy`, assigned as in `stratified_score.grouped`. A test is one
(platform, depth, phaser, group): 96 (platform, depth, phaser) combinations
times 6 groups gives 576 tests.

## Statistic

Within one test, a unit is `(chromosome, position // 100 kb)` holding at
least one pair of the group, and it is wrong when any of its pairs is wrong
(all pairs, no threshold). Two units are adjacent when they are on the same
chromosome and their indices differ by 1. The statistic is the number of
adjacent pairs with both units wrong.

The null keeps each chromosome's units and its number of wrong units, and
permutes the wrong labels among that chromosome's units. The p-value is
`(1 + #{null >= observed}) / (1 + B)`, with B = 2,000 permutations, raised to
B = 40,000 when the first p is below 0.01. Seed 0.

Reported with each test: units, wrong units, adjacent pairs, observed and
null mean joins, their ratio, the lag-1 correlation `rho` of the wrong
indicator over adjacent pairs, and the effective sample size
`n (1 - rho) / (1 + rho)` (`n` when `rho <= 0`).

## Decision rules

The level is 0.05 / 576 = 8.68e-5 (Bonferroni over every test, including
tests with too few wrong units to reach it). A test is clustered when its p
is at most that level.

**Merge check.** For every group, the unit is widened to 200 kb and 300 kb
(index `// 2` and `// 3`) and the thresholds of `learn_then_test_by_group`
(betas 3e-4, 1e-3, 3e-3, delta 0.05) are recomputed on two pools:

- the K1b pool (`benchmarks/k1b_score.py`): HG002 and the five Pedigree
  members, HiFi 30x, agreeing pairs and all pairs;
- the W1b pool (`benchmarks/w1b_score.py`): HG002 v5.0q and 20 HPRC genomes
  on their consensus truth, depths 30, 20, 10.

A merged cell is unchanged or more conservative when its threshold is None or
at least the 100 kb threshold (trusts no more pairs).

**L5 is supported** when either

- (a) no test is clustered, or
- (b) every clustered test lies in a set of named groups, and in each of those
  groups every merged cell of both pools, both phasers and every beta is
  unchanged or more conservative.

Otherwise L5 is not supported, and the certificates of the clustered groups
rest on an assumption the data contradict. The worst effective sample size
ratio `n_eff / n` over all tests is reported whatever the verdict.

## Limits

A unit is wrong on all pairs here, not on trusted pairs at a certified
threshold; trusted pairs have fewer wrong units, so clustering among them is
harder to detect and is not tested. Adjacency at lag 1 misses dependence
that skips a unit. The permutation keeps each chromosome's wrong count, so
dependence between chromosomes is not tested.
