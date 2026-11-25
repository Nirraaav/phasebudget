# Q2: flags for hard-region errors from the sample's own reads

Written 17 October 2026, before any Q2 computation ran. Development gate,
not registered. Direction Q2 of the October directions memo (development repository)
(round 6). It follows W4 (`docs/design/w4-reference-gate.md`), whose
confidently wrong relations survived a change of reference, and Z7, whose
flip flag found 6 to 16 per cent of genotype errors.

## Question

In the hard groups (`segdup`, `lowmap`, `tandem`, `difficult`, `indel`) the
per-group certificates abstain or violate. A flag is a rule computed from
the sample's own reads with no truth. Do flags find most of the wrong
relations in a hard group at a small cost in pairs, so that a certificate
calibrated and scored on unflagged pairs certifies the group without
violation?

## Data

A relation is a consecutive phased pair (the rows of
`benchmarks/confidence_score.py`): two heterozygous variants the phaser puts
in one block, both phased by the truth. Groups are those of
`benchmarks/stratified_score.py`: `indel` when either variant is an indel,
else the stratum of `stratified_gate.strata`.

| Platform | Reads | Depths | Truth | Role |
|---|---|---|---|---|
| `hprc_HG002` | HG002 HPRC HiFi, winnowmap | 30, 20, 10 | v5.0q | gate |
| `w4_grch38` | the same reads, minimap2 | 30 | v5.0q | descriptive; W4's confidently wrong relations |
| `cal_v5` | HG002 ONT sup v5.0.0 | 30, 20, 10 | v5.0q | descriptive; tuning source for the transfer |
| `hg001_v5` | HG001 ONT sup v5.0.0 | 30, 20, 10 | Platinum Pedigree | transfer |

Phasers: `whatshap phase` and LongPhase, as stored. Odd confirmatory
chromosomes tune and calibrate; even ones (chr20 excluded) are scored.
No W1 or W1b test genome is read.

## Features

All are computed from the stored fragments of the depth being phased
(`fragments_<depth>_realign.json.gz`: reads with mapping quality at least 20,
no supplementary alignments, alleles as `whatshap phase` reads them), except
the k-mer feature, which reads the whole HG002 HiFi read set.

1. **Depth.** At a variant, the reads whose span covers it (gaps included),
   divided by the median of that count over the chromosome's variants at the
   same depth. Pair value: the larger `|log2 ratio|` of its two variants. A
   collapsed duplication raises it; reads lost below the mapping-quality
   filter to an unlisted paralog lower it.
2. **Allele balance.** At a variant, reads calling each allele; the
   two-sided binomial p-value of the smaller count under a fair coin. Pair
   value: the larger `-log10 p` of its two variants. A paralogous variant
   seen through a collapsed copy gives about one quarter or three quarters.
3. **Window MEC.** The optimum wMEC over a window of columns: the pair's two
   variants and the two variant columns on each side outside the pair
   (`W = 2`, at most 6 columns, fixed now), every read with at least two
   called alleles there, weights as stored, complementary model. Divided by
   the total weight of those reads' calls in the window, so it is the
   smallest weighted share of calls that no two complementary haplotypes
   explain. Pair value: that share. The definition of the optimum is
   `exact.solve_bruteforce` on the window's matrix; the vectorised
   implementation answers to it (see Tests).
4. **K-mer copy number** (HG002 HiFi platforms only). For each variant, the
   reference 31-mer ending just before it and the one starting just after
   its reference allele. Each is counted, canonically, in all reads of
   `data/w4/reads.fq.gz` and in `GRCh38_no_alt_analysis_set.fna` with
   jellyfish 2.3.0 (static binary, `--if` restricted to these k-mers). Copy
   excess of a k-mer: read count over (the median read count of all such
   k-mers times its reference count). Pair value: the largest `log2` excess
   of its four k-mers. A k-mer containing `N` is missing. Only gains are
   used: losses are dominated by nearby variants inside the k-mer. This
   feature cannot transfer to HG001 here, because HG001's ONT reads are not
   held locally; if jellyfish fails, the gate runs with features 1 to 3.

A missing value is treated as not extreme.

## Flag rule and tuning

Per platform, depth, phaser and hard group, on odd chromosomes only:

- a pair's **tail rank** under feature `f` is the share of the group's odd
  pairs whose value is at least as large;
- for a feature subset `S`, a pair is **flagged** when its smallest tail rank
  over `S` is at most `tau`; `tau` is the largest value at which at most
  `A = 0.09` of the group's odd pairs are flagged (0.09 leaves room under the
  0.10 limit);
- the subset is the one with the highest odd-chromosome **recall of wrong
  relations** among the non-empty subsets of the available features, ties
  to the smaller subset, then to the earlier feature order above.

On even chromosomes the same subset and `tau` apply, with tail ranks taken
against the odd pairs' values (a fixed cutoff per feature). Tuning uses the
truth of odd chromosomes; the flag is then a fixed function of the reads.
The `easy` group is never flagged.

**Transfer.** `hg001_v5` uses the subset and `tau` that `cal_v5` chose for
the same depth, phaser and group; tail ranks are against HG001's own odd
pairs, which needs no truth.

## Certificates

Per cell (platform, depth, phaser): `learn_then_test_by_group` at
`delta = 0.05` and `beta` in {3e-4, 1e-3, 3e-3}, 100 kb units, grid from 0.5
in steps of 0.0025, as in every earlier gate.

- **Baseline:** calibrated on all odd pairs, scored on all even pairs.
- **Flagged:** calibrated on unflagged odd pairs, scored on unflagged even
  pairs. Flagged pairs are abstentions: the trusted share is reported over
  all of the group's even pairs.

The statement is R3's: calibration and test see the same population (the
unflagged pairs), so Learn then Test's guarantee covers unflagged pairs of
exchangeable units and nothing about flagged ones.

**Unit floor.** With no wrong unit, certifying at `beta = 1e-3` and
`delta = 0.05` needs at least 2,995 calibration units (`(1 - beta)^n <=
delta`). A group with fewer units holding its pairs on the odd chromosomes
cannot certify at 1e-3 whatever the flags do. The output prints each group's
odd units; I expect `segdup` to fall below the floor (W3 pooled had about
17,000 `segdup` units over six ONT genomes).

## Gate

Fixed before looking. On `hprc_HG002`, even chromosomes, the six depth and
phaser cells:

- **(a) the flag finds most errors cheaply**, per hard group: pooled over the
  six cells, recall of wrong relations at least 0.5; recall of W4-style
  confidently wrong relations (at least 10 reads, none disagreeing) at
  least 0.5 where the group has at least 10 of them pooled (otherwise
  reported, not judged); abstention at most 0.10 in every cell;
- **(b) valid:** no violating flagged cell over the hard groups, three betas
  and six cells (90 cells), at level 0.05/90;
- **(c) useful:** per hard group, among the depth and phaser cells where the
  baseline certified nothing at `beta = 1e-3`, the flagged procedure
  certifies in at least half.

Q2 passes when (b) holds and at least one hard group meets both (a) and
(c). The output reports, per group and cell, what the tuning optimised
(odd recall at the chosen `tau`, and the chosen subset) separately from what
was achieved on even chromosomes (recall, abstention, violations, certified
share), with the worst cell beside the pooled value. It also reports each
single feature's even recall at its own odd-tuned cutoff, to say which flag
carries the signal.

**Transfer, reported, not gated:** HG001 recall and abstention with
`cal_v5`'s settings, and HG001's baseline and flagged certificates calibrated
on its own odd chromosomes.

## Tests and code

`tests/test_q2_hard_flags.py`: window MEC against `exact.solve_bruteforce`
on fixed matrices and on random small windows; the window's columns at
chromosome ends and for non-adjacent pairs; tail ranks and the choice of
`tau` under ties; the allele-balance p-value on fixed counts; the k-mer
positions for an SNV and an indel. Implementation in
`benchmarks/q2_hard_flags.py`.

## Limits

One genome tunes and is scored (HG002); the transfer genome is on ONT and
lacks the k-mer feature. Tuning and calibration share the odd chromosomes,
so the odd-chromosome numbers are optimistic; only even-chromosome numbers
are evidence. The flags are judged against v5.0q and the Pedigree truth;
where a truth is wrong (R1 could not reach most hard-group relations), a
"wrong" relation the flags miss may be right.
