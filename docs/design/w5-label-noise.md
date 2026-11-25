# W5 part 3: calibrating on noisy transmission labels

Written 15 October 2026, after W5 part 2 (development notebook, 15 October). This is
development work on HG002 (`cal_v5`), whose phase is known from v5.0q. No
registered test reads it.

## The problem part 2 found

A child's heterozygous site gets a transmission label when a parent is
homozygous there. Part 2 took the parents' genotypes from ONT's Clair3 calls
and read a site with no record as `0/0`. In `easy`, at 30x whatshap, 2.1e-4
to 2.8e-4 of labelled pairs disagree with v5.0q, while whatshap's own error
on those pairs is 1.8e-5. Calibrating against the labels then counts label
errors as phasing errors, and abstains at beta 1e-3 and below.

## Why the error cannot be subtracted

Let a unit's label loss be 1 when it holds a trusted labelled pair whose
relation differs from the label, and its true loss 1 when it holds a
trusted pair that is wrong. Label loss exceeds true loss by the false
alarms (label wrong, phaser right) and falls short by the hidden errors
(label and phaser wrong together). Calibration on labels is valid when
hidden errors are rare, which they are when label and phaser errors are
independent; it is conservative by the false alarms. A unit in `easy` holds
about 40 labelled pairs, so a per-pair label error of 2e-4 puts a false
alarm in about 1 per cent of units. Certifying beta 1e-3 means telling a
unit risk of 1e-3 from 1e-2 plus a correction, and the correction would have
to be known to a few per cent of its value. A per-pair error rate estimated
from the data is not that precise. The labels have to be made cleaner
instead.

## What a label rests on

For a child's heterozygous site, with `f` and `m` the parents' calls:

- **double**: both parents homozygous for different alleles. The child's
  paternal allele follows from either parent alone (from the mother's call,
  the maternal allele is fixed and the child's other allele is paternal).
  The label is wrong only if both calls are wrong, or the child is not
  heterozygous, which is outside every guarantee here.
- **called**: one parent has a PASS `1/1` call and the other is
  heterozygous or uninformative. The label rests on one positive call.
- **absent**: one parent has no record, read as `0/0`, and the other is
  heterozygous or uninformative. A heterozygous site the caller missed in
  that parent makes the label wrong half the time.

Two flowcells of each parent were called separately (replicates A and B).
A label that both replicates give is wrong only when both callings err the
same way at that site.

## Filters

`benchmarks/w5_noise.py` reports, per group, the share and error of each
kind of site label, and then calibrates on the labelled odd-chromosome
pairs kept by each filter, scoring all even-chromosome pairs against v5.0q:

| filter | keeps a site label when |
|---|---|
| `all` | replicate A labels it (part 2) |
| `agree` | both replicates label it, with the same allele |
| `positive` | replicate A labels it and the kind is not `absent` |
| `agree_positive` | both replicates label it, the same, and neither rests on an absent record |
| `double` | replicate A labels it as `double` |

A filter that keeps few pairs has a second problem. The loss counts a unit
as wrong when any trusted pair in it is wrong, and calibration sees only the
labelled pairs. With most pairs labelled (part 1: 0.78) the two agree; with
a few per unit, calibration sees a lower unit risk than the test does, and
the bound can fail. The filters' labelled shares are reported for that
reason.

## Gate (development)

A filter passes when, pooled over depths 30, 20 and 10 and both phasers, no
cell violates at 0.05 over the number of filter cells, and at beta 1e-3 the
`easy` group certifies in all 6 cells with a trusted share within 0.01 of the
truth-calibrated oracle. A passing filter is a candidate for a registered
test on a family not seen here. Part 1's five CEPH children, whose parents
have ONT reads in the Platinum Pedigree bucket, are that family.
