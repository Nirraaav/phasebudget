# K1a and K4: certificates that require an assembly to agree

Written 19 October 2026, before the scorer ran. Development gates, not
registered. Directions K1 and K4 of the October directions memo (development repository),
round 7.

## Question

Q6 found that a HiFi-only hifiasm assembly of HG002 is right on most of the
hard-group relations whatshap gets wrong, and that its own errors are
switches in `easy`. K1a asks whether a certificate that trusts a relation only
when the read score passes and the assembly phases the pair the same way
certifies the hard groups at `beta = 1e-3`. K4 asks whether such links make
certified cis or trans calls in clinical genes, which Q4 could not because
nearly every chain holds an indel link.

## Material

Everything Q6 used, nothing new: HG002 HPRC HiFi realigned with minimap2
(`w4_grch38`), 30x, whatshap and LongPhase; the assembly tables of
`benchmarks/q6_assembly_phase.py` (`data/r6_q6/tables/`); the cached W4 rows
(frozen score, `wrong` against v5.0q); v5.0q and GIAB trio transmission. The
data were seen in Q6, so this is development.

## The certificate

A pair of consecutive phased variants **agrees** when the assembly relates the
two variants (Q6's contig rule) and gives the same relation as the phaser. The
agreement certificate runs `learn_then_test_by_group` on the agreeing pairs
only; a pair that does not agree is never trusted. The filter uses no truth,
so the Learn then Test guarantee holds for the trusted pairs as before: a unit
is wrong when it holds a wrong trusted pair, and pairs that are never trusted
cannot make it wrong. The comparison is the frozen-score certificate
calibrated on all pairs, as in Q6's hybrid table.

## K1a gate

Calibration on the 10 odd confirmatory chromosomes, test on the 11 even ones,
`beta` in (3e-4, 1e-3, 3e-3), `delta = 0.05`, 100 kb units, groups as in Q6.

- **(a) validity:** among all test cells (phaser, beta, group) with a
  threshold and at least one trusted unit, none has `P(X >= k | n, beta)`
  below 0.05 divided by the number of such cells.
- **(b) usefulness:** in each of `segdup`, `tandem` and `indel`, at least one
  phaser certifies at `1e-3` and trusts at least half of that group's even
  pairs (all pairs the phaser relates and v5.0q phases, agreeing or not).

K1a passes when (a) and (b) hold. Reported, not gated: on all 21 chromosomes,
per phaser and group, the pairs wrong in both the phaser and the assembly
(for a two-way relation, both wrong means they agree), and the share of the
group's units that hold one; above `1e-3` no threshold can certify that group
at `1e-3` for that phaser. Also the error of agreeing pairs against the trio
truth, since v5.0q derives from an assembly.

## K4 gate

Q4's call and certificate (`docs/design/q4-clinical-cis-trans.md`) on
`w4_grch38` 30x, with one change: a link the assembly does not agree with gets
a score below every threshold, and thresholds come from the agreeing pairs
only. Two folds by chromosome parity, targets `r` in (1e-3, 1e-2). The same
run without the change is the baseline on the same alignments.

- **(a) validity:** in each (phaser, target) cell, one certified pair per gene
  chosen by CRC as in Q4; no cell has `P(X >= k | n, r)` below 0.05/4.
- **(b) usefulness:** at `r = 1e-2`, at least one phaser certifies at least
  half of the within-gene pairs.

K4 passes when (a) and (b) hold.

## Code

`benchmarks/k1_assembly_agreement.py`; `tests/test_k1_assembly_agreement.py`
checks the agreement filter and the link rewrite.
