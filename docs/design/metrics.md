# Metric definitions

Switch error, flip error, and the N50 family are all reported under more than
one convention in the phasing literature. This is the one used here.

## Comparability

Phase is only comparable inside a block. `phasebudget.metrics` compares
adjacent phased variants within a block and nothing else. A block of one
variant contributes nothing, and neither does a variant no selected read
covers.

## Switch error

For consecutive phased variants `p` and `q` in a block, the predicted relation
is `alleles[p] XOR alleles[q]` and the true relation is `truth[p] XOR truth[q]`.
A switch is a disagreement, and it is recorded at `q`. The switch error rate is
the number of switches divided by the number of comparable adjacent pairs.

A block phased entirely opposite to the truth has no switch errors. That is the
point of measuring relations rather than alleles: which haplotype is called
first is not part of the answer.

A real truth set does not phase every variant; the transmission truth of
`docs/design/truth-sets.md` phases 79.5 per cent of HG002's heterozygous SNVs
on chromosome 20. A variant the truth leaves as `None` is skipped: the pairs
compared are consecutive variants of a predicted block that the truth also
phases, so a skipped variant can neither cause nor hide a switch between its
neighbours. `comparable_pairs(phasing, truth)` is the denominator under that
rule. This is how `whatshap compare` counts, and with a complete truth it
changes nothing.

## Flip error

Two switches at consecutive comparable pairs put one variant on the wrong
haplotype and put the phase back afterwards. Reporting that as two switch
errors overstates the damage, so it is reported once, as a flip. Pairing is
greedy from the start of each block and non-overlapping, so a run of three
consecutive switches yields one flip and one leftover switch rather than two
overlapping flips. `switch_errors_excluding_flips` reports the switches that no
flip consumed.

## Phased fraction

The fraction of variants that lie in a block of at least two phased variants. A
variant alone in its block carries an allele but no phase information.
Counting it would let a method raise its phased fraction by emitting singleton
blocks.

## N50, AN50, QAN50

`block_n50` is the N50 of block spans: the largest `L` such that blocks of
length at least `L` cover half the total phased length. Lengths are variant
counts unless variant positions in base pairs are supplied.

`an50` scales each block's span by the fraction of the variants inside that
span which the block actually phases, so a block that skips half the variants
between its endpoints is credited with half its span.

`qan50` is `an50` computed after every block has been cut at each of its switch
errors. It is the number that stops a method from buying N50 by extending
blocks through regions it cannot phase: a long block with a switch in the
middle scores as two short ones.

## The four levels, and why they stay apart

A read-selection method is measured at four levels, and a method that wins one
can lose another. They are reported in separate columns and never combined.

1. **Selection objective.** What the selector scored, under the base weights.
   Reported as a fraction of the certified optimum from `select_optimal`, so
   the quotient is an optimality gap. A selector that maximizes something else
   also gets a column for its own objective, `min_coverage` for the
   maximum-minimum-coverage selector and `distinct_variants_covered` for the
   Fischer and Marschall slice, so it is not scored only against a target it
   was not aiming at.
2. **Runtime.** Selection time and solve time, separately. Exact selection is
   worth nothing if it costs more than the solve it enables, and the two are
   different orders of magnitude at different caps.
3. **wMEC cost**, scored against the full read set.
4. **Switch error**, with the phased fraction beside it.

### wMEC cost is comparable only on a fixed read set

The cost of a haplotype on the reads a selector retained is not a quality
measure and must never be reported as one. Discarding a read removes its
corrections from the sum, so keeping fewer reads lowers the cost by itself: in
`benchmarks/benchmark_heuristic_gap.py` random selection keeps 26 reads and
scores 19.27 on them, against 21.94 for exact selection on 30 reads, while
phasing three times worse. Every cross-selector wMEC comparison here scores the
reported haplotype against every read, by `FragmentMatrix.cost` on the full
matrix. The retained-set number is printed next to it so the difference is
visible rather than argued.

The same rule applies to the selection objective under reweighting. The
adaptive loop changes the weights, so its selections are scored under the base
weights whenever it is compared with anything else.

## What is deliberately absent

There is no composite accuracy-completeness score. A method that phases half
the variants perfectly and a method that phases all of them with a few switches
are different answers to different questions, and averaging them into one
number decides in advance which one the reader wanted.
