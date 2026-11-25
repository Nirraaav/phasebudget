# Accepting the adaptive loop's rounds window by window

## The problem

`adaptive_selection` keeps the round whose haplotype has the lowest wMEC
against the full matrix. At 60 variants that is the right rule. At chromosome
scale a reweighting round repairs some regions and disturbs others, and a
round is kept only if the repairs outweigh the disturbances over the whole
chromosome. The development notebook records the consequence: with the scale-free fragile
rule the loop improved on static selection in 6 of 9 cells, by small amounts,
and a traced instance showed round 1 lowering switch error while later rounds
raised the global cost.

## Stitched acceptance

Let `H_0, ..., H_R` be the haplotypes of the rounds, each the optimum of its
own selected subset. Cut the variants into consecutive windows of `W`
variants. For each window `w` and round `r`, the local cost `L_w(r)` is the
full-matrix wMEC, minimized over read assignment, of the reads whose span lies
inside `w`, scored against `H_r` restricted to `w`. The stitched haplotype
takes, in each window, the round with the smallest local cost, ties to the
lower round index.

Pieces from different rounds meet at window junctions, and each piece is
defined only up to complementation, so a junction needs a relation. Windows
are joined left to right: the piece for window `w+1` is taken as is or
complemented, whichever gives the reads that span the junction the lower
full-matrix cost over the variants already built plus window `w+1`. Reads
that span the junction carry the only evidence about it, so this is the
evidence-based choice.

Last, the stitched haplotype's full-matrix wMEC is compared with the best
round's, and the lower is returned. Stitching therefore never reports a
haplotype worse under the full matrix than the global rule would, which
`tests/test_adaptive.py` checks; whether it phases better is a measurement.

## What is reported

The result's `solution` holds the returned haplotype, its full-matrix cost,
and the full-matrix best assignment. `support` is the union of the reads
selected by the rounds that supplied a window; phase blocks for evaluation
come from `support`, since those are the reads whose evidence the haplotype
rests on. With the global rule, `support` is the kept round's selection, so
the two rules are evaluated the same way.

## Limits

The windows are fixed, not chosen by the data, and `W` is a parameter. It
was 50 until 26 September. The milestone-4 sweep found an interior optimum
near 25 on the tuning seeds, and a held-out run on fresh instances (seed base
20000) confirmed it in three of four cells and tied in the fourth, so the
default is now 25. Both runs are in the development notebook. The
stitched haplotype is not the wMEC optimum of any read set; it is a
composition of optima scored against all the reads, which is the same
yardstick the global rule uses.

## A failure found on HG002

Stitching orients pieces only at window junctions. The haplotype of each
round is the optimum of its own selected reads, and wherever those reads
leave a boundary unspanned, the round's relation across it is arbitrary.
When such a break lies inside a window, the junction test cannot flip there;
the reads that cross the next junction also cross the break, and the cheaper
orientation for them can be to flip at the junction. The full-matrix cost
falls, the relation at the break stays wrong, and a new wrong relation
appears at the junction, inside a block of the reads the result reports.

On HG002 chromosome 20 at cap 3 (development notebook, M5 weeks 19 to 21) every one
of the loop's 15 extra switch errors sat on a junction, and the break was
found within the two windows before each of the four examined. Every window
had taken the same round, so the stitched haplotype differed from that round
only by junction flips. Simulation did not show this because the loop's
gains there were larger than the damage, and at caps of 5 and above the
selected reads rarely leave breaks.

The repair, orienting at every unspanned boundary of the chosen rounds as
well as at junctions, is a new rule. It has to be stated, checked against
an oracle and measured on simulation before it is used on real data, and it
does not change the registered HG002 result.
