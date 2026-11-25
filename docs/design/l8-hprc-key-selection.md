# L8: L1's key-selection prediction on 19 held-out HPRC genomes

Rules written 24 November 2026, after the development run on the 20 W1b
pool genomes (`2026-11-22_l8_dev.txt`) and before any held-out genome's
phase was compared with its dipcall truth alone. The held-out genomes were
scored in W1 and W1b against the consensus truth (dipcall equal to 1000
Genomes transmission) only; the quantities below, which need the dipcall
truth on links the consensus does not label, have not been computed for
them.

## Claim

L1 found that the consensus key labels easy links, so thresholds calibrated
on it hold on key-labelled links and fail on all links. If that is the
mechanism, it holds on new genomes of other ancestries.

## Data and computation

`benchmarks/l8_key_selection.py`, `MODE=test`: the 19 held-out HPRC genomes
of W1b (`w1b_score.HELD_OUT`), HiFi, depths 30, 20, 10, whatshap and
LongPhase, 21 chromosomes. Truth: each genome's HPRC release 2 dipcall alone
(`PHASEBUDGET_TRUTH=dipcall`). Links where dipcall and transmission both
phase but disagree are left out and counted. Classes: `all` (dipcall-labelled
links), `key` (links the consensus labels), `random` (links kept at random
with the key's share per cell, seed 0). Thresholds: W1b's pooled thresholds.
Cells: (genome, depth, phaser, beta, group), 19 x 108 = 2,052, level
0.05/2,052 = 2.44e-5. Ratio: wrong rate on all links over wrong rate on key
links, over `segdup`, `lowmap` and `tandem` together.

## Development numbers (in sample: W1b's thresholds were fitted on these genomes)

Ratio at 30x above 3 in 20 of 20 genomes for both phasers (smallest 15.1
whatshap, 30.5 LongPhase). Violating cells at 0.05/2,160: all 1,247, key 0,
random 1,195; every genome had a violating hard or indel cell on all links.
Links where dipcall and transmission disagree are 0.03 to 0.10 per cent of
links over the 120 (genome, depth, phaser) rows, and the phaser is wrong
against dipcall on 8.2 per cent of them, so transmission, not dipcall,
carries most of that disagreement.

## Decision rules (held-out genomes)

- **L8-1.** At 30x the ratio is above 3 in at least 17 of 19 genomes, for
  each phaser.
- **L8-2.** On all links, at least 15 of 19 genomes have a violating cell in
  `segdup`, `lowmap`, `tandem` or `indel` at 2.44e-5.
- **L8-3.** On key links, no cell violates at 2.44e-5.

L8 is supported when all three hold. The random control and the
disagreement counts are reported and decide nothing.

## Known limit

On links transmission does not label, the dipcall truth is the only truth,
and its own errors would count as phaser errors there. The development
disagreement counts bound how often dipcall and transmission differ where
both phase; they do not bound dipcall's error where transmission is
missing. The random control separates selection from lost power, not from
truth error.

## Result (24 November 2026, `2026-11-24_l8_test.txt`)

All three rules hold, so L8 is supported.

- **L8-1 holds.** At 30x the ratio is above 3 in 19 of 19 genomes for both
  phasers; smallest 14.28 (whatshap, HG02257) and 31.92 (LongPhase, NA19240).
- **L8-2 holds.** On all links 1,221 cells violate at 2.44e-5, and every one
  of the 19 genomes has a violating `segdup`, `lowmap`, `tandem` or `indel`
  cell.
- **L8-3 holds.** On key links no cell violates.

The random control, which keeps the key's share of links per cell at random,
has 1,159 violating cells and 19 of 19 genomes with a violating hard or indel
cell. So the key's 0 comes from which links it keeps, not from how many.
Links where dipcall and transmission both phase but disagree are 0.0586 per
cent of dipcall-labelled links over the 114 (genome, depth, phaser) rows
(range 0.038 to 0.084 per cent), and the phaser is wrong against dipcall on
7.5 per cent of them, close to the 8.2 per cent in development. The known
limit above still applies.
