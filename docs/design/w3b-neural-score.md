# W3b: a neural per-pair score from read-level evidence

Written 15 October 2026, after W3 failed its gate on one genome and with
pooled calibration (development notebook, 15 October). Development work on data
already scored in earlier tests; nothing here is registered.

## What it optimizes

A score for each consecutive phased pair, `s = -p`, where `p` is a
network's estimate of the probability that the pair's relative phase is
wrong. The network is trained with class-weighted binary cross-entropy, so
it is fitted to rank wrong pairs below right ones. Learn then Test then
picks a threshold per group on `s` exactly as for the frozen Wilson score.
The guarantee needs only that the score is fixed before calibration, so the
network trains on chromosomes 1 to 6 and never sees the chromosomes that
calibrate or test it. Ties on `p` are broken by the Wilson bound, then the
number of reads, as W3's learned score does.

## Why it might help

W3's gradient-boosted score saw 22 per-pair summaries. In segmental
duplications a phase error often comes from reads of another copy aligned to
the wrong place. Such reads disagree with both haplotypes along several
variants, while a read with one sequencing error disagrees at one variant.
Summaries of the pair alone cannot tell these apart; the per-read pattern
across nearby variants can.

## Input

For a pair `(a, b)`, take the window of 8 phased variants around it in the
phasing's order: 3 before `a`, `a`, `b`, and 3 after `b`. Slots past the end
of the chromosome, or in another phase block, are empty. Take the reads that
span `a` or `b`; when there are more than 32, a sample of 32 with a seed
fixed by the pair's position. For each read and window slot:

- the call: +1 when the read's allele is the phaser's first-haplotype
  allele there, -1 when it is the other allele, 0 for no call;
- the call's weight from the fragment file, clipped to 255.

Per read, also the number of variants it calls in total. Per pair, W3's 22
features and the group one-hot, as in `benchmarks/w3_gate.py`.
`benchmarks/w3_reads.py` writes these arrays in the same pair order as
`benchmarks/w3_features.py`.

## Model and training

One model per depth and phaser. Each read goes through a two-layer
perceptron (width 64); the read vectors are pooled by mean and by maximum
over the reads present; the pooled vector and the pair features go through
a second two-layer perceptron to one logit. Adam, learning rate 1e-3, batch
4,096, 3 epochs, positive-class weight the ratio of right to wrong pairs
capped at 100, fixed seeds. Training data: chromosomes 1 to 6 of the seven
genomes of `benchmarks/w3_pooled.py`, pooled. `benchmarks/w3_neural.py`
trains on a GPU and writes `p` for every pair.

## Evaluation and gate

`benchmarks/w3_pooled.py` adds the neural score beside the frozen and the
gradient-boosted ones: each genome in turn is held out, thresholds are
calibrated on chromosomes 7 to 22 of the other six pooled, and scored on the
held-out genome against its truth. The gate is W3's with the neural score in
place of the learned one: at beta 1e-3, depths 30 and 20, both phasers, the
neural score trusts at least 1.5 times the frozen score's share of pairs in
`segdup` and in `indel` (and more than none), at most 0.005 less in `easy`,
and no cell violates at 0.05 over the cells scored.

## Limits

The seven genomes are HG002, HG001 and five of HG001's children, so the pool
is mostly one family, and HG002's reads come from a newer basecaller than
the Pedigree genomes'. W3 pooled already showed that pooling across those
basecallers breaks the bound when HG002 is held out; a better score cannot
fix that, so the gate is read with that failure in mind. The model trains
on the same genomes it is calibrated on, on other chromosomes.
