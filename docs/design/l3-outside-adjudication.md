# L3: are errors outside the benchmark the phaser's or the answer key's?

Written on 22 November 2026, before any pair was compared with a parental
label. Development analysis with rules stated in advance, not a registered
test: Q3 (HG002, 17 October) and K3b (HG005, 27 October) are already scored,
so their error counts are known. Round 8, idea 3 of the five-idea plan.

## Claim

Q3 and K3b counted a pair as wrong when the phaser's relation differs from an
assembly-based key: Q100 for HG002, the HPRC release 2 dipcall for HG005.
Outside GIAB's benchmark BED most wrong pairs lie (K3b: 24,818 of 25,962
WhatsHap and 6,346 of 6,747 LongPhase wrong pairs on HG005), but outside the
BED the key is less certain, so some of these may be key errors. Phase labels
from the parents' reads are independent of the child's reads and of the
child's assembly. Where they label a disputed pair, they decide which side is
right. L3 claims that most outside disputes are phaser errors.

## Data

- **HG002.** Platform `q3_q100`: HPRC HiFi at 30x, W4's minimap2 GRCh38
  alignments, phased on Q100's heterozygous variants, rows cached by Q3 in
  `data/r6_q3/rows/`. Key: Q100 phase. Inside: both variants in
  `HG002_GRCh38_v5.0q_smvar.benchmark.bed`.
- **HG005.** Platform `k3_hg005`: HPRC HiFi at 30x phased on the release 2
  dipcall variants, rows cached by K3b in `data/k2/rows/`. Key: the same
  dipcall phase. Inside: both variants in
  `HG005_GRCh38_1_22_v4.2.1_benchmark.bed`.
- **Phasers.** WhatsHap 2.8 and LongPhase 2.0.2, as stored.
- **Pairs.** The cached rows: consecutive phased pairs whose two variants the
  key phases, 21 confirmatory autosomes. Group: `indel` when either variant
  is an indel, else the stratum of `stratified_gate.strata` (`segdup`,
  `lowmap`, `tandem`, `difficult`, `easy`), as in Q3 and K3b.
- **Labels.** ONT `wf-human-variation` sup calls (Clair3) of both parents,
  two flowcell replicates each (`data/w5/`, `w5_called.REPLICATES`): HG003
  and HG004 for HG002, HG006 and HG007 for HG005. A child site gets the
  paternal allele by `w5_called.paternal`'s rule; a site is labelled under
  R3's `agree` filter (both replicates give the same paternal allele). A
  pair is labelled when both sites are; its label relation is whether the
  two paternal alleles are equal.

Why these labels and not GIAB's transmission truth for HG002: the GIAB
transmission VCFs are restricted to the v4.2.1 benchmark, which lies inside
v5.0q's BED, so they label almost no outside pair. The parents' ONT calls are
genome-wide.

## Computation

Per genome, phaser, side (inside, outside) and group:

- `pairs`, `disputed` (phaser differs from the key), labelled share of each.
- **Label error where both sources agree, `e`.** Among labelled undisputed
  pairs, the share whose label differs from the phaser and the key. Two
  independent sources agree there, so `e` estimates the label's own error.
- **Verdicts on disputes.** Among labelled disputed pairs, `k_key` (label
  sides with the key, a phaser error) and `k_ph` (label sides with the
  phaser, a key error); `f = k_key / (k_key + k_ph)`.
- **Corrected phaser share.** If the label is wrong at rate `e`
  independently of which side is right, `f = p(1 - e) + (1 - p)e`, so
  `p = (f - e) / (1 - 2e)`, clipped to [0, 1]. The 95% interval combines
  Wilson intervals: lower `(f_lo - e_lo) / (1 - 2 e_lo)`, upper
  `(f_hi - e_hi) / (1 - 2 e_hi)`, each clipped. For `f > 1/2`, a smaller `e`
  gives a smaller `p`, so `e_lo` is the conservative end.
- **Check of the verdict where the key is trusted.** HG002 inside v5.0q's
  BED, Q100 is backed by GIAB's curation. There the share of labelled
  disputed pairs whose label sides with Q100 is the label's accuracy on
  disputed pairs, `a_in`. It tests the independence assumption behind `p`:
  if labels fail more often at disputed sites, `a_in` falls well below
  `1 - e`.
- **Coverage.** Labelled share of disputed outside pairs per group, and a
  group-reweighted phaser share: per-group `p` weighted by the group's
  disputed outside pairs, with groups of fewer than 20 labelled disputed
  pairs counted at 0 (the worst case for the claim).

## Decision rules

- **L3-0, the verdict is usable.** On HG002 inside v5.0q's BED, pooled over
  groups, `a_in >= 0.9` for both phasers. If L3-0 fails, L3-1 is reported as
  not decidable.
- **L3-1, primary.** On HG005 outside the GIAB BED, pooled over groups, the
  lower 95% bound of the corrected phaser share `p` is at least 0.5 for both
  phasers.
- **L3-2.** The same on HG002 outside v5.0q's BED (Q3's outside errors).
- **Descriptive.** Per-group `p`, `e`, coverage, the worst-case reweighted
  share, and `p` restricted to SNV pairs (the `indel` group's labels depend
  on Clair3 writing the indel the same way as dipcall; a missing record
  reads as homozygous reference).

## Risks stated now

- A child site that is not truly heterozygous (a paralog collapsed in the
  reference, often in `segdup`) has no true phase. The label then reflects
  how the parents' reads map, which may share the reference's error with
  the phaser. Such sites push `f` toward the phaser's side and lower `p`.
- Labels exist only where a parent is homozygous, so they reach a selected
  subset of disputes. The reweighted share bounds the effect of that
  selection by group, not within a group.
- Indel labels can fail by representation: a parent's indel written
  differently reads as absent. `a_in` measures this on HG002's indels.

## Script and order

`benchmarks/l3_outside_adjudication.py`, run once on the server after this
file is committed. Output
`benchmarks/results/runs/2026-11-22_l3_outside_adjudication.txt`.
Deviations are dated amendments below.

## Amendments

### Amendment 1, 22 November 2026, before any scoring

The interval in "Computation" assumed `f > 1/2`. The script takes, for each
end, whichever Wilson bound of `e` moves `p` outward: the lower bound of `p`
is the smaller of `(f_lo - e) / (1 - 2e)` at `e_lo` and at `e_hi`, the upper
bound the larger of the same at `f_hi`. For `f > 1/2` this is the stated
rule.

### Amendment 2, 22 November 2026: the run

The script ran once on the server (412 s), output
`2026-11-22_l3_outside_adjudication.txt`. L3-0 holds (a_in 0.968 WhatsHap,
0.955 LongPhase), L3-1 is supported (HG005 outside p 0.936 [0.927, 0.945]
WhatsHap, 0.642 [0.595, 0.686] LongPhase) and L3-2 is supported (HG002
outside 0.840 [0.814, 0.863], 0.717 [0.664, 0.765]). Labels reach 1.9 to
14.5 per cent of outside disputes. Seen after the run and not covered by any
rule: inside HG005's GIAB BED the labels side with the phaser against the
HPRC dipcall on many SNV disputes (WhatsHap p 0.317, 41 of 60), so part of
K3b's inside error is key error.
