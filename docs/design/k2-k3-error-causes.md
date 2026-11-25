# K2 and K3: where phasing errors come from, and the benchmark gap on HG001

Written 19 October 2026, before the scorer ran. Development gates, not
registered. Directions K2 and K3 of the October directions memo (development repository),
round 7.

## K2. Search, data, or no evidence

Q3's cause split (`docs/design/q3-outside-benchmark.md`, section "Causes") is
applied to the consecutive phased pairs of every stored 30x run whose truth
phases both variants. For each wrong pair, a window of 10 variants on each
side is cut from the run's fragments, reads are chosen by the library's port
of WhatsHap 2.8's selection at cap 15, and the exact weighted MEC optimum is
computed with the relation forced each way. A wrong pair is `none` (no read
calls both variants), `search` (the window optimum strictly prefers the true
relation), `data` (most spanning reads back the wrong relation or tie, and the
optimum does not prefer the truth) or `other` (most spanning reads back the
truth, the optimum does not).

**Correction to the directions entry.** The entry said Q3 measured the split
outside the benchmark only. Q3's output also has it inside, on Q100's
variants scored against Q100: whatshap's `data` share inside is 0.618
(`segdup`, 330 wrong), 0.583 (`lowmap`, 60), 0.694 (`tandem`, 1,068) and
0.659 (`indel`, 11,160). This was read before this note was written, so the
HiFi half of the gate below is expected to fail in `lowmap`. The gate is
kept as the directions entry stated it.

**Gate (as stated on 19 October).** On HG002 inside the v5.0q BED, HiFi
(`w4_grch38`, the v5.0q variants and truth) and ONT (`cal_v5`), the `data`
share of whatshap's and LongPhase's wrong pairs is at least 0.6 in each of
`segdup`, `lowmap`, `tandem` and `indel` holding at least 20 wrong pairs.

**Descriptive.** The same split for every other stored 30x run: `cal_r10`,
`hprc_HG002`, `hg001_v5`, the five `ped_*` and five `phifi_*` Pedigree
members, and the 39 scored W1 genomes against their consensus truth. The
claim the paper can use, if the numbers allow it, is the `search` share: the
part a better solver could remove.

## K3. The benchmark gap on a second genome

Q3 found on HG002 that relations outside v5.0q's BED hold half of whatshap's
and 94 per cent of LongPhase's wrong relations. K3 asks the same of HG001
(NA12878) ONT sup v5.0.0 (`hg001_v5`, 30x), phased on the Platinum Pedigree
truthset's heterozygous variants inside its high-quality regions and scored
against that truthset. Inside means both variants lie in GIAB's HG001 v4.2.1
benchmark BED (`data/giab_hg001/HG001_GRCh38_1_22_v4.2.1_benchmark.bed`);
outside means at least one does not, though both are inside the Pedigree
regions.

**Change from the directions entry.** The entry asked for an assembly truth.
The Pedigree truthset is inheritance-checked across four generations and does
not come from either phaser, which is the property the gate needs; no
assembly-based truth for a second genome is on the server. Stated before any
count was made.

**Gate.** For both whatshap and LongPhase, the share of wrong pairs outside
the GIAB BED is at least 0.3. Reported with it: the error per pair inside and
outside and their ratio, per group, and the cause split on each side.

The HG001 data were seen in X8 and the round 6 checks, but no inside and
outside split against the GIAB BED was computed.

## Code

`benchmarks/k2_k3_error_causes.py`, reusing `q3_outside_benchmark.py`'s
tested window, cost and cause functions (`tests/test_q3_outside_benchmark.py`).
