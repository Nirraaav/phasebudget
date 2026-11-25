# The real-data truth set

## Why a third truth

The GIAB v4.2.1 benchmark VCF for HG002 is unphased. Two phased versions sit
in its `SupplementaryFiles/`, and neither is clean for this project:
`StrandSeqANDTrio` was phased partly with WhatsHap, which favours
WhatsHap-like methods when used to score them, and `hifiasm_v11_phasetransfer`
is described by Wagner et al. 2022 as "not evaluated for accuracy". OQ-12
records the choice. The third option is phasing by transmission from the
parents' own benchmark calls, which uses no reads and no phasing tool.

## Phasing by transmission

The inputs are the v4.2.1 benchmark VCF and BED of the child, HG002, and of
the father, HG003, and mother, HG004.

A child site is used when the child's call is a biallelic SNV (one-base REF,
one one-base ALT) with genotype `0/1`. Call its two bases `r` and `a`.

A parent's genotype at the site is a multiset of two bases:

* the bases of the parent's call, when the parent's VCF has a call at the same
  position with the same REF; a base outside `{r, a}` makes the site
  **inconsistent**;
* `{r, r}`, when the parent has no call and the position lies inside the
  parent's benchmark BED, because inside the BED a missing call is a
  confident homozygous reference;
* **unknown** otherwise, and the site is skipped. A parent call with a
  different REF is also skipped, as a representation difference.

Let `F` and `M` be the father's and mother's multisets. The candidate phasings
are `(p, m)` in `{(r, a), (a, r)}` with `p` in `F` and `m` in `M`, where `p`
is the paternal base. One candidate: the site is **phased**, and the output
genotype is paternal|maternal as allele indices, `0|1` or `1|0`. Two
candidates, which happens only when both parents are heterozygous:
**ambiguous**. None: **inconsistent**, a Mendelian violation or a de novo call.

Every phased site is oriented against the parents, so any two phased sites are
comparable, whatever their distance. There are no phase blocks.

## Comparing a phased file against it

A site counts when it is phased by transmission, and the other file's call at
that position has the same REF and ALT, a `|` genotype that is `0|1` or `1|0`,
and a phase set that asserts paternal|maternal order: `PATMAT` in the hifiasm
file, and `1` in `StrandSeqANDTrio`, whose README states that order. The
discordance rate is the fraction of counted sites whose paternal allele
differs. Local phase sets in the hifiasm file are not parent-oriented and are
reported separately, not counted.

## What this truth can and cannot do

It phases only sites where at least one parent is homozygous, so it misses
every site where both parents are heterozygous. Its errors are the parents'
genotyping errors, which GIAB's Mendelian filtering already reduces. It is
independent of every read-based phaser, which is the reason to use it.
Adjacent-pair switch error against it has to be computed over the sites it
phases, skipping the rest, as `whatshap compare` does within a phase set.

## Where the inputs come from

Downloaded on 25 September 2026 from the GIAB release tree,
`https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/`,
keeping the header and chromosome 20 of each file:

| local file in `data/giab/` | source under the release tree |
|---|---|
| `HG002_chr20_v4.2.1.vcf.gz` | `HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` |
| `HG003_chr20_v4.2.1.vcf.gz` | `HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` |
| `HG004_chr20_v4.2.1.vcf.gz` | `HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` |
| `HG00{2,3,4}_chr20_v4.2.1.bed` | the matching `..._benchmark_noinconsistent.bed` |
| `HG002_chr20_hifiasm_phasetransfer.vcf.gz` | `HG002_NA24385_son/NISTv4.2.1/GRCh38/SupplementaryFiles/HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.vcf.gz` |
| `HG002_chr20_StrandSeqANDTrio.vcf.gz` | `HG002_NA24385_son/NISTv4.2.1/GRCh38/SupplementaryFiles/HG002_GRCh38_1_22_v4.2.1_benchmark_phased_MHCassembly_StrandSeqANDTrio.vcf.gz` |

Each VCF was streamed through `zcat | awk '/^#/ || $1=="chr20"' | gzip`, and
each BED through `awk '$1=="chr20"'`. The chromosome-20 record counts were
85,951 for HG002 and both phased files, 84,252 for HG003 and 86,456 for
HG004; the BEDs have 10,192, 11,601 and 11,438 intervals. `data/` is not
tracked, so these two paragraphs are the record of what the numbers in
The development notebook were computed from.

## An error in the transmission truth, found on HG002 reads

At chr20:12,558,514 and 12,558,534 the transmission truth puts the child's
alternate alleles on the maternal haplotype, because HG004 is heterozygous
there and HG003 has no call inside its BED, which the rule reads as
homozygous reference. All 74 phasings of the HG002 run, from every pipeline
and both external tools, put them on the paternal haplotype, and so does
the hifiasm phase transfer. The father's own reads confirm the cause. In GIAB's HG003 PacBio HiFi
alignment (`HG003.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam`),
reads of mapping quality 20 or more show 41 T and 42 C at both positions,
so HG003 is heterozygous there, and his v4.2.1 benchmark has two false
negatives inside its own benchmark interval (chr20:12,556,741-12,558,536,
which ends 22 and 2 bp after the sites). Both parents are then heterozygous,
the site is ambiguous by transmission, and the truth's rule that a missing
call inside the BED is homozygous reference is what failed.
The two sites produce the two switch errors that nearly every phasing has.
The registered scoring keeps them; a post-hoc scoring that sets them to
unphased (`score_hg002.py --exclude 10793-10820`) is reported beside it.
