# Literature re-check for the measurement paper, 10 October 2026

The paper's literature was last swept on 30 September and 1 October. This
note records a re-check for work that appeared since then. It changes no row
of `01-prior-art-matrix.md` and no claim in the earlier read-selection paper (development repository).

## What was searched

Europe PMC REST, first publication date 25 September to 31 December 2026,
which includes bioRxiv and medRxiv preprints:

| Query | Hits |
|---|---:|
| `"read selection" AND (phasing OR haplotype)` | 0 |
| `WhatsHap` | 6 |
| `"haplotype assembly"` | 4 |
| `"minimum error correction" AND haplotype` | 0 |
| `(phasing OR haplotype) AND "switch error"` | 0 |
| `"read-based phasing" OR "read based phasing"` | 1 |
| `(phasing OR haplotype) AND (calibration OR confidence) AND "long reads"` | 10 |
| `"registered report" AND (benchmark OR sequencing)` | 2 |
| `(preregistration OR preregistered OR "pre-registered") AND (benchmark OR benchmarking) AND (genomics OR bioinformatics)` | 2 |
| `(trio OR pedigree) AND phasing AND "long reads"` | 1 |

Web search: WhatsHap read selection against optimal selection in 2026, and
arXiv haplotype assembly or minimum error correction in October 2026. Neither
returned anything dated after 1 October. The arXiv API returned HTTP 429 and
was not retried.

## What was found

Every hit except one uses WhatsHap or phasing as a step in an applied study
(kidney allele-specific expression, *Anopheles* gene flow, NF1 somatic
mutations, repeat expansions, a cell-line assembly) or is off topic.

Heinz, Marin, Meyerson and Li, "HipHap: Haplotype Assignment and Confidence
Scoring for Diploid Reference Genomes", bioRxiv, posted 29 September 2026,
DOI 10.64898/2026.09.24.754153. Status: `abstract`. It aligns reads to each
haplotype of a diploid assembly separately, assigns each read to a
haplotype, and reports a per-read assignment quality (HapQ). It does not
phase variants from reads, select reads under a coverage cap, or calibrate a
per-region phasing error. It bears on no claim of the measurement paper. It
is a possible related-work row for the follow-up paper's confidence scores,
where it would be cited as per-read, not per-region, confidence.

These are searches, not proof of absence. No deliverable rests on them.

## HipHap's code, read 10 October

Repository `github.com/jheinz27/hiphap`, last push 1 October 2026. GitHub
reports no license file, so its code cannot be copied into this repository.
What the README and `src/lib.rs` show:

- HapQ is `6.02 * (AS_winner - AS_loser) / A`, capped at 60, where `A` is
  minimap2's match score; reads with more than three split alignments are
  scaled by `3 / n_splits`. 6.02 is about `10 log10 4`, so on a phred scale
  each base of score difference counts as a factor of four in odds. The source
  shows no fit to a truth, so HapQ is a heuristic phred scale.
- Exact ties get HapQ 0 and are assigned by a seeded hash of the read name,
  so tied reads are flagged and the choice is reproducible. Any other read
  gets at least 1.
- HapQ is turned off when comparing different references (GRCh38 against
  CHM13), where the README says it is not meaningful.

The tie flag is the same practice the measurement paper recommends for tied
sites (report abstentions beside errors). This is a README fact, not a
finding of the preprint, and is not cited.
