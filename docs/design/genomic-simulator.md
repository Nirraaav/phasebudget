# Simulating in genomic coordinates

## Why the variant-count simulator cannot vary heterozygosity

`phasebudget.simulate.simulate` takes `n_variants` and `read_span`, where
`read_span` counts variants. Its input is already a list of heterozygous sites,
so the distance between them does not exist as a quantity. Adding a
`heterozygosity` parameter to it would change nothing measurable: whatever the
density, a read still touches `read_span` variants, which is the only thing the
fragment matrix records. The development notebook recorded this as a stated limit of the
objective-truth grid, and it is why that grid varies coverage, span, and error
rate only.

## What density actually controls

A sequencing read has a physical length in bases. The number of heterozygous
variants it touches is a consequence of that length and of how densely the
variants are spaced:

```text
E[variants per read] = heterozygosity * read_length.
```

Heterozygosity is therefore not an independent axis alongside read span. It is
one of the two factors that produce read span, and the quantity a phasing
method actually feels is their product. A sweep over heterozygosity at fixed
physical read length is a sweep over phase linkage: halving the density halves
the expected variants per read, which lengthens the physical span of a phase
block while making each block rest on fewer observations.

## The parameterization

`simulate_genomic` fixes the number of variants and draws the gaps between
them, so a density sweep keeps the matrix width constant and changes the
genome length instead:

- variant positions start at `first_position` and each successive gap is drawn
  from a geometric distribution with mean `1 / heterozygosity`, which is the
  discrete waiting time between sites of a Bernoulli process at that rate.
- the genome spans `positions[0]` to `positions[-1]`, and read start points are
  drawn uniformly from `positions[0] - read_length + 1` to `positions[-1]`, so
  a variant near either end is as likely to be covered as one in the middle.
  Without the negative offset the first variants would be systematically
  under-covered.
- the read count is `round(coverage * span_in_bases / read_length)`, so
  `coverage` keeps its usual meaning of expected reads over a base.
- a read's fragment is the variants inside its half-open base interval. Reads
  containing no variant are dropped, and reads containing exactly one are kept:
  they carry an allele, they occupy coverage, and they contribute no phase
  link. `selection.is_phase_informative` is the filter for callers that want
  them gone, and it is deliberately not applied here.

The instance reports `positions`, so the base-pair forms of `block_n50`,
`an50`, and `qan50` in `phasebudget.metrics` mean what they say on these
instances. On `simulate` output those positions are a fixed grid with no
relation to the reads.

## What stays the same

Error rate, base-quality jitter, the interior gap rate, and the weight model
are the same parameters with the same meaning as in `simulate`, and both
generators return a `SimulatedInstance`. Coverage in `simulate_genomic` is
physical coverage over bases rather than over variants; at a density of one
variant per `1 / heterozygosity` bases the two agree in expectation once
`read_length * heterozygosity` variants fit in a read.

`simulate_genomic` is not a sequencing simulator. Reads have a fixed length
rather than a length distribution, error is independent across bases, variant
sites are known in advance and are all truly heterozygous, and there is no
reference bias, no mapping error, and no indel. It exists to make one axis
measurable, not to be realistic.

## Errors the base qualities do not report

Every error above is independent and reported: a call with error probability
`e` carries the weight `log((1-e)/e)`. That is the most favourable case for
every method here, and for read merging most of all, since two reads that
agree are then unlikely to share an error. Real long reads also make errors
their qualities do not describe. Two kinds are simulated, both off by default,
and when off the simulator draws no extra random numbers, so every earlier
instance is reproduced exactly; `tests/test_simulate.py` pins three of them
by fingerprint.

**Systematic sites.** A fraction `systematic_rate` of variant sites is chosen,
and each is given a biased allele. A read covering such a site reports the
biased allele with probability `systematic_strength`, whichever haplotype it
came from, and otherwise follows the ordinary error process. The weight is the
ordinary one, so the qualities are overconfident exactly there. This is the
signature of mapping artefacts and context-specific sequencing errors: reads
from both haplotypes agree on the same wrong allele, which makes two reads
from different haplotypes look alike to a merging rule and pulls MEC toward a
wrong column consensus.

**Bursts.** With probability `burst_rate` a read carries one burst: a
contiguous run of `burst_length` of its variants, placed uniformly, over which
the true error probability is `burst_error` while the weights stay those of
the ordinary process. Errors inside a burst are independent of each other but
all unreported, which is how a low-quality stretch or a local misalignment of
one read looks.

`SimulatedInstance.systematic` lists the systematic sites, so a benchmark can
tell errors the model could have known about from ones it could not.
