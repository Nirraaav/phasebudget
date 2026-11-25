# Calibrating the per-boundary phase posterior

Written on 30 September 2026, before any correction was fitted, and before
any ONT phasing of a confirmatory chromosome had been compared with a truth.

## The problem

G8 of the genome-wide HiFi run found the posterior's per-boundary error
probability under-predicting errors: 23 observed against 8.1 predicted
under the registered truth, and 41 against 9.0 under v5.0q. G9 found it
ranking the wrong calls almost perfectly, with an AUC of 0.961 and 0.997.
Two causes are measured:
- *The allele weights are wrong.* The `realign` alleles carry WhatsHap's
  constant quality 30, an allele error of 0.001. On chromosome 20 the
  measured error is 0.00013 for HiFi and 0.0050 for ONT
  (`ont-prediction.md`).
- *Errors are not independent.* On HiFi chromosome 20, 16 per cent of the
  allele errors sit at 3 of 43,827 sites.

## The correction

For a posterior error probability `p` at a boundary:

```text
p' = eps + (1 - eps) * logistic(T * logit(p))
```

- `T` rescales the evidence, which absorbs a wrong allele weight;
- `eps` is a floor, the rate of errors that no amount of read evidence
  removes: site-level errors, truth errors, errors outside the vote.

Both are fitted by maximum likelihood on a calibration set of boundaries
with known outcome. A boundary counts when the selected reads span it and
v5.0q phases both of its variants; it is wrong when its called relation
differs from v5.0q's. The signals are those the cap-5 `optimal` runs
record (`hg002_pipelines.py --signals`), on `realign` fragments with
log-odds weights.

## Evaluation

On a test set, as in G8:
- the ratio of observed to predicted wrong calls, with a 95 per cent
  interval from 2,000 resamples of 1 Mb windows (seed 0);
- the expected calibration error over ten equal-width bins.

Because nearly every `p'` is below 0.1, the table also gives observed
against predicted in bins of `log10 p'`, and the log-loss. **Calibrated**
means the ratio's interval contains 1 and the expected calibration error is
below 0.01, as in G8. A correction that is calibrated only in the pooled
ratio, and not in the log bins, is reported as such.

## HiFi, exploratory

Fit on the odd confirmatory chromosomes and test on the even ones, at full
depth. Both were seen in G8, so this supports no claim.

## ONT, registered

Fit on chromosome 20 of the ONT run, which enters no test. Test on the 21
confirmatory chromosomes at full depth, and report the other depths beside
them. The correction is supported if it is calibrated at full depth by the
rule above. The uncorrected posterior is reported beside it.

## Fitted values, frozen on 30 September 2026

- **HiFi** (odd chromosomes, full depth; exploratory): `T = 0.7464`,
  `eps = 4.255e-06`. On the even chromosomes the corrected ratio is 0.71
  [0.40, 1.06], against 2.09 [1.21, 3.29] uncorrected, so it passes the rule
  at full depth. The log bins show where the rule is lenient:
  - 8 of the 31 errors sit where the posterior gives below 1e-7;
  - the bin 1e-3 to 1e-2 is over-predicted, 30.4 against 12.
  At 10-fold the uncorrected posterior is already calibrated (1.06), and the
  full-depth correction over-predicts (0.58). Output:
  `benchmarks/results/runs/2026-09-30_calibration_hifi.txt`.
- **ONT** (chromosome 20, full depth): `T = 0.7264`, `eps = 6.9e-287`,
  effectively zero. The fit rests on 5 wrong calls in 44,254 boundaries. The
  registered test uses these values unchanged:
  `PHASEBUDGET_PLATFORM=ont ... calibrate_posterior.py --T 0.7264 --eps 0
  --test <21 chromosomes> --depths full 30 20 10`.
