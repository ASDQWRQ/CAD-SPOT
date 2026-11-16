# Ovarian-primary / PDAC-transfer task evaluation

Specification fixed on 16 November 2026 before observing new ovarian CAD-SPOT
outcomes. This extends the 15 November benchmark; PDAC data and the pure
simulation outcomes were previously examined. It is not an untouched holdout.
No formulas, cutoffs, calipers or proxy panels are retuned on ovarian outcomes.

## Fixed accuracy comparisons

- Preserve seven scenarios and 40 paired seeds starting 881001: no nuisance,
  linear compact, linear dispersed, noisy proxy, nonlinear compact,
  architecture-coupled and no-program, on two 24-by-28 score-level grids.
- Preserve six methods: raw product, technical-only product, composition OLS
  product, adjusted target-only, adjusted minimum and additive spline product.
  All-spot AP is primary; top-5% precision/recall and outlier burden secondary.
  No-program AP remains undefined. All outcomes are reported.
- Real backgrounds: ovarian GSE211956 and PDAC GSE274103, first two sections in
  accession order, every third barcode-sorted spot from each. This is a fixed
  computational subsample, not selection by score or tissue phenotype.
- Compact/dispersed injections: 32 paired seeds starting 771001, ceil(5%) spots
  per section, target/reference increments 1.2/1.1 original section score SDs.
  Locations are chosen without score inspection. True labels are injected score
  perturbations, not validated biological regions.
- AP differences use 10,000 paired replicate resamples, seed 551001. They are
  conditional on the two backgrounds, not patient-population confidence bounds.

## Implementation checks

- Compare frozen v2.1.0 with v2.2.0 and an independent reduced-output OLS product.
  Version 2.2.0 retains v2.1.1 buffered ranking arithmetic; its new features are
  the Matrix Market ovarian adapter and section-pitch bootstrap geometry.
- Exact full ranking equality is required on both real tables and all three
  ranking modes on tie/missing-coordinate fixtures. The independent reference
  should match default priorities within 1e-10 and all candidate flags.
- Time five randomized serial subprocess repeats per build on ovarian 10,000
  spots (seed 661001 subset), full ovarian and full PDAC. Random job order uses
  seed 991001. Set OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=1. Run performance after
  other numerical jobs finish.
- Ranking-kernel time excludes downloads, preparation, bootstrap and matched
  sets. Also record process wall time and sampled aggregate process-tree RSS;
  shared pages may be counted more than once. The reduced-output reference does
  less work and is not an equivalent-function speed comparator.

## Reporting

No wet-lab validation, patient prediction, user-study result, or accuracy/speed
superiority over Giotto/Squidpy/latent spatial models is claimed. Geometry,
nonlinearity, proxy error and biological composition coupling define applicability.
The companion Study_Design.md documents cohort roles, QC and source reconstruction.
