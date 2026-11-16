# Decision-utility extension fixed before execution

This extends the completed v2.2.0 ovarian-primary/retrospective-PDAC analysis.
It does not redefine the original scores, candidates, estimator or main cohorts.
All prespecified outcomes, including losses against simple baselines, are kept.
Software integration is a contribution to auditable use, not novel OLS algebra.

## Fixed-budget known-signal evaluation

Reuse every seed and generator from Utility_Benchmark/benchmark.py: seven
40-seed score-level scenarios and four 32-seed real-background injection cases.
Compare all six existing methods, without retuning. Within each section select
the first 10, 25 (primary) and 50 spots, breaking score ties by barcode without
truth labels. Report per-section precision/recall, positive-score selections,
and the number of inspections needed to recover 80% of known positives. The
last metric is undefined with no positives. Fixed-budget lists do not use the
descriptive outlier cutoff. Aggregate sections within replicate before paired
replicate inference; do not treat spots as independent replicates. Use 10,000
paired replicate resamples, seed 20261117, and retain all baseline comparisons.
Truth refers only to the specified injected perturbations, not real pathology.

## Model-dependence audit

For each full cohort, fit raw, technical OLS, composition OLS and composition
additive-spline products. Raw means within-section standardization without
regression; other fits include section indicators. Export per-section ranks
with deterministic barcode tie breaking. Audit-list membership also requires
positive priority, so it can contain fewer than K positions. Conditional-model agreement concerns
technical, composition and spline lists, not raw. At K=25 label spots selected
by all three, some, or none; retain membership at K=10/50 and flags showing
technical-versus-composition or OLS-versus-spline changes. This is descriptive
agreement, not FDR, biological accuracy or automatic model choice.

## Input perturbations

Use the first two accession-ordered sections of each cohort, all their retained
spots (not the every-third-spot injection subset). Fit new baseline models on
this same two-section universe. Twenty paired seeds (2026111700 through
2026111719) each evaluate independent binomial count thinning at 80% and 50%
retention, and 10% panel-member loss (ceil(10% of each used panel), at least
three remaining genes). Use the original fixed panel membership otherwise;
recompute normalization, gene z-scores, technical covariates and proxies. Keep
the original QC universe to isolate input perturbation; report the number now
below original QC thresholds, and fail if any library becomes zero or a fixed
gene becomes constant. Such a replicate is recorded as not computable, including
its reason and denominator; other scheduled replicates still run. Do not silently
change membership or exclude new rows. This clarifies failure handling, not the
predeclared scoring or perturbation settings.
Separately refit after omitting each of the five composition proxies. There is
no truth-based choice of the best proxy subset. Export top-K retention/Jaccard,
rank correlation and baseline audit-label retention for every section/method.
Stability is not correctness; baseline rank can influence stratum retention.

## Second, non-metabolic score pair

Freeze these illustrative, disjoint panels before examining associations:

- IFN-response target: IFIT1 IFIT2 IFIT3 ISG15 MX1 MX2 OAS1 OAS2 OAS3 RSAD2 IFI6 IFI44 IFI44L.
- Antigen-processing/presentation reference: B2M HLA-A HLA-B HLA-C TAP1 TAP2 TAPBP PSMB8 PSMB9 PSMB10 ERAP1 NLRC5.

Biological motivation: Reactome R-HSA-909733 (interferon alpha/beta signaling)
and R-HSA-983170 (class-I MHC folding/assembly/peptide loading), accessed
2026-11-16. These manually declared illustrative panels are not claimed to be
complete Reactome gene sets, independently validated signatures or activation
assays. They do not overlap each other, the original score pair, or proxies.
Use all 13 sections with each cohort's common nonconstant coverage, the existing
QC universe and scoring rule. Keep the original technical/five proxy covariates.
Run association/model-dependence audits regardless of direction or significance.
No second matched-null, biological accuracy or untouched external-validation
claim is planned. All missing/constant genes and exact weights are exported.

## Deliverables and constraints

Add a documented CLI producing a coordinate audit TSV, model summary, JSON
manifest and self-contained printable HTML decision report. No learned winner,
consensus superiority claim, or automatic suppression of unfavorable models.
Keep original v2.2.0 performance measurements identified as that version; new
audit functions are not included in those timings. Update manuscript, legends,
ordered tables, reviewer response, Word/LaTeX and archives together.
