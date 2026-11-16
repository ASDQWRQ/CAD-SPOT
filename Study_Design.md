# Ovarian primary example and PDAC transfer assessment

Analysis specification fixed on 16 November 2026, before inspecting the new
ovarian CAD-SPOT outcomes. This is a retrospective extension, not prospective
registration or a claim that the existing PDAC data were untouched.

## Roles and fixed choices

- Primary worked example: all eight public Visium samples GSM6506110-GSM6506117
  from GEO GSE211956. Retrieve missing public count/spatial files rather than use
  historical image-prediction targets, deconvolution scores or selected barcodes.
- Cross-cancer transfer: the five GSE274103 PDAC sections already analyzed in the
  previous revision. Their prior use is disclosed. Do not call them an untouched
  external validation set or claim that CAD-SPOT was developed exclusively on
  ovarian data.
- Keep the v2.1.1 mathematical estimators, positive-part product, optional
  six-neighbor annotation, QC (100 UMIs/100 detected symbols/in-tissue/finite
  coordinates), intended 23/21 gene panels and five disjoint RNA proxy panels.
  Do not optimize panels, thresholds or covariates on the new outcomes.
- Candidate threshold remains an adjustable median-centered outlier statistic
  >=3 with actual MAD/SD fallback exported; it is not coordinate-level FDR.
- Reconstruct counts to log1p CP10K, within-section gene z scores (ddof=0), then
  equal-weight panel averages. Before fitting any ovarian association, use the
  intersection of nonconstant available panel genes across its sections, with
  all exclusions and revised weights exported. This input-coverage rule avoids
  comparing a pooled observed score with differently defined random scores;
  membership is never selected by association strength or candidate outcomes.
- Preserve full-resolution image x/y in coordinate outputs. For resampling only,
  divide coordinates by each section's median nearest-other-spot pitch; use
  block widths 6, 12 (primary), 18 pitches, 1,000 draws/width, seed 20261116.
  Export section-specific transformations. Scores/proxies stay fixed in draws;
  residual models/design scaling are refitted. No patient-population inference.
- Generate 500 strict matched sets per cohort, using the existing calipers,
  Ensembl-113 exon-union lengths, exclusions and seed 20261116. Fail explicitly
  when full distinct-gene matching is impossible. Report such failures rather
  than relax calipers or swap favorable target genes.
- Compare cohorts using their declared available-gene scores and a post hoc
  common-used-gene sensitivity. Do not interpret cohort differences as causal
  cancer-type effects; assay, treatment and sampling differ.

## Task evaluation and reporting

- Preserve all seven existing 40-seed score-level simulation scenarios, their
  six fixed methods and all favorable/unfavorable outcomes.
- New real-background injection comparison: first two barcode-sorted sections
  per cohort, every third spot, selected without score inspection. Thirty-two
  paired seeds starting 771001 per compact/dispersed geometry; ceil(5%) positives
  per section; amplitudes 1.2/1.1 original score SDs. This evaluates score-injection
  recovery, not biological ground truth. Paired AP intervals use 10,000 draws.
- Repeat implementation-equivalence and randomized single-thread performance
  checks on ovarian, ovarian 10,000-spot subset (if sufficient) and full PDAC
  inputs. No accuracy superiority over mathematically equivalent OLS is expected.
- Reuse PDAC spatial-model fits only as supplementary context with their fixed
  scope/convergence limitations; do not imply an ovarian model comparison.
- Main narrative: a reusable workflow illustrated thoroughly in ovarian cancer
  and assessed for cross-cancer technical transfer to PDAC. No prognostic,
  mitochondrial activation or independently validated candidate-biology claim.
- Rebuild all figures, captions, editable tables, manuscript/response/availability
  text, Word, LaTeX and final archives from these identified results.

## Provenance

The incoming ovarian raw-data directory is D:/hgsoc_visium_clean_raw. It contains
additional ovarian studies, but they are not substituted based on outcomes.
GSE211956 spatial inputs missing there are retrieved from NCBI GEO into the
project Data cache. E: is unavailable in this execution environment; original
files and previous submission archives are preserved.
