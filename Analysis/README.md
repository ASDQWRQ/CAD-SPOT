# CAD-SPOT manuscript reproduction

Revision software version: **2.0.0-review**. Run commands from the project root.
The core workflow is score-agnostic. Optional expression preprocessing and
third-party model comparisons are separate from the core runtime.

## 1. Install and test the core

Use Python 3.12 in a new virtual environment. The tested environment is recorded
in each output manifest; the lower bounds in requirements are not a claim that
every possible dependency combination has been tested.

```bash
python -m venv .venv
# Activate .venv using the command appropriate for your operating system.
python -m pip install -r Analysis/requirements.txt
python -m unittest discover -s Software/tests -v
python Analysis/run_reanalysis.py --stages core
```

This reproduces the **frozen score-table** results using the supplied input:
23,013 spots; pooled raw rho 0.539512, technical rho 0.228645, composition rho
0.199254; 755 freshly calculated candidates. The supplied input is frozen, not
claimed to be reconstructed by downloading expression counts.

## 2. Spatial block uncertainty

```bash
python Analysis/run_reanalysis.py --stages bootstrap --reps 1000
```

Three widths (1,500/3,000/4,500 full-resolution image pixels) each use 1,000
section-stratified occupied-block draws. All models and numeric standardization
are refitted in every draw. The output intervals are conditional on the observed
sections, not patient-population confidence intervals. Blocks, valid-draw
counts and all draw-level correlations are released in `Source_Data/`.

## 3. Download and prepare the separate expression audit

```bash
python Analysis/download_data.py --datasets pdac annotation
python Analysis/prepare_pdac.py --root .
python Analysis/finalize_metadata.py --root . --gtf Data/Homo_sapiens.GRCh38.113.chr.gtf.gz
python Analysis/run_reanalysis.py --stages null
```

Public raw data are not duplicated in the software archive. The GEO processed
FFPE Visium matrices have 19/23 target and 20/21 reference genes. This audit
recomputes both observed and random scores from the same matrix, yielding an
observed adjusted rho approximately 0.324, not the frozen-score rho 0.199.
The 500 matched sets use fixed calipers and seed 20261115. Matching is not a
spatial-independence or spot-level significance test. Exon-union length does
not substitute for FFPE probe-target coverage or GC matching.

## 4. External technical-only audit

```bash
python Analysis/download_data.py --datasets external
python Analysis/run_reanalysis.py --stages external --n-external-sets 500
```

The ten selected human GSE203612 expression files supply library size and
detected-feature covariates, but not a composition table or registered spatial
coordinates. These samples are not validation of composition adjustment or
coordinate ranking. Each sample has its own gene-coverage list, null sets and
empirical correlation distribution. Source H5 checksums are recorded.

## 5. Limited model comparisons

See `Analysis/comparators/README.md` for official source commits, compatibility
patches, candidate-blind inputs and training budgets. These optional GPU-heavy
models require their own dependencies. A short smoke test is not reported as
a completed fit. POPARI's budget-limited final checkpoint is explicitly labeled;
no claim of optimization convergence or biological superiority is made.
Restore the pinned official sources and recorded patches with
`python Analysis/comparators/setup_sources.py` before optional model execution.

## 6. Figures and tables

```bash
python Analysis/make_figures.py
python Analysis/make_comparator_figure.py --summary Analysis/comparators/summary_all
python Analysis/build_tables.py
```

Each figure has corresponding TSV source data. Supplementary table numbering is
defined in the editable Excel `Contents` sheet. Candidate maps use GEO's
registered hires H&E images with the supplied scaling, not unscaled pixel
coordinates. They are localization views, not independent pathology validation.

## Files and estimators

- `Software/cadspot.py`: validation, OLS score-residual Spearman, genuine partial
  Spearman, additive-spline sensitivity, six-neighbor ranking, signed-product
  baseline comparisons and refitted spatial block bootstrap.
- `Software/cadspot_null.py`: strict gene-set generation and paired calibration.
- `Source_Data/archived_coordinate_table.tsv.gz`: preserved frozen input.
  Historical derived fields are renamed or retained as archival columns and
  do not override fresh outputs.
- `Source_Data/fresh_coordinates.tsv.gz`: fresh residuals, neighbors, components
  and candidate flags; the primary candidate column is `cadspot_candidate`.
- `Source_Data/fresh_run_manifest.json`: actual core arguments and hashes.
- `Source_Data/score_reconstruction_audit.tsv`: explicit input-score versus
  count-reconstructed score concordance; these are not identical estimators.

Use `python Software/cadspot.py --help` and
`python Software/cadspot_null.py --help` for parameters. Input and coordinate
conventions, missing-value behavior and all formulas are documented in
`Software/README.md`. Defaults are configurable; changing them creates a new
analysis rather than reproducing the reported manuscript run.
