# CAD-SPOT 2.2.0

CAD-SPOT takes user-supplied target/reference scores and covariates, reports
association layers, and exports an interpretable coordinate ranking. This local
release defaults to **positive residual co-enrichment**. Six-neighbor reference
support remains visible as an annotation, not a required positive filter for the
default ranking. Optional support-weighted methods remain available.

The count adapters are [prepare_ovarian.py](../Analysis/prepare_ovarian.py) and
[end_to_end.py](../Analysis/end_to_end.py).
Its prepared primary input is [Primary/score_table.tsv.gz](../Primary/score_table.tsv.gz).
Numerical analysis results belong in [Primary/Results](../Primary/Results/), not
in hardcoded example values in this README.

## Version and Installation

The designated Python 3.12 runtime passes 28 core unit/integration tests, with
additional count-adapter and pitch-scaling tests under Analysis. Version 2.2.0
retains the v2.1.1 core ranking/model arithmetic and adds upstream adapters and
explicit bootstrap-geometry transformation. The local
version is `2.2.0`, exported as `cadspot.__version__` and the manifest's
`cadspot_version`. Identify the build by this version, source hash and selected
ranking method. Source hashes are recorded in the release manifest; the exact
release is tagged `v2.2.0` in the public repository.

Version 2.1.1 buffers six-neighbor output writes per section instead of assigning
pandas fields per spot. Arithmetic, tie handling and candidate rules are
unchanged. Full outputs are exactly equal to 2.1.0 on both real score tables and
edge-case fixtures. Utility_Benchmark contains all prespecified accuracy/timing
comparisons, including cases favoring baselines; it is not a universal accuracy
or speed claim over other spatial tools.

Version 2.2.0 retains that ranking arithmetic. It adds Matrix Market ovarian
reconstruction, coverage-only cohort panel harmonization and section-pitch
bootstrap geometry, while preserving image pixels in coordinate exports.

Run the examples below from the **package root**, one level above Software:

```bash
python -m pip install -r Software/requirements.txt
python -B -m unittest discover -s Software/tests -v
```

The analysis driver additionally requires `h5py` for 10x H5 count matrices.
It is not needed by the table-only core. Run the preparation tests with
`python -B -m unittest discover -s Analysis -p 'test_*.py' -v`.

## Primary Example

The primary ovarian example uses fresh scores from public GSE211956 counts and registered
spatial metadata, rather than an archived target/reference table. The current
panel definitions and QC decisions are in
[Primary/pipeline_config.json](../Primary/pipeline_config.json); public input
locations and checksums are in [Primary/input_manifest.json](../Primary/input_manifest.json).
Executable reconstruction and analysis commands are in
[Analysis/README.md](../Analysis/README.md).

The preparation driver aggregates duplicate gene symbols, retains in-tissue
spots with at least 100 UMIs, at least 100 detected unique symbols and finite
coordinates, then library-normalizes to 10,000 counts and applies log1p.
Each available nonconstant panel gene is standardized within section using
population SD; its panel score is the unweighted mean. Missing/constant genes
are recorded and omitted, not zero-imputed. The driver requires at least three
usable genes in each panel. Section folders retain per-gene scoring parameters,
panel coverage, QC records and matrices.

The ovarian adapter fixes cohort membership using the intersection of available
nonconstant genes before any association fitting. Noncommon exclusions and
revised equal weights are exported. This harmonizes observed/random panel
definitions; it is not selection by correlation or candidate outcome.

The five disjoint composition-proxy panels are:

| Column | Intended genes |
| --- | --- |
| `Epithelial_proxy` | EPCAM KRT7 KRT8 KRT18 KRT19 |
| `ECM_proxy` | COL1A1 COL1A2 COL3A1 DCN LUM |
| `Contractile_proxy` | ACTA2 TAGLN FAP PDGFRB |
| `Myeloid_proxy` | LST1 TYROBP FCER1G AIF1 CD68 |
| `Tcell_proxy` | CD3D CD3E TRAC CD247 LCK |

These are transcriptional proxies, not malignant-cell calls, validated cell
fractions or functional activation measures. They exclude target/reference
genes. Actual used membership is recorded in
[Primary/panel_coverage.tsv](../Primary/panel_coverage.tsv), not inferred from
the intended lists. The generic core is score-agnostic; the panels are decisions
of this example, not mandatory software inputs for other applications.

Run the prepared primary table through the CLI:

```bash
python -B Software/cadspot.py --input Primary/score_table.tsv.gz --target target_score --reference reference_score --technical log_nCount log_nFeature --composition Epithelial_proxy ECM_proxy Contractile_proxy Myeloid_proxy Tcell_proxy --coordinate-units full_resolution_pixels --ranking-method positive_product --output Primary/Results/cli_associations.tsv --coordinate-output Primary/Results/cli_coordinates.tsv.gz --ranking-comparison-output Primary/Results/cli_ranking_comparison.tsv
```

The `cli_` output prefix avoids overwriting the analysis driver's outputs. Omitting
`--ranking-method` also selects `positive_product`.

## Input Contract

The core reads tab-delimited tables, optionally gzip-compressed. A minimal
primary-example header is:

```text
sample_id	barcode	x	y	target_score	reference_score	log_nCount	log_nFeature	Epithelial_proxy	ECM_proxy	Contractile_proxy	Myeloid_proxy	Tcell_proxy
```

- `(sample_id, barcode)` identifies one spot uniquely. Repeated barcodes in
  different sections are allowed; missing/empty IDs and within-section duplicates
  fail. Column names are configurable.
- Target/reference scores and covariates must be numeric, finite and complete.
  There is no implicit median imputation or silent target-row deletion.
- Technical covariates in this example are `log1p(total UMI)` and
  `log1p(detected unique symbols)`. Composition covariates are optional in the
  generic CLI. Targets/reference cannot also be covariates.
- Coordinates are physical x/y in common Euclidean units within section; this
  example uses full-resolution pixels, x=column and y=row. `--coordinate-units`
  records a declaration, not a unit inference. Never mix pixel and array axes.
- A shared bootstrap width assumes comparable coordinate scales. Scale choices
  across studies must be declared; the analysis driver exposes width settings.

## Association Layers

`associations.tsv` contains pooled and separate-section results. `raw` is ordinary
Spearman without section/covariate adjustment. Adjusted layers are
`technical_adjusted` and, when proxies are supplied, `composition_adjusted`.
Pooled adjusted models include section fixed effects; separate-section fits use
only that section's spots.

`--methods` defaults to these three estimators:

- `score_residual_spearman`: regress each unranked score on standardized numeric
  covariates, intercept and section dummies, then correlate residuals by Spearman.
- `partial_spearman`: average-rank both scores and numeric covariates, retain
  unranked section dummies, fit OLS, then Pearson-correlate rank residuals.
- `spline_score_residual_spearman`: additive cubic B-splines with quantile knots,
  intercept and section dummies, then Spearman of score residuals. Duplicate
  knots are collapsed; constant covariates are omitted.

Undefined constant/too-small correlations have an explicit status, not an
invented zero. These layers describe covariation under supplied proxies; they
do not establish causal adjustment or biologically truer residual signals.

## Ranking and Candidate Rule

Ranking uses the final **linear score residuals**, independently of sensitivity
estimator choices. Residuals are standardized within section using population SD.
Write `a = target_resid_z` and `b = reference_resid_z`. The default score is:

```text
P = max(a, 0) * max(b, 0)
```

For each section, center and scale use all spatial-valid spots, not just the
positive quadrant. The exported metadata make the actual rule explicit:

| Column | Value |
| --- | --- |
| `priority_center` | section median of `P` |
| `priority_scale` | `1.4826 * MAD(P)` if positive; otherwise population SD (`ddof=0`) if positive; otherwise 1 |
| `priority_scale_method` | `MAD`, `SD`, or `constant`, identifying the chosen scale |

The existing `priority_robust_z` column is calculated as
`(P - priority_center) / priority_scale`. The positive-part product can contain
many exact zeros, so MAD can be zero even for a nonconstant score. The default
therefore may use the **SD fallback**, rather than MAD scaling; the retained
column name does not guarantee a MAD-based score. Missing/excluded spatial rows
have missing center/scale metadata, not an invented `constant` label. A default
candidate requires all of:

- `spatial_status == "ok"`;
- `a > 0` and `b > 0`;
- `priority_robust_z >= 3`, adjustable with `--robust-z-cutoff`.

There is **no** required `h > 0` filter for this method. Zero/negative neighborhood
support does not by itself disqualify a positive-product candidate. This is a
heuristic outlier rule, **not a Gaussian three-sigma rule**, calibrated spot-level
significance or FDR control. No Gaussian distribution or tail probability is
assumed by the cutoff, including when SD is used. It does not fix a candidate
percentage or desired count. These metadata annotate the existing computation;
they do not change numerical scores or candidate selection.

`priority_rank` ranks every eligible spot within section, descending with minimum
rank for ties, before applying the robust-z cutoff. The coordinate table records
the selected `ranking_method`; CLI manifests record `arguments.ranking_method`
even when the method is selected by default. Unknown API methods raise
`ValueError`; unknown CLI methods fail argparse validation with exit code 2.

## Six-Neighbor Annotation and Eligibility Fields

Each spatial-valid spot has exactly six other Euclidean nearest neighbors in
its own section. Distance ties are resolved deterministically by spot ID; the
focal row is excluded even when positions coincide. Reference residual z values
are averaged over these neighbors and those means standardized within the
section's spatial-valid subset. This annotation is `h`, exported as
`neighbor_reference_resid_z`. `local_support = a*h` is also exported. Neither
term weights the default product.

The masks are deliberately separate:

| Exported mask | Meaning |
| --- | --- |
| `positive_product_eligible` | spatial-valid, `a > 0`, `b > 0` |
| `positive_enrichment_eligible` | historical support mask: product eligibility and `h > 0` |
| `ranking_eligible` | eligibility of the selected method |

The support mask is **not** reused as the default product candidate mask.
`ranking_eligible` is used for selected-method ranks and candidates.

Missing-coordinate spots retain nonspatial residuals/products but have
`excluded_missing_coordinates`, no spatial score/rank and candidate false.
Sections with fewer than seven coordinate-valid spots have
`excluded_insufficient_neighbors`. These existing spatial-validity rules apply
to all methods so six-neighbor annotations are available; missing support is
never silently converted to zero.

## End-to-End Driver and Results

Use [Analysis/end_to_end.py](../Analysis/end_to_end.py) for preparation and the
two adjustment-layer coordinate outputs, associations, block bootstrap and
matched-set audit. Raw input provisioning is separate: `prepare-primary` expects
the public Matrix Market/H5 files and extracted spatial metadata at the layout its
`--source-root` specifies. Their actual public URLs and checksums are documented
in the input manifest linked above; the core does not download or guess inputs.

For already prepared inputs, the following runs associations, both coordinate
layers and spatial bootstrap without requesting matched-set calibration:

```bash
python -B Analysis/end_to_end.py --stage run --folder Primary --reps 1000 --bootstrap-in-pitches --seed 20261116 --n-sets 0
```

The integrated release includes the annotation table at
`Provenance/ensembl113_exon_union_lengths.tsv`. Its public-GTF derivation is
implemented by `Analysis/build_gene_lengths.py`. A legacy `--source-root` may
also provide this table under `Source_Data/` when the release-local file is
absent. Do not
mistake the `--n-sets 0` command for a matched-null analysis. A full driver run
can overwrite its corresponding Results files; archive completed runs first.

The authoritative primary output directory is [Primary/Results](../Primary/Results/).
Expected driver outputs include:

| Output | Content |
| --- | --- |
| `associations.tsv` | all association layers and sensitivity estimators |
| `technical_coordinates.tsv.gz`, `composition_coordinates.tsv.gz` | separate coordinate rankings with selected-method metadata |
| `candidate_counts.tsv` | section and adjustment-layer counts |
| `*_baseline_comparison.tsv` | descriptive ranking agreements on common universes |
| `bootstrap_*_ci.tsv`, `bootstrap_*_draws.tsv.gz`, `bootstrap_*_audit.tsv.gz` | intervals, all draws and sampling audit |
| `blocks_*.tsv.gz` | section-specific block membership |
| `matched_summary.tsv`, `matched_rhos.tsv`, `matched_assignments.tsv.gz` | matched-set calibration when requested |
| `matching_metadata.tsv.gz`, `matching_universe.tsv` | matching inputs when requested |
| `run_manifest.json` | driver settings, source table hash and runtime versions after successful completion |

Results may appear incrementally while the driver is running. A partial output
directory does not establish successful completion or validate later stages.
Do not substitute statistics or candidate counts from the previous revision.

## Spatial Bootstrap and Matched Nulls

The core CLI exposes `--bootstrap-output`, explicit `--block-size`,
`--bootstrap-reps`, seed and auxiliary draw/audit/block-map outputs. Blocks are
non-overlapping squares with section-specific coordinate origins. Each draw
resamples occupied blocks with replacement within each fixed section and refits
the score models, numeric scaling and section effects. Requested partial-rank or
spline sensitivities are also refitted, not applied to cached residuals.
Intervals are conditional on the observed sections and chosen spatial scheme,
not patient-population inference. Missing-coordinate exclusions and undefined
draws remain explicit; an entire coordinate-missing section fails.

`cadspot_null.py generate` requires explicit expression/spot axes, metadata and
a gene universe. `calibrate` rebuilds observed and random scores with the same
scoring/estimator choices. Use `--help` on either subcommand for all arguments.
The analysis driver selects `mean_gene_z_by_section` to match primary scoring.
Its strict calipers are 0.25 on `log1p(mean_expression)`, 0.05 on detection rate
and 0.25 on log gene length. Target/reference and proxy genes are excluded;
unique-gene matching fails if infeasible and does not relax tolerances.

The one-sided plus-one value is `(1 + count(null >= observed)) / (1 + n_controls)`.
This is a gene-set audit, not a spatial-dependence test or coordinate-level FDR.
GC/probe-coverage matching is not implemented. Available gene membership and
scoring must match the primary model; the driver checks the observed statistic
against its main analysis rather than calibrating a different score.

## Migration and Optional Neighborhood Methods

Earlier scripts that omitted `ranking_method` used the four-term heuristic.
Replaying it now requires explicit selection:

```python
ranked = cadspot.spatial_ranking(table, target_residuals, reference_residuals,
                                 ranking_method="legacy")
```

| Method | Priority score | Eligibility |
| --- | --- | --- |
| `positive_product` (new default) | `max(a,0) * max(b,0)` | spatial-valid, positive `a`/`b` |
| `legacy` | `(a + b + h + a*h) / 4` | spatial-valid, positive `a`/`b`/`h` |
| `positive_support` | `max(a,0) * max(b,0) * max(h,0)` | spatial-valid, positive `a`/`b`/`h` |

The optional methods preserve their prior numerical rules and support filter.
Neither older counts nor archived score correlations are this follow-up's primary
example. Saved derived fields are recomputed, with protected collisions renamed
`archived__...`; saved eligibility never overrides the newly selected method.
The protected archive list also includes `priority_center`, `priority_scale`
and `priority_scale_method`, so old metadata cannot supply fresh scale values.

The completed [ranking ablation summary](../Ranking_Ablation/summary.tsv) shows
that relative performance depends on signal geometry: the legacy method ranks
better in the compact setting, the simple product in the dispersed setting.
Both select candidates in the no-program simulation. These experiments support
retaining neighborhood weighting as an optional sensitivity, not a universally
better default. They do not establish biological accuracy or FDR control.

The signed baseline `a*b` remains exported; negative-negative pairs can score
high and are not positive co-enrichment. The positive baseline remains exported
as `baseline_positive_product`; it equals the new priority on spatial-valid
rows by construction, not independent validation. `ranking_comparison` uses the
same spatial-valid universe for both scores, with separately declared all-spot,
positive-quadrant and support-positive restrictions. Top fractions use that
common universe and ties use deterministic IDs. Agreement is descriptive.

## Tests and Scope

The tests cover default/explicit product equivalence; zero/negative `h`
eligibility; legacy numerical regression; optional support filtering; inclusive
cutoff and ranking rules; MAD/SD/constant scale metadata and their archival;
spatial exclusions; API/CLI
unknown-method rejection; method metadata; archive separation; and unchanged
associations across ranking methods. Existing tests also cover partial/spline
estimators, bootstrap refitting, null generation/calibration and sparse keyed I/O.
Temporary test outputs stay inside Software and are removed afterward.

Unit fixtures check implementation contracts, not superiority on simulations or
real data. Ovarian preparation and PDAC transfer are documented in
`Analysis/README.md` and `Study_Design.md`. Successful execution and unit tests do not
establish histopathologic or biological validation.
