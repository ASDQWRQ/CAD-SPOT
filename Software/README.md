# CAD-SPOT Score-Table Revision

This directory is a self-contained revision implementation. It does not modify
or import the old software package. Python, NumPy, pandas, SciPy, and
scikit-learn are required. Install the core dependencies using
`python -m pip install -r requirements.txt`. Optional comparison and plotting
dependencies are documented separately in `../Analysis/comparators/README.md`.
The verified core runtime uses Python 3.12.14.
Every CLI run writes a JSON manifest with arguments, input hashes, and exact
runtime versions. All random generation uses NumPy's recorded seed.

## Portable Quick Start

From this directory, run the released score table without installing a spatial
object framework:

```bash
python cadspot.py --input ../Source_Data/archived_coordinate_table.tsv.gz --target target_score --reference reference_score --technical log_nCount log_nFeature --composition Epithelial_malignant_proxy ECM_Stroma CAF Macrophage T_cell --coordinate-units full_resolution_pixels --output associations.tsv --coordinate-output coordinates.tsv.gz --ranking-comparison-output ranking_comparison.tsv
python -m unittest discover -s tests -v
```

`../Analysis/run_reanalysis.py` contains the complete manuscript analysis
commands, including the three block widths and the fresh gene-set calibration.

## Two Separate Analyses

1. **Score-table demonstration:** the supplied `target_score` and
   `reference_score`, two technical columns, five composition columns, and
   section fixed effects reproduce pooled linear score-residual Spearman
   `0.1992540732014225` (verified on 23,013 supplied rows). This reproduces a
   table-level association, not the upstream construction of its scores.
2. **Fresh expression null:** observed and random gene-set scores are both
   recomputed from the supplied expression matrix with the same scoring and
   residual model. The explicit available sets can differ from the intended
   sets. This is not a reproduction of the archived scores or their old null.
   Never compare the old null's observed statistic to `.199` as though they
   were the same model. Fresh observed rho is exported, not forced to `.199`.

The supplied intended target/reference lists have 23/21 genes. The explicit
available lists supplied for this revision have 19/20 genes: target ALDOA,
TPI1, PKM, LDHA and reference GAPDH are omitted, not zero-imputed. CLI validation
checks used genes against the matrix; intended-list audits report missing and
present-but-unused genes. These list files record the supplied specification,
not an independent biological validation of those signatures.

## Score-Table CLI

Run from this directory; replace output paths and block size as appropriate:

```powershell
& 'D:\hgsoc_cuda_env\Scripts\python.exe' -B cadspot.py `
  --input 'D:\CADSPOT_REVISION_20261115\Source_Data\archived_coordinate_table.tsv.gz' `
  --target target_score --reference reference_score `
  --technical log_nCount log_nFeature `
  --composition Epithelial_malignant_proxy ECM_Stroma CAF Macrophage T_cell `
  --sample-col sample_id --id-col barcode --x-col x --y-col y `
  --coordinate-units pixels --ranking-method legacy --robust-z-cutoff 3 `
  --top-k 100 --top-fractions 0.01 0.05 0.10 `
  --output example_output/associations.tsv `
  --coordinate-output example_output/recomputed_coordinates.tsv `
  --ranking-comparison-output example_output/baseline_comparison.tsv `
  --archived-comparison-output example_output/archived_candidate_comparison.tsv
```

Input is a tab-delimited table, optionally gzip-compressed. `(sample_id,
barcode)` must uniquely identify each row; repeated barcodes across sections
are allowed. IDs must be nonempty. Scores and covariates must be numeric,
finite, and complete. This implementation does not median-impute covariates
or silently drop missing targets. Prepare and archive an explicit complete
analysis table if necessary. The target and reference cannot be covariates.

`--methods` defaults to all three adjusted estimators below. `raw` is always
reported once as ordinary Spearman with no covariate/section adjustment.
`associations.tsv` includes pooled and individual-section rows:

- `score_residual_spearman`: OLS of each **unranked score** on standardized
  numeric covariates, intercept, and section dummies; then Spearman of the
  two residual vectors. This is the primary `.199` statistic.
- `partial_spearman`: average-rank both scores and each numeric covariate,
  retain unranked section dummies, fit OLS, then **Pearson** correlation of
  rank residuals. It is not Spearman of score residuals.
- `spline_score_residual_spearman`: additive cubic B-spline covariate terms
  with `--spline-knots 4` quantile knots, intercept, and section dummies;
  then Spearman of score residuals. Duplicate quantile knots are collapsed;
  constant covariates are omitted. There are no interactions or penalties.

Layers are `raw`, `technical_adjusted`, and (when supplied)
`composition_adjusted`. Per-section fits use that section's observations;
pooled fits include fixed effects for the observed sections. Undefined
constant/too-small correlations have an explicit status, not an invented zero.
No nominal independent-spot p value is supplied.

## Spatial Ranking and Archive Separation

Coordinates must be finite or explicitly missing, with x and y in the same
Euclidean units **within each section**. `--coordinate-units` records the
user's declaration, not a unit inference. Never combine a pixel x axis with
an array-index y axis. If scales differ across sections, harmonize them before
using one shared block size. Coordinate values alone cannot verify units.

The ranking uses the final linear score model (technical plus composition,
or technical only), independent of sensitivity estimator choices:

1. Standardize target/reference residuals within section (population SD).
2. For each spatial-valid spot, find exactly **six other Euclidean nearest
   neighbors within that section**, with no radius-based neighborhood or
   cross-section pooling. Seven valid spots suffice. Distance ties use stable
   spot-ID ordering; duplicate positions still exclude the focal row itself.
3. Average reference residual z over those six neighbors, then standardize
   the neighbor means within the section's spatial-valid spots.
4. Set `local_support = target_resid_z * neighbor_reference_resid_z`.

`--ranking-method legacy` is the explicit default existing heuristic:
`priority_score = mean(target_z, reference_z, neighbor_reference_z,
local_support)`. `positive_support` is an optional sensitivity:
`max(target_z,0) * max(reference_z,0) * max(neighbor_reference_z,0)`.
Within-section robust priority z uses median and `1.4826 * MAD`, falling back
to population SD, then 1 for a constant priority. Candidates require robust
z >= cutoff and **all three** z terms > 0. Candidate `priority_rank` is
restricted to this positive-support eligibility (without the cutoff).

Missing-coordinate spots retain score residuals but have
`spatial_status=excluded_missing_coordinates`, no spatial priority/rank,
neighbor count 0, and candidate false. A section with fewer than seven
coordinate-valid spots gets `excluded_insufficient_neighbors`; other sections
continue. Missing coordinates never silently become zero neighborhood support.

Fresh values are always recomputed. Colliding saved output columns, including
`priority_score` and `cadspot_candidate`, are renamed `archived__...`. Other
original columns such as `MMS_resid_z` remain original archived input columns;
they are not renamed into fresh target/reference residuals. Legacy override
flags (`--target-z-col`, etc.) are intentionally unsupported. Archive presence
does not change fresh residuals or neighbors. `ranking_provenance` records this.

`--archived-comparison-output` reports actual archived/fresh candidate counts,
overlap, differences, and Jaccard at pooled and section levels. It does not
assume the archived count is 507, label archived replay de novo reproduction,
or claim saved rank-residual z values equal fresh score-residual z values.

## Product Baselines

`baseline_signed_product = target_resid_z * reference_resid_z` is the exact
requested signed residual-product baseline. It is exported even when a spot
lacks coordinates, because that product is nonspatial. A negative-negative
pair can have a large positive product; it is not joint-high enrichment.
`baseline_positive_product = max(target_resid_z,0) * max(reference_resid_z,0)`
is retained as a separately named sensitivity, **not a replacement**.
Raw-score within-section z products are additional, clearly named comparators.

`ranking_comparison` exports pooled and per-section comparisons for top K and
each requested top fraction. For the signed residual product, restrictions are:

- `all_spots`: both continuous scores ranked on all common spatial-valid spots,
  with no hidden candidate or neighborhood-positivity filter.
- `positive_quadrant`: both scores ranked on the same target_z>0, reference_z>0
  subset; rank correlation is calculated on that restricted subset too.
- `cadspot_positive_support`: both scores ranked on the same positive-quadrant
  plus neighbor_reference_z>0 subset.

Both top lists use the **same universe** and deterministic ID tie-breaking.
Fractions use `ceil(fraction * n_comparison_universe)` (minimum request 1;
an empty universe produces empty lists). Top K is capped by that common
universe. Counts, actual selections, coordinate exclusions, intersection,
Jaccard, overlap fractions relative to each actual list and requested K, and
negative-negative top selections are explicit. `overlap_fraction` is overlap
divided by the actual CAD-SPOT top-list size. Comparison rank correlation uses
the unrestricted continuous priority score, not its candidate-filtered rank.
These are descriptive ranking agreements, not predictive accuracy or validation.

## Block Bootstrap

Add the following to the score-table command to request spatial uncertainty:

```powershell
  --bootstrap-output example_output/block_ci.tsv `
  --block-size 1500 --bootstrap-reps 1000 --seed 20261115 `
  --bootstrap-draws-output example_output/block_draws.tsv `
  --bootstrap-audit-output example_output/block_audit.tsv `
  --block-map-output example_output/block_map.tsv
```

Block size is an explicit scientific sensitivity choice in declared units,
not an estimated correlation range. Use several plausible widths. Blocks are
non-overlapping square grid cells `floor((xy - section_min_xy)/block_size)`;
each occupied cell belongs to exactly one section. For every replicate, draw
the original number of occupied blocks **within each fixed section**, with
replacement. Keep every selected block's full spot multiplicity. Section
block counts are fixed; section spot counts may vary. Sections themselves are
not resampled. A section with fewer than two occupied blocks fails explicitly.

Missing-coordinate spots are excluded from this bootstrap only; ordinary
score associations still include them. Bootstrap point estimates are refitted
on the coordinate-complete subset, so their n and rho may differ from the
nonspatial table result. Excluding an entire section for lack of coordinates
fails rather than silently changing section strata.

Every draw refits all residual models and section fixed effects. Requested
partial-rank transforms and spline knot locations are also fitted anew on the
resample. The default `--bootstrap-methods score_residual_spearman` reports all
available layers for the primary method. Specify any/all of the three method
names to bootstrap sensitivities too. `--confidence .95` gives percentile CIs;
undefined draw statistics remain explicit and are not used as zeros. A CI is
withheld if fewer than 80% (or fewer than two) requested draws are finite.
The draw/audit files expose failures and the exact section block/spot counts.

**These CIs are conditional on the observed sections and the chosen block
scheme. They do not establish population-level inference across sections,
remove spatial dependence by assumption, or validate candidate locations.**

## Fresh Matched Sets and Null Calibration

`cadspot_null.py generate` requires an explicit universe and numeric matching
metadata: nonnegative mean expression, detection probability in [0,1], and
positive gene length. Length construction and metadata provenance belong in
the supplied metadata audit (e.g. a versioned exon-union annotation); this
program does not fabricate lengths. Universe genes must exist in both matrix
and metadata. Missing/invalid values fail, rather than silently narrowing the
universe. Observed target/reference genes are always excluded from controls;
`--exclusions` adds a further explicit exclusion list.

Sparse input is supported without a full dense gene-by-spot export:

```powershell
& 'D:\hgsoc_cuda_env\Scripts\python.exe' -B cadspot_null.py generate `
  --expression INPUT/logcp10k_all.npz --genes INPUT/genes.tsv `
  --expression-spots INPUT/spots_all.tsv.gz --matrix-orientation spots-by-genes `
  --expression-scale log1p_CP10K --spot-key-cols sample_id barcode `
  --target-genes gene_sets/target_available19.tsv `
  --reference-genes gene_sets/reference_available20.tsv `
  --intended-target-genes gene_sets/target_intended.tsv `
  --intended-reference-genes gene_sets/reference_intended.tsv `
  --metadata INPUT/length_complete_metadata.tsv --length-col gene_length `
  --universe INPUT/explicit_universe.tsv --exclusions INPUT/excluded_genes.tsv `
  --mean-tolerance 0.25 --detection-tolerance 0.05 --length-tolerance 0.25 `
  --n-sets 500 --seed 20261115 --output example_output/matched_sets.tsv
```

All gene-list TSVs need a `gene` header (or explicit `--gene-col`). The NPZ
gene axis must exactly follow `--genes`. The matrix spot axis must exactly
follow `--expression-spots`. Orientation is required, never guessed. Dense
gene-by-spot TSVs are also accepted; optional `expression_column` in the spot
map verifies header order. The default composite keys are `(sample_id,
barcode)`, never a potentially repeated barcode alone. Calibration aligns
matrix/analysis order by these keys and requires identical spot membership.

Each target gene is a matching slot. Absolute calipers are applied on
`log1p(mean_expression)`, untransformed detection probability, and natural
`log(gene_length)`. The numbers above are explicit example calipers, not an
automatic recommendation. A randomized maximum bipartite matching assigns
one eligible unique control gene to every slot. It fails if no full unique
assignment exists and **never relaxes tolerances**. The output records matched
slot, chosen gene, all three differences, candidate count, and seed. Genes are
unique within each set; complete sets may repeat between draws. This matching
algorithm is not claimed to sample uniformly over all feasible gene sets.

Calibrate only after successful generation:

```powershell
& 'D:\hgsoc_cuda_env\Scripts\python.exe' -B cadspot_null.py calibrate `
  --expression INPUT/logcp10k_all.npz --genes INPUT/genes.tsv `
  --expression-spots INPUT/spots_all.tsv.gz --matrix-orientation spots-by-genes `
  --expression-scale log1p_CP10K --spot-key-cols sample_id barcode `
  --target-genes gene_sets/target_available19.tsv `
  --reference-genes gene_sets/reference_available20.tsv `
  --intended-target-genes gene_sets/target_intended.tsv `
  --intended-reference-genes gene_sets/reference_intended.tsv `
  --matched-sets example_output/matched_sets.tsv --spot-table INPUT/spots_all.tsv.gz `
  --technical log_nCount log_nFeature `
  --composition Epithelial_malignant_proxy ECM_Stroma CAF Macrophage T_cell `
  --scoring mean_gene_z --method score_residual_spearman `
  --output example_output/fresh_null_summary.tsv `
  --null-output example_output/fresh_null_rhos.tsv `
  --score-output example_output/fresh_observed_scores.tsv
```

Generation and calibration verify matrix, axis-map, gene-list, and generated
set provenance via the matching manifest. Early generation manifests without
`spot_key_cols` are accepted only for the original hardcoded
`[sample_id, barcode]` key; that inference is explicitly recorded. No generated
manifest is rewritten by calibration.

Scoring takes the matrix values **as supplied**, with no implicit normalization:
`mean_gene_z` averages gene-wise z-scores across all supplied keyed spots
(population SD); `mean_gene_z_by_section` fits each gene's z within section;
`mean_expression` averages supplied values directly. Constant genes fail
z-scoring rather than silently becoming invented standardized values. Dense
arrays are materialized only for the selected genes, not the entire NPZ.

`score_gene_set` is used for both the observed and every random target set.
The reference set stays fixed. The same `core.fit_pair`, covariates, section
fixed effects, method, and spline options are used for observed and random
statistics. Existing `target_score`/`reference_score` columns cannot replace
these fresh scores. New columns are `fresh_target_score` and
`fresh_reference_score`. An existing column with either fresh name is rejected
to prevent ambiguity. No failed/null-NaN control is silently dropped.

The default one-sided plus-one value is
`(1 + count(null_rho >= fresh_observed_rho)) / (1 + n_controls)`.
`--alternative two-sided` compares absolute rho magnitudes instead. This is
**gene-set empirical calibration under the stated matching scheme, not a
spatial p value**, and has Monte Carlo resolution limited by control count.
It neither replaces block uncertainty nor tests candidate-location enrichment.

For archival audit only, the core CLI accepts `--null-rho-file TABLE` with a
named `rho` column and optional `--null-output`. Its output is explicitly
`external_archived_model_equivalence_unverified_not_fresh_reproduction`.
That path does not generate controls or claim matching/model equivalence.

## Public Functions and Tests

`cadspot.py`: `read_table`, `validate_table`, `design`, `resid`, `fit_pair`,
`correlation`, `associations`, `z_by_sample`, `positive_product`,
`spatial_ranking`, `ranking_comparison`, `archived_candidate_comparison`,
`make_spatial_blocks`, `sample_blocks`, `block_bootstrap`,
`empirical_calibration`, `write_tsv`, `file_sha256`, `write_manifest`, `main`.

`cadspot_null.py`: `gene_list`, `ExpressionData.aligned`, `load_expression`,
`audit_gene_sets`, `score_gene_set`, `generate_matched_sets`,
`calibrate_gene_sets`, `main`. Function signatures/docstrings define the same
contracts used by the CLIs; matching/calibration are also callable in Python.

```powershell
& 'D:\hgsoc_cuda_env\Scripts\python.exe' -B -m unittest discover -s tests -v
```

Focused synthetic tests cover ID/numeric/missing validation; score-residual
versus partial Spearman; tied spline knots; exact six-neighbor/self exclusion;
missing coordinates and insufficient neighbors; section isolation and
deterministic ties; negative-negative products and requested baselines; fair
restricted top lists/rank correlation; archive separation; bootstrap refitting,
fixed section block strata and seed repeatability; strict matching/exclusions
and infeasible assignments; available-gene audits/no imputation; keyed sparse
alignment; observed/random model equivalence for all three estimators; finite
empirical controls; and sparse generation/calibration CLI integration.
Temporary test artifacts are confined to this directory and removed afterward.
