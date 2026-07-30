# CAD-SPOT

CAD-SPOT (Composition-Aware Decomposition for Spatial Prioritization of Transcriptomic coordinates) is a lightweight command-line workflow for auditing spot-level spatial transcriptomic score associations and exporting ranked candidate coordinates. The workflow is score-agnostic: the example uses a stress-correlate target and hypoxia/glycolysis reference, but users can provide other target-reference score pairs if appropriate covariates are available.

Repository: https://github.com/ASDQWRQ/CAD-SPOT. A Zenodo DOI should be added here after archival deposition if this repository is used for manuscript submission.

## Installation

CAD-SPOT requires Python 3.10 or later. Install dependencies with:

```bash
pip install -r requirements.txt
```

The demonstration script uses only `pandas`, `numpy`, `scipy` and `scikit-learn`.

## Input

CAD-SPOT reads a tab-delimited spot-level table. The example file is `example_spot_table.tsv`.

Required columns in the example:

- `sample_id`: section or sample identifier.
- `barcode`: Visium spot barcode or another stable coordinate identifier.
- `x`, `y`: pixel-space or array-space spot coordinates from the upstream spatial object.
- `Core_MRS`: example target score.
- `HG_score`: example reference score.
- `log_nCount`, `log_nFeature`: technical covariates.
- `Epithelial_malignant_proxy`, `ECM_Stroma`, `CAF`, `Macrophage`, `T_cell`: example composition covariates.

Users may replace `Core_MRS`, `HG_score` and the covariate columns with their own score names by changing the command-line arguments.

Possible score-pair examples include ligand-receptor co-activity, immune-exclusion signatures, transcription-factor target modules or any other user-defined target-reference score pair.

## Quick start

```bash
python cadspot.py --input example_spot_table.tsv --target Core_MRS --reference HG_score --technical log_nCount log_nFeature --composition Epithelial_malignant_proxy ECM_Stroma CAF Macrophage T_cell --output example_cadspot_associations.tsv --coordinate-output example_cadspot_coordinates.tsv
```

The default candidate rule uses `--robust-z-cutoff 3.0`. This parameter can be changed by the user, but the manuscript analyses used the default value.

## Outputs

`example_cadspot_associations.tsv` reports raw, technical-adjusted and composition-adjusted Spearman correlations.

`example_cadspot_coordinates.tsv` retains the original input columns and adds:

- `target_resid`, `reference_resid`: residualized target and reference scores.
- `target_resid_z`, `reference_resid_z`: within-sample residual z-scores.
- `neighbor_reference_resid_z`: local neighborhood support for the reference residual.
- `local_support`: local joint support based on neighboring target/reference residuals.
- `priority_score`: average of residual and neighborhood support terms.
- `priority_robust_z`: within-sample robust z-score for the priority score.
- `cadspot_candidate`: binary candidate-coordinate flag.

## Interpretation

CAD-SPOT is a coordinate-prioritization aid, not a causal adjustment model and not a histologic validation assay. A candidate coordinate means that the spot has high residual target-reference co-enrichment relative to its section-specific score distribution and local neighborhood support. Candidate coordinates should be treated as follow-up targets for pathology review, multiplex imaging, spatial proteomics or targeted in situ assays.

In the simulation terminology used in the manuscript, the architecture field represents broad composition-linked mixture and density gradients, whereas the local microfield represents spatially restricted coupling that remains after composition-aware decomposition.

## Reproducibility

The script is table-based and does not require a Seurat, Scanpy, Giotto or Squidpy object. This is intentional: users can export scores and covariates from any upstream workflow, run the same command-line audit, and archive the flat association and coordinate tables with the manuscript.

Processed manuscript output tables are provided in `data/processed/`. These tables include density-adjustment summaries, matched-null calibration outputs, residual-method sensitivity analyses, simulation summaries, external-audit results and the compressed candidate-coordinate table.
