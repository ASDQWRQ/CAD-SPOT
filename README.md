# CAD-SPOT

Composition-aware association auditing and spatial coordinate prioritization.

## Revision source: 2.0.0-review

CAD-SPOT starts from a table of target/reference scores, technical covariates,
optional composition proxies and spatial coordinates. It reports raw,
technical-adjusted and composition-adjusted associations and exports the
residuals, six-neighbor support and ranking components for every coordinate.
It does not estimate cell fractions or make causal or clinical predictions.

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s Software/tests -v
python Analysis/run_reanalysis.py --stages core
```

The supplied frozen score table reproduces 23,013 PDAC spots, pooled adjusted
rho 0.199254 and **755 freshly computed candidates**. The earlier 507-candidate
list reused archived ranking fields; these no longer override fresh residuals.

## Included analyses

- OLS score-residual Spearman, genuine partial Spearman and additive-spline
  sensitivity, with explicit estimator and adjustment-layer labels.
- Section-stratified spatial block bootstrap, with all models refitted per draw.
- Continuous-priority versus signed residual-product comparisons on common
  spot universes, including the positive residual quadrant.
- Strict expression/detection/exon-union-length matched gene-set generation
  and internally paired observed/null calibration.
- Limited contextual comparisons with the official POPARI, STAMP and SpaHDmap
  implementations. These are overlap analyses, not accuracy benchmarks.

**Two distinct score inputs:** the frozen-score demonstration (rho 0.199) and
the independently reconstructed 19/20-gene expression audit (rho 0.324) are
not interchangeable. The gene-set empirical P value does not calibrate the
frozen-score estimate or individual candidate coordinates.

## Documentation

- [Full reproduction commands](Analysis/README.md)
- [Core CLI, formulas and input contract](Software/README.md)
- [Official-model comparisons and limitations](Analysis/comparators/README.md)

The core tests use Python 3.12. Optional raw-data preparation requires h5py;
plotting and Excel export have separate requirements. Model comparisons require
their own scientific/CUDA environment and pinned official sources, restored by
`python Analysis/comparators/setup_sources.py`. They are not core dependencies.

Public expression counts and registered images are downloaded from NCBI GEO
GSE274103 and GSE203612 using `Analysis/download_data.py`; they are not bundled
in this repository. Source outputs, coverage lists and run manifests accompany
the workflows. Original pre-revision example files are retained for historical
reference and are not the inputs for the revised manuscript analyses.

## Interpretation

The default robust priority-Z cutoff of 3 is adjustable and heuristic.
Block-bootstrap intervals are conditional on the observed sections. Matching
does not control FFPE probe-target coverage or GC effects. POPARI's reported
50-update result is budget limited, not claimed converged. SpaHDmap zero-vector
spots remain unassigned, and the three-method comparison uses their explicitly
reported valid-output intersection.

CAD-SPOT source is MIT licensed. Third-party methods retain their own licenses.
