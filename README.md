# CAD-SPOT 2.1.0

Composition-Aware Decomposition for Spatial Prioritization of Transcriptomic
coordinates: a table-first workflow connecting explicit score provenance,
association sensitivities, spatial uncertainty and inspectable coordinate ranks.

## Install and run

```bash
python -m pip install -r requirements.txt
python -B -m unittest discover -s Software/tests -v
python -B -m unittest discover -s Analysis -p 'test_*.py' -v
python -B Software/cadspot.py --input Primary/score_table.tsv.gz --target target_score --reference reference_score --technical log_nCount log_nFeature --composition Epithelial_proxy ECM_proxy Contractile_proxy Myeloid_proxy Tcell_proxy --coordinate-units full_resolution_pixels --output replay_associations.tsv --coordinate-output replay_coordinates.tsv.gz
```

Core-only installation uses `Software/requirements.txt`. Python 3.10+; MIT license;
no restriction on nonacademic use. The recorded analyses used Python 3.12.14.
See [Software/README.md](Software/README.md) for the input schema and estimators,
and [Analysis/README.md](Analysis/README.md) for all public-count reconstruction,
spatial bootstrap, matched-set, ranking simulation and figure commands.

## Active reproducible example

The primary example is reconstructed from GSE274103 public filtered UMI counts
and registered spatial metadata, without archived scores/candidate flags. It has
23,392 QC-passing spots and 650 positive-product candidates. The raw/technical/
composition-adjusted pooled correlations are 0.465/0.351/0.339. Available-gene
coverage, score parameters, source checksums and complete derived outputs are
recorded. Counts and correlations are results, not fixed selection targets.

The separate-study GSE203612/GSM6177618 run retains 1,759 spots and 48 candidates.
**Disease identity is unresolved:** PDAC sample naming conflicts with ovarian
metadata. It is a qualified technical transfer example, not confirmed external
PDAC biological validation. See [External/SOURCE_ISSUE_NYU_PDAC1.md](External/SOURCE_ISSUE_NYU_PDAC1.md).

## Design decisions

- Default rank is `max(target_residual_z,0)*max(reference_residual_z,0)`.
  Six-neighbor reference support is an annotation; neighborhood-weighted modes
  are explicit optional sensitivities. Simulations show geometry-specific tradeoffs.
- Candidate threshold is adjustable and its actual center/scale are exported.
  The positive-part scores use SD fallback here. This is not Gaussian significance
  or coordinate-level FDR control; null simulations still nominate coordinates.
- Spatial block intervals refit score residual models but hold upstream scores/
  proxies fixed; they are conditional on observed sections, not patient inference.
- Strict matched-gene calibration uses the SAME observed/control scoring and
  residual model. Its P value concerns gene-set choice, not spatial independence.
- Proxy adjustment describes dependence on transcript proxies, not causal
  correction, cell deconvolution or biologically truer residuals.

Small optional model exports and scripts contextualize new candidates against
preserved STAMP/POPARI/SpaHDmap fits. They are not a superiority benchmark and
were not refitted on the expanded primary universe; see
[Analysis/comparators/README.md](Analysis/comparators/README.md).

## Version and historical files

The manuscript analysis release is tagged `v2.1.0`. Older tagged versions and
legacy root examples/Source_Data files are retained for history, not the active
public-matrix primary example. Follow the explicit `Software/` and `Analysis/`
commands above rather than mixing old frozen-score outputs into this release.
Large public matrices/images, optional checkpoints and external runtime/source
trees are not bundled; download/reconstruction commands are supplied.
