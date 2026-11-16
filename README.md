# CAD-SPOT 2.3.0

Composition-Aware Decomposition for Spatial Prioritization of Transcriptomic
coordinates: an auditable workflow connecting declared scores, composition/model
sensitivities, spatial uncertainty and physical-coordinate priorities.

For the software-oriented manuscript scope and reproducibility entry points, see
[SOFTWAREX.md](SOFTWAREX.md). [Licence.txt](Licence.txt) is an identical MIT-license
alias of [LICENSE.txt](LICENSE.txt), provided for the SoftwareX repository checks.

## Applications and scope

- **Primary worked example:** ovarian HGSOC, GEO GSE211956, eight Visium samples
  GSM6506110-GSM6506117. `Primary/` contains newly reconstructed count-derived
  scores, not previous image-prediction/deconvolution targets.
- **Cross-cancer transfer assessment:** PDAC GSE274103, five sections
  GSM8443449-GSM8443453, under `PDAC_Transfer/`. These data were examined before;
  this is not an untouched external holdout.
- Seven score-level benchmark scenarios and ovarian/PDAC background injections
  test conditional recovery. All favorable/unfavorable results are retained.
- The interface accepts other score pairs, but the present empirical examples
  do not validate every possible biological program.

## Install and replay

```bash
python -m pip install -r requirements.txt
python -B -m unittest discover -s Software/tests -v
python -B -m unittest discover -s Analysis -p 'test_*.py' -v
python -B Software/cadspot.py --input Primary/score_table.tsv.gz --target target_score --reference reference_score --technical log_nCount log_nFeature --composition Epithelial_proxy ECM_proxy Contractile_proxy Myeloid_proxy Tcell_proxy --coordinate-units full_resolution_pixels --output replay_associations.tsv --coordinate-output replay_coordinates.tsv.gz
```

Python 3.10+; MIT license; no nonacademic restriction. Recorded runs used Python
3.12.14. Core-only installation uses `Software/requirements.txt`. See
[Software/README.md](Software/README.md) for schemas and estimators and
[Analysis/README.md](Analysis/README.md) for both public-count reconstruction
routes, bootstrap, matching, figures and complete execution commands.

## Design decisions

Default priority is the positive-part product of within-section standardized
target/reference residuals. Neighborhood support is an annotation, with explicit
optional weighted modes. Candidate center/scale/fallback are exported; the
adjustable outlier cutoff does not control coordinate FDR.

RNA proxies are not measured cell fractions. Adjustment conditions on genuine
abundance-coupled biology as well as mixture effects. Keep technical-only and
composition layers together. Matched-gene probabilities concern the chosen panel,
while spatial-block intervals are conditional on observed sections/fixed scores.
H&E registration localizes hypotheses, not histopathological validation.

## What changed in 2.3.0

The new decision-audit CLI exports model-specific budget memberships, ranks,
coordinates and a printable report. `Decision_Utility/README.md` gives one-command
use and complete budget, input-perturbation and second-program-pair evaluations.
Agreement is descriptive, not biological validation or automatic model choice.
All losses against alternatives and noncomputable perturbations are retained.

## Earlier additions and mathematical scope

Ranking arithmetic is unchanged from buffered v2.1.1. New additions are the
ovarian Matrix Market adapter, coverage-only panel harmonization and section-pitch
bootstrap geometry. Full-resolution image coordinates stay in the exports.
`Study_Design.md` records choices before observing the new ovarian results.

`Utility_Benchmark/` supplies a standard OLS reference, frozen v2.1.0 comparison,
paired recovery tests and single-machine implementation measurements. Equivalent
OLS reproduces the default, and alternative baselines can perform better.
The contribution is the reproducible analysis chain, not novel OLS mathematics
or universal accuracy/speed superiority over spatial platforms.

Optional PDAC-only contextual fits and their scope/convergence limitations remain
in the supplement. They are not ovarian model fits. Large matrices/images and
third-party environments are retrieved through documented public links.

## Versions and historical material

The current software release is tagged `v2.3.0`. Earlier tags, legacy root
examples and old source files are retained for history. Active ovarian sections
are identified by `Primary/input_manifest.json`; PDAC is separately identified by
`PDAC_Transfer/input_manifest.json`. Follow these manifests and current commands,
not stale files from previous PDAC-primary or ambiguous-sample examples.
