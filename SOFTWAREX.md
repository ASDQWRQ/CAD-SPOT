# SoftwareX manuscript companion

CAD-SPOT 2.3.0 is described in the manuscript "CAD-SPOT: Model-aware coordinate
prioritization for spatial transcriptomics". This companion note records the
software-oriented scope without changing the released estimator or candidate rule.

## Intended task

Given declared target/reference scores, suitable covariates and within-section
coordinates, the software exports source-keyed priorities and exposes how a
fixed inspection list depends on technical, composition and nonlinear models.
It does not estimate cell fractions, learn spatial topics, or validate tissue
activation states. The independent OLS-product implementation can reproduce the
default priorities; the contribution is the integrated audit and export workflow.

## Entry points

- `Software/README.md`: installation, table schema, association and ranking CLI.
- `Software/decision_audit.py`: model-specific budget memberships and report.
- `Decision_Utility/README.md`: fixed-budget tests, perturbations and second pair.
- `Analysis/README.md`: public-count reconstruction and analysis drivers.

Run from the repository root:

```bash
python -m pip install -r requirements.txt
python -B -m unittest discover -s Software/tests -v
python -B -m unittest discover -s Analysis -p 'test_*.py' -v
```

The root requirements include count-preparation dependencies such as `h5py`.
For table-only use, `Software/requirements.txt` provides the core dependencies.
The prepared examples allow table-based replay. Raw count reconstruction and
count perturbations require separately retrieved public matrices, whose locations
and checksums are recorded in the input manifests.

## Evidence scope

The complete code-and-results package is the named
[Reproducibility_Archive.zip asset](https://github.com/ASDQWRQ/CAD-SPOT/releases/download/v2.3.0/Reproducibility_Archive.zip)
on the [fixed v2.3.0 release](https://github.com/ASDQWRQ/CAD-SPOT/releases/tag/v2.3.0).
It is hosted externally rather than uploaded as a journal ZIP attachment.
The journal supplements provide full technical text and editable tables; a
separate reproducibility access note gives the release link and checksum.
This hosting/documentation update does not alter the tag, algorithms or results.

The primary example is ovarian GSE211956; PDAC GSE274103 is a retrospective
cross-cancer application, not an untouched validation cohort. Injection recovery
is computational evidence, not pathology ground truth or measured experimental
cost saving. Conditional-model agreement is not an automatic consensus winner.
All unfavorable comparisons and noncomputable perturbation repetitions remain
part of the recorded evidence. Existing v2.2.0 kernel timings exclude the new
v2.3.0 audit report and are not a full-pipeline speed claim.

`Licence.txt` and `LICENSE.txt` contain the same MIT license. The documentation-only
alias does not alter release v2.3.0 or its frozen source.
