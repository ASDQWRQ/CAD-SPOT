# Reproduce the integrated CAD-SPOT 2.1.0 example

Run commands from the repository root. Python 3.10+ is required; the recorded
analyses used Python 3.12.14. Install `requirements.txt`, or only
`Software/requirements.txt` for the table-based core. No archived residuals,
candidate labels or old composition columns enter the current score generation.
`Primary/` and `External_EndToEnd/` contain derived tables, not new experiments.

## 1. Tests and fast table execution

```bash
python -B -m unittest discover -s Software/tests -v
python -B -m unittest discover -s Analysis -p 'test_*.py' -v
python -B Software/cadspot.py --input Primary/score_table.tsv.gz --target target_score --reference reference_score --technical log_nCount log_nFeature --composition Epithelial_proxy ECM_proxy Contractile_proxy Myeloid_proxy Tcell_proxy --coordinate-units full_resolution_pixels --output replay_associations.tsv --coordinate-output replay_coordinates.tsv.gz
```

The default ranking is the positive residual product. An optional neighborhood
term is NOT required by that default. See `Software/README.md` for input column
definitions, missing coordinates, six-neighbor rules and exported scale metadata.
Outputs beginning `replay_` avoid overwriting the recorded analysis.

## 2. Reconstruct scores from public count matrices

The primary archive is GEO GSE274103. `download-primary` fetches its official
archive and safely extracts it; it is a large download, not part of the small
software archive. `--source-root .` means the current repository root.

```bash
python -B Analysis/end_to_end.py --stage download-primary --source-root .
python -B Analysis/end_to_end.py --stage prepare-primary --source-root .
python -B Analysis/build_gene_lengths.py
python -B Analysis/end_to_end.py --stage run --folder Primary --reps 1000 --widths 1500 3000 4500 --n-sets 500
```

`build_gene_lengths.py` uses public Ensembl GRCh38 release 113. Its derived table
is already included under `Provenance/`. Rebuilding needs the full public GTF.
The preparation stage aggregates duplicate symbols, applies declared QC,
normalizes counts and computes available-gene scores/proxies. It exports
gene-specific means/SDs/weights, coverage and source checksums. Missing panel
genes are omitted and reported, not silently set to zero. Downstream full
analysis needs its generated `logcp10k_all.npz`, omitted from small archives.
This preparation starts from deposited filtered UMI matrices, not FASTQ reads.

## 3. Separate-study example

```bash
python -B External/prepare_external_inputs.py
python -B External/finalize_source_audit.py
python -B Analysis/end_to_end.py --stage prepare-external --prepared-external External/prepared/GSM6177618_NYU_PDAC1_Vis --folder External_EndToEnd/NYU_PDAC1
python -B Analysis/end_to_end.py --stage run --folder External_EndToEnd/NYU_PDAC1 --reps 1000 --widths 2394 4788 7182 --n-sets 500
python -B Analysis/common_gene_sensitivity.py
```

See `External/SOURCE_ISSUE_NYU_PDAC1.md`: sample GSM6177618 has conflicting
PDAC/ovarian source annotations. The publisher supplement labels PDAC NYU1 but
does not supply pathology. This is a qualified separate-study technical transfer,
NOT a confirmed independent PDAC biological validation or an untouched holdout.
The source audit remains part of the release even when upstream metadata changes.

## 4. Ranking ablation and figures

```bash
python -B Analysis/ranking_ablation.py
python -B Analysis/make_figures.py
python -B Comparators/recompute_overlap.py
python -B Analysis/make_comparator_figure.py
```

Simulations use 50 paired seeds in each of three score-level geometries; no
positive labels are injected in the null geometry. AP is undefined there.
These are controlled ranking tests, not fitted gene-count or histology models.
Figure generation requires registered images from stages 2 and 3; recorded
source manifests contain acquisition-time absolute paths. If relocated, rerun
preparation so the image paths point to the new cache.

The comparator step reuses candidate-blind fits on a fixed 23,013-spot scope,
then recomputes overlap using NEW candidate labels. It does not claim those
models were refitted on the expanded 23,392-spot universe. Source scripts,
configurations, compatibility patches and stopping records are in
`Analysis/comparators/`; small factor exports are bundled for replay. Training
them from scratch requires the original model input and separate third-party
environments. POPARI is a budget-limited checkpoint, not converged. The overlap
null is approximate block-label exchangeability, not comparative accuracy.

## Result interpretation and provenance

- Raw, technical and composition layers are retained separately. Residualization
  is compositional sensitivity analysis, not causal correction.
- Bootstrap draws refit residual models/design scaling while holding original
  gene scores/proxies fixed; inference is conditional on the observed sections.
- Matched sets use the same expression scoring and residual model as observed.
  Calipers never relax. The gene-set P value is not coordinate-level FDR.
- The outlier cutoff is user-adjustable. Its scale method is exported; the
  positive product uses the SD fallback in these examples, not a Gaussian test.
- Run manifests are completion records; partial output folders do not establish
  successful completion. Full commands overwrite the corresponding Results
  files, so preserve previous runs before executing them.

The main examples are 23,392/650 primary spots/candidates and 1,759/48 external
spots/candidates. These counts are outputs, not fixed selection targets. The
repository preserves earlier tagged analyses separately; do not substitute their
frozen-score counts or statistics into this release.
