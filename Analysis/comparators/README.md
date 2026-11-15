# Limited overlap/enrichment comparisons

This directory executes official methods, not substitutes. No classification accuracy or winner ranking is justified. POPARI and STAMP run on five PDAC sections; SpaHDmap runs the first section with verified registered H&E in fixed sample order. The primary freeze is `Source_Data/fresh_coordinates.tsv.gz` (755 candidates, not archived 507).

## Reproduction

Run from this directory using a scientific Python3.12 environment with the core numeric packages and CUDA-capable PyTorch already installed. Exact executed package versions are recorded in `requirements-core.txt`, `requirements-staged.txt`, and the run records. Additional packages are isolated through `COMPARATOR_DEPS`; the shared core environment is not modified. Source repositories are in `vendor`; restore the recorded official commits and listed compatibility patches before running.

```powershell
$env:COMPARATOR_DEPS = (Join-Path $PWD 'deps')
python setup_sources.py
python -m pip install --target deps --no-deps -r requirements-staged.txt
python prepare_inputs.py --raw ../../Data/GSE274103/raw --spatial-root ../../Data/GSE274103/prepared --candidates ../../Source_Data/fresh_coordinates.tsv.gz --out inputs
python test_comparators.py
python run_methods.py stamp --batch-size 512 --out results/stamp/reproduction_seed20261115
python run_methods.py popari --popari-dtype float64 --popari-iterations 50 --out results/popari/reproduction_seed20261115
python run_methods.py spahdmap --batch-size 4 --out results/spahdmap/reproduction_seed20261115
python enrichment.py --outputs results/stamp/reproduction_seed20261115 results/popari/reproduction_seed20261115 results/spahdmap/reproduction_seed20261115 --out summary_reproduction
```

The environment above is the executed package freeze, not a claim that every optional distribution constraint is satisfied; `dependency_resolution.json` records remaining declared version conflicts. `common.py` adds the isolated dependency directory to `sys.path`. Existing nonempty method directories are never overwritten; use a fresh results directory for a new run. The `retry_numpy2` suffix identifies the actual completed STAMP artifact, not a different algorithm. Execution on the reference8GB GPU used float64 POPARI and SpaHDmap batch4; peak combined VRAM was approximately6.4GB. The source image-feature inference helper separately uses batch32 and allocates a large host-memory buffer. Sequential execution is suitable for reproduction.

The commands deliberately use new `reproduction_seed20261115` directories so the
released computed results remain unchanged. Choose another unused directory for
a second reproduction. `finalize_comparison.py` audits the shipped final run
names and summaries; it is not required to summarize a newly named run.

## Fixed Protocol

A fixed50-spatial-update POPARI budget was used after5 nonspatial warmup updates, with the official SVD initialization, shared metagenes, section-specific affinities, and float64. The exact50-update checkpoint was frozen from a longer parent run; `export_popari_checkpoint.py` verifies the checkpoint interval before exporting it. The direct50-update command above reproduces that budget. This is a budget-limited exploratory comparison, NOT a converged fit: final relative objective change was4.09e-7, but relative loading change0.007573 exceeded the1e-3 criterion. The longer-run candidate-blind stability rule required both objective change<1e-4 and loading change<1e-3 for10 consecutive iterations after at least50 spatial updates; it was not met at the reported checkpoint.

STAMP used800 maximum epochs and the official early-stopping rule after100 epochs, batch512. Completed epoch counts were140/127/177/150/154. SpaHDmap used its full official sequence: NMF initialization, image pretraining (nominal maximum5000, official early stop at2600), all5000 GCN steps, Voronoi extension, and nominal2000 final fusion steps. Only final fusion-derived `SpaHDmap_spot` scores are comparator evidence. K=10, seed=20261115, and3000 pooled top log-normalized-variance genes were fixed independently of candidate enrichment. This gene rule is not claimed to be Seurat variance-stabilized HVG selection. `--smoke` is explicitly marked non-evidence and rejected by enrichment.

## Final Run Directories

- `results/stamp/k10_seed20261115_retry_numpy2`: five independent official STAMP fits.
- `results/popari/k10_seed20261115_budget50`: the actual five-section joint POPARI50-update model and exports; not the initialization checkpoint.
- `results/spahdmap/k10_seed20261115_production`: one representative registered H&E section and full official fusion pipeline.
- `summary_popari_stamp`: five-section two-method comparison.
- `summary_all`: final three-method representative-section extension plus the five-section two-method summaries, with explicit scope and fit-budget columns. No manuscript or figure files are modified.

## Outputs and interpretation

- `inputs/universe_audit.tsv.gz`: every evaluated spot, inclusion flags, raw count QC; `eligible_spots.tsv.gz` and per-section h5ad freeze the common universe.
- `inputs/manifest.json`: exact raw H5 paths/hashes, registered image/scalefactor paths, representative section and gene-selection rule.
- Each method run: `run.json` includes official git SHA, all local source patches, device/packages/settings/status; `console.log` retains failures and training losses. POPARI and STAMP additionally export machine-readable objectives. A run reaching a fixed iteration budget is not automatically proof of convergence.
- `spot_factors.tsv.gz`: raw barcode plus sample ID, posterior proportions only for STAMP; descriptive row-normalized final embeddings for POPARI and SpaHDmap. `raw_factor_scores.tsv.gz` retains native scales. Dominant factors are descriptive argmax domains; SpaHDmap also exports downstream fixed-K KMeans domains separately. Nonnegative embeddings are not calibrated probabilities.
- `candidate_factor_enrichment.tsv`: effects, observed/expected overlap and approximate block p/q per factor/section. All candidate labels come from the frozen new table. `within_section_summary.tsv` exposes common-universe attrition; `pairwise_domain_overlap.tsv` reports all domain pairs without asserting factor correspondence.

When three methods are supplied, only representative-section results have all three in their common scope. Five-section inference must not silently be described as a three-method comparison. Independent STAMP topic numbers are arbitrary across sections. Primary spatial blocks span 8 median spot pitches; widths 4/12 are sensitivity analyses. Equal-occupancy exchange is an approximate null and can be underpowered on irregular tissue. P values are deliberately undefined if fewer than 8 blocks or fewer than half the spots can move. No naive spot Fisher test is substituted.

Fold enrichment uses the overall candidate rate as a descriptive denominator. The conditional spatial-block null may have a different expectation because unmatched edge blocks remain fixed; `block_null_expected_overlap` exposes that distinction. All methods use the same expression data as CADSPOT, so enrichment is a descriptive relationship, not independent biological validation. In the final representative-section joint universe, width12 permits only37.12% of spots to move, below the50% threshold; its P/Q values are therefore undefined for all three methods. Primary width8 and sensitivity width4 remain testable.

SpaHDmap produced343 all-zero native vectors out of4804 spots, including36 of148 candidates. All raw vectors are retained; zero vectors are explicitly unassigned and are never given pseudocounts or artificial domains. Normalized SpaHDmap outputs contain4461 valid spots/112 candidates. The representative-section three-method comparison uses this joint valid intersection for every method; the other four sections retain two-method scope. Candidate loss is24.32%, versus7.14% overall spot loss, so this conditional comparison has nonuniform coverage and cannot establish method superiority. `all_spot_assignments.tsv.gz`, `raw_factor_scores.tsv.gz`, and `coverage.json` expose the attrition. Model training was completed; `recover_spahdmap_export.py` only repaired export of the existing final native scores, without retraining.

For source-only distribution, `setup_sources.py` restores the pinned official repositories and recorded compatibility patches; source redistribution is not required. Retain `source_versions.json` and the completed run records containing exact commits and source diffs. Restored packages must include runtime assets, particularly SpaHDmap's `SpaHDmap/utils/color.csv`; Python files alone are insufficient. `RunLog` falls back to the recorded official commit and exact local patch when `.git` or the Git executable is unavailable. Package-version freezes do not bundle the shared environment. Seeds and versions are recorded; bitwise equivalence across hardware is not asserted.

## Minimal compatibility patches

- POPARI: defer `pymetis` import to the unused partition-downsampling call. The chosen single-level fit never partitions or bins; no substitute partition algorithm is introduced.
- STAMP: `np.Inf` becomes `np.inf` for numpy 2; append the already-computed epoch objective to a list. No model likelihood, optimizer, or stopping condition is changed.
- SpaHDmap: use equivalent `plt.get_cmap` for the removed matplotlib `cm.get_cmap` visualization API; model math is unchanged.

Any further patches and runtime constraints are recorded in the feasibility note and each `run.json`.
