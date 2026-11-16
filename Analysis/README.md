# Reproduce CAD-SPOT 2.2.0: ovarian worked example and PDAC transfer

Run commands from the extracted release root with Python 3.10+ and
`requirements.txt`. Analyses were recorded with Python 3.12.14. Use
`Software/requirements.txt` for table-core-only installation. `Study_Design.md`
records fixed settings and the retrospective roles of the two cohorts.

## 1. Tests and fast execution

```bash
python -B -m unittest discover -s Software/tests -v
python -B -m unittest discover -s Analysis -p 'test_*.py' -v
python -B Software/cadspot.py --input Primary/score_table.tsv.gz --target target_score --reference reference_score --technical log_nCount log_nFeature --composition Epithelial_proxy ECM_proxy Contractile_proxy Myeloid_proxy Tcell_proxy --coordinate-units full_resolution_pixels --output replay_ovarian_associations.tsv --coordinate-output replay_ovarian_coordinates.tsv.gz
python -B Software/cadspot.py --input PDAC_Transfer/score_table.tsv.gz --target target_score --reference reference_score --technical log_nCount log_nFeature --composition Epithelial_proxy ECM_proxy Contractile_proxy Myeloid_proxy Tcell_proxy --coordinate-units full_resolution_pixels --output replay_pdac_associations.tsv --coordinate-output replay_pdac_coordinates.tsv.gz
```

The core is score-agnostic; the 23/21 intended panels and five RNA proxies belong
to this example. Schemas, estimators, missing positions, tie rules and all export
columns are documented in `Software/README.md`. Replay names avoid overwriting
recorded Results. Default positive-product ranking does not require positive
neighborhood support. Its candidate cutoff is descriptive, not coordinate FDR.

## 2. Rebuild ovarian scores from public matrices

```bash
python -B Analysis/download_ovarian.py
python -B Analysis/prepare_ovarian.py
python -B Analysis/build_gene_lengths.py
python -B Analysis/end_to_end.py --stage run --folder Primary --reps 1000 --bootstrap-in-pitches --n-sets 500 --seed 20261116
```

GEO GSE211956 supplies all eight Visium samples GSM6506110-GSM6506117. The
downloader supports byte-range resumption, validates response positions/length,
and safely extracts the complete TAR. Spatial ZIP CRC and count gzip readers
check files during preparation. Source checksums are retained. Missing large
raw inputs are fetched; no image-prediction targets, historical deconvolution
scores or selected-barcode list is used. Source single-cell data are not newly
analyzed by these commands.

The adapter reads Matrix Market gene/barcode lists, sums duplicate symbols,
aligns positions by barcode, applies declared QC, and computes log1p CP10K scores.
Before fitting, each panel uses its common nonconstant genes across sections;
all exclusions and per-gene means/SDs/weights are exported. `logcp10k_all.npz`
and per-section count matrices are needed for full matching/reconstruction and
are rebuilt from public files, not included in small replication archives.

Ensembl-113 exon-union lengths are already supplied in `Provenance/`; rebuilding
them needs the public GTF. They do not represent assay probe-region length.

## 3. Rebuild the PDAC transfer application

```bash
python -B Analysis/end_to_end.py --stage download-pdac --source-root .
python -B Analysis/end_to_end.py --stage prepare-pdac --source-root . --folder PDAC_Transfer
python -B Analysis/end_to_end.py --stage run --folder PDAC_Transfer --reps 1000 --bootstrap-in-pitches --n-sets 500 --seed 20261116
python -B Analysis/cross_cancer_coverage.py
```

GSE274103 supplies the five PDAC sections GSM8443449-GSM8443453. The folder flag
is important: it keeps PDAC outputs separate from ovarian `Primary/`. These data
were explored in earlier versions; the run is retrospective cross-cancer
assessment, not untouched validation. The ambiguous GSE203612 sample from an
earlier version is not the current main transfer case.

Both cohorts retain full-resolution image x/y. For bootstrap only, coordinates
are divided by each section's median nearest-other-spot pitch. Widths 6/12/18
are capture-grid units, not micrometers; the 12-pitch analysis is primary. Exported
transforms make the conversion explicit. Conditional intervals hold original
scores/proxies fixed and refit residual/design models each draw.

## 4. Task-level evaluation and optional context

```bash
python -B Utility_Benchmark/benchmark.py prepare --data-root .
python -B Utility_Benchmark/benchmark.py accuracy
python -B Utility_Benchmark/benchmark.py equivalence
python -B Utility_Benchmark/benchmark.py performance
python -B Utility_Benchmark/make_report.py
python -B Analysis/make_figures.py
python -B Analysis/make_comparator_figure.py
```

Run performance after other numerical jobs finish; it uses five serial repeats
per build/size and one BLAS/OMP thread. All simulation and injection outcomes
are retained, including cases favoring alternatives. The OLS reference is
mathematically equivalent for default ranking and has fewer outputs. Timing
does not include download, preprocessing, bootstrap or matching.

The existing 50-seed neighborhood experiments can be replayed with
`python -B Analysis/ranking_ablation.py`. The optional PDAC-only model-overlap
replay is `python -B Comparators/recompute_overlap.py`. It reuses candidate-blind
STAMP/POPARI/SpaHDmap fits on 23,013 PDAC spots, not ovarian tissue. Source,
compatibility patches and stopping records are under `Analysis/comparators/`;
retraining requires public inputs and separate third-party environments.
POPARI is budget-limited/nonconverged. Factor overlap is context, not a detection
accuracy benchmark. No package imports archived model labels as ovarian inputs.

## Presentation and relocation

Figures 1/3/4 show the ovarian worked example, Figure 2 task evaluation and
Figure 5 PDAC transfer. Figure S1 is ablation; Figure S2 is PDAC model context.
Figure creation needs registered images. Source manifests preserve actual
acquisition-time paths; when relocating the project, rerun preparation to point
to the new public-data cache. Tables and source records identify cohort roles
explicitly; do not mix earlier PDAC-primary releases into current outputs.
