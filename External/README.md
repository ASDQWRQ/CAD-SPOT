# Separate-study input reconstruction

## Ready input

`prepared/GSM6177618_NYU_PDAC1_Vis/` contains:

- `counts.npz`: SciPy CSR, spots by unique gene symbols, nonnegative integer
  deposited counts; duplicate symbols are summed. No normalization or scoring.
- `genes.tsv`: one `gene` column, in exact count-matrix column order.
- `positions.tsv`: exact count-matrix row order, `barcode`, `x`, `y`,
  `in_tissue`, original array/image coordinates, technical count summaries,
  `passes_locked_qc`, and scaled image coordinates. There is **no `sample_id`
  column**, so the shared `prepare_section` may insert it.
- `positions_all.tsv`: all 4,992 capture-array locations with expression-presence
  flags; do not treat this as the count-matrix row axis.
- `barcodes.tsv`, `source_features.tsv`, `feature_to_gene.tsv`: axis/provenance
  audits. Only `barcodes.tsv` separately records sample identity.

The matrix has 1,789 rows and 33,514 unique gene symbols. All expression rows
join unique registered coordinates. The locked QC rule flags 1,759 rows as
passing. No additional filtering is applied to the supplied matrix. Apply the
same row mask to counts and positions in the shared main preparation function.

Coordinates are full-resolution image pixels: `x=pxl_col_in_fullres` and
`y=pxl_row_in_fullres`, image origin top left. Array row/column are not
micrometers. Multiply by the appropriate scale factor for the corresponding
hires or lowres image. All 1,789 scaled centers are inside both supplied images.
The median nearest-other-spot distance is 399 full-resolution pixels.

```python
from pathlib import Path
import pandas as pd
from scipy import sparse

folder = Path("External/prepared/GSM6177618_NYU_PDAC1_Vis")
counts = sparse.load_npz(folder / "counts.npz").tocsr()
genes = pd.Index(pd.read_csv(folder / "genes.tsv", sep="\t")["gene"])
positions = pd.read_csv(folder / "positions.tsv", sep="\t")
assert counts.shape == (len(positions), len(genes))
assert "sample_id" not in positions
```

Use the existing `Analysis/end_to_end.py::prepare_section` for all QC and scoring.
This folder defines no new scoring or composition-proxy panels. Later module
proxies remain transcript-derived covariates, not validated cell proportions.

## Source status

Read `SOURCE_ISSUE_NYU_PDAC1.md` before writing biological validation claims.
Publisher TableS1 lists `PDAC NYU1` but pathology is `NA`; inconsistent GEO and
BioSample ovarian annotations remain unresolved. This is suitable as a
qualified independent-study full-workflow technical example, not a
disease-confirmed external PDAC cohort.

**Identity decision: unresolved.** TableS1 cell A14 is `PDAC NYU1` and B14 is
`NA`. This supports the published grouping label, not a verified resolution of
the GEO/BioSample disease conflict. The study is separate from GSE274103, but
this expression matrix has already appeared in the earlier external audit:
it is **not an untouched holdout**, nor a newly verified independent patient.

- `source_manifest.json`: official URLs, SHA-256, matrix axes, coordinate units,
  QC flags and identity scope.
- `source_checksums.tsv`: source files and publisher supplement checksums.
- `Input_Verification.json`: independent saved-file integrity and alignment tests.
- `metadata/NYU_PDAC1_identity_audit.json`: published table and exact row.
- `metadata/Barkley2022_publisher_Supplementary_Tables.xlsx`: valid original
  publisher workbook; its extracted TableS1 is also provided as TSV.
- `spatial/GSM6177618_NYU_PDAC1_Vis/`: registered images, positions and scale factors.

## Scope and reproducibility

Only public source download, lossless input formatting and input validation were
performed here. No scores, composition proxies, association models, candidate
ranking, bootstrap or random controls were computed by these source-preparation
scripts. The sample was selected by availability and study identity, not by a
score result. The previous manuscript, software and submission archives were
not modified.

The source preparation can be rerun with `prepare_external_inputs.py
--skip-documentation-fetch`; do not rerun while downstream analyses are reading
the inputs. `finalize_source_audit.py` validates saved inputs without changing
their matrix/axis files. It requires numpy, pandas, scipy, h5py and openpyxl.

The HD backup download is stopped and not part of this analysis. Incomplete
files remain in `raw/Visium_HD_3prime_Human_Pancreatic_Cancer/`; do not package
this backup, its partial file, or the unused runtime dependencies. The HD
downloader now requires explicit opt-in and should not be run for this round.

The 1,817-byte PMC download response was HTML, not an XLSX workbook. It has been
renamed `metadata/PMC_supplement_blocked_response.html` and is explicitly
excluded from packaging. Only the validated publisher workbook named above
is used as supplementary-table evidence. All `*.partial` files, the entire
unused HD backup directory and unused `runtime_deps/` are excluded as well.
No new archive has been created in this source-preparation task.
