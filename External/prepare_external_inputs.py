"""Fetch and align public NYU_PDAC1_Vis inputs without computing any scores."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import sys
import time
import urllib.request

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree


HERE = Path(__file__).resolve().parent
SAMPLE = "GSM6177618_NYU_PDAC1_Vis"
PREFIX = SAMPLE + "_processed_"
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6177nnn/GSM6177618/suppl/"
SERIES_SOFT = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE203nnn/GSE203612/soft/GSE203612_family.soft.gz"
COORDINATE_DOC = "https://www.10xgenomics.com/support/software/space-ranger/2.1/analysis/outputs/spatial-outputs"
FILES = [
    PREFIX + "filtered_feature_bc_matrix.h5",
    PREFIX + "spatial_tissue_positions_list.csv.gz",
    PREFIX + "spatial_scalefactors_json.json.gz",
    PREFIX + "spatial_tissue_hires_image.png.gz",
    PREFIX + "spatial_tissue_lowres_image.png.gz",
    PREFIX + "spatial_aligned_fiducials.jpg.gz",
    PREFIX + "spatial_detected_tissue_image.jpg.gz",
]


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    status = "reused_existing"
    headers = {}
    if not destination.exists():
        temporary = destination.with_name(destination.name + ".partial")
        for attempt in range(3):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "CADSPOT-input-provenance/1.0"})
                with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as out:
                    headers = {k: response.headers.get(k) for k in ("Content-Length", "Last-Modified", "ETag", "Content-Type")}
                    shutil.copyfileobj(response, out)
                if headers.get("Content-Length") and temporary.stat().st_size != int(headers["Content-Length"]):
                    raise ValueError("Content-Length mismatch")
                temporary.replace(destination)
                status = "downloaded_from_official_source"
                break
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(2 * (attempt + 1))
    record = {
        "url": url, "file": str(destination.relative_to(HERE)),
        "bytes": destination.stat().st_size, "sha256": sha256(destination),
        "retrieved_or_verified_utc": datetime.now(timezone.utc).isoformat(),
        "status": status, "http_headers_this_run": headers,
    }
    print(json.dumps(record), flush=True)
    return record


def decode(values):
    return np.asarray([v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in values])


def extract_sample_soft(soft_path, output):
    with gzip.open(soft_path, "rt", encoding="utf-8") as handle:
        text = handle.read()
    selected = []
    active = False
    for line in text.splitlines():
        if line.startswith("^SAMPLE = "):
            active = line == "^SAMPLE = GSM6177618"
        elif line.startswith("^"):
            active = False
        if active:
            selected.append(line)
    if not selected:
        raise ValueError("GSM6177618 absent from official family SOFT")
    output.write_text("\n".join(selected) + "\n", encoding="utf-8")
    return [line for line in selected if any(key in line for key in
        ("!Sample_title =", "!Sample_source_name_ch1 =", "!Sample_characteristics_ch1 =", "!Sample_relation ="))]


def prepare(skip_documentation_fetch=False):
    import h5py
    from PIL import Image

    raw = HERE / "raw" / SAMPLE
    spatial = HERE / "spatial" / SAMPLE
    prepared = HERE / "prepared" / SAMPLE
    metadata = HERE / "metadata"
    for folder in (raw, spatial, prepared, metadata):
        folder.mkdir(parents=True, exist_ok=True)
    sources = [fetch(BASE, metadata / "GSM6177618_official_ftp_listing.html")]
    for name in FILES:
        sources.append(fetch(BASE + name, raw / name))
        if name.endswith(".gz"):
            short_name = name.removeprefix(PREFIX + "spatial_").removesuffix(".gz")
            with gzip.open(raw / name, "rb") as source, (spatial / short_name).open("wb") as out:
                shutil.copyfileobj(source, out)
    sources.append(fetch(SERIES_SOFT, metadata / "GSE203612_family.soft.gz"))
    sample_lines = extract_sample_soft(metadata / "GSE203612_family.soft.gz", metadata / "GSM6177618_official_sample.soft.txt")
    if not skip_documentation_fetch:
        try:
            sources.append(fetch(COORDINATE_DOC, metadata / "10x_spatial_outputs.html"))
        except Exception as error:
            sources.append({"url": COORDINATE_DOC, "status": "not_saved", "error": str(error)})

    h5path = raw / FILES[0]
    with h5py.File(h5path, "r") as handle:
        matrix = handle["matrix"]
        counts = sparse.csc_matrix((matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]),
                                   shape=tuple(matrix["shape"][:])).T.tocsr().astype(np.int64)
        barcodes = decode(matrix["barcodes"][:])
        symbols = decode(matrix["features"]["name"][:])
        feature_ids = decode(matrix["features"]["id"][:])
        feature_types = decode(matrix["features"]["feature_type"][:]) if "feature_type" in matrix["features"] else np.repeat("Gene Expression", len(symbols))
    if counts.data.size and (not np.isfinite(counts.data).all() or (counts.data < 0).any()):
        raise ValueError("Invalid count matrix values")
    if not pd.Index(barcodes).is_unique:
        raise ValueError("Duplicate matrix barcodes")
    keep = feature_types == "Gene Expression"
    if not np.all(symbols[keep] != ""):
        raise ValueError("Empty gene symbols must not be silently replaced")
    features = pd.DataFrame({"feature_id": feature_ids, "gene": symbols, "feature_type": feature_types, "retained_gene_expression": keep})
    features.to_csv(prepared / "source_features.tsv", sep="\t", index=False)
    counts = counts[:, keep]
    symbols = symbols[keep]
    # Preserve the first appearance of each symbol; combine duplicate symbol counts.
    genes = pd.Index(pd.unique(symbols))
    mapping = sparse.coo_matrix((np.ones(len(symbols), dtype=np.int64),
        (np.arange(len(symbols)), genes.get_indexer(symbols))), shape=(len(symbols), len(genes))).tocsr()
    original_totals = np.asarray(counts.sum(axis=1)).ravel()
    counts = (counts @ mapping).tocsr()
    counts.eliminate_zeros()
    counts.sort_indices()
    total = np.asarray(counts.sum(axis=1)).ravel()
    if not np.array_equal(total, original_totals):
        raise ValueError("Duplicate symbol aggregation changed library totals")
    sparse.save_npz(prepared / "counts.npz", counts)
    pd.DataFrame({"gene": genes}).to_csv(prepared / "genes.tsv", sep="\t", index=False)
    duplicate_map = pd.DataFrame({"feature_id": feature_ids[keep], "gene": symbols})
    duplicate_map.to_csv(prepared / "feature_to_gene.tsv", sep="\t", index=False)

    columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    positions = pd.read_csv(spatial / "tissue_positions_list.csv", header=None, names=columns)
    if not positions.barcode.is_unique:
        raise ValueError("Duplicate coordinate barcodes")
    positions[columns[1:]] = positions[columns[1:]].apply(pd.to_numeric, errors="raise")
    if not positions.in_tissue.isin([0, 1]).all():
        raise ValueError("Invalid in_tissue flags")
    positions["x"] = positions.pxl_col_in_fullres
    positions["y"] = positions.pxl_row_in_fullres
    positions["has_expression"] = positions.barcode.isin(barcodes)
    positions.to_csv(prepared / "positions_all.tsv", sep="\t", index=False)
    order = pd.Index(positions.barcode).get_indexer(barcodes)
    if (order < 0).any():
        raise ValueError("Matrix barcode missing from coordinate table")
    spots = positions.iloc[order].reset_index(drop=True)
    if not np.array_equal(spots.barcode.to_numpy(), barcodes):
        raise ValueError("Count-position row order mismatch")
    if not np.isfinite(spots[["x", "y"]].to_numpy()).all() or spots.duplicated(["x", "y"]).any():
        raise ValueError("Invalid or duplicate coordinates")
    spots["total_umi"] = total
    spots["n_genes"] = np.diff(counts.indptr)
    spots["passes_locked_qc"] = (spots.in_tissue == 1) & (spots.total_umi >= 100) & (spots.n_genes >= 100)
    scale = json.loads((spatial / "scalefactors_json.json").read_text(encoding="utf-8"))
    for label in ("hires", "lowres"):
        factor = scale[f"tissue_{label}_scalef"]
        spots[f"x_{label}"] = spots.x * factor
        spots[f"y_{label}"] = spots.y * factor
    spots.to_csv(prepared / "positions.tsv", sep="\t", index=False)
    pd.DataFrame({"sample_id": SAMPLE, "barcode": barcodes}).to_csv(prepared / "barcodes.tsv", sep="\t", index=False)

    image_audit = {}
    for label in ("hires", "lowres"):
        with Image.open(spatial / f"tissue_{label}_image.png") as image:
            image.load()
            width, height = image.size
        x, y = spots[f"x_{label}"], spots[f"y_{label}"]
        valid = (x >= 0) & (x < width) & (y >= 0) & (y < height)
        image_audit[label] = {"width": width, "height": height, "n_centers_in_bounds": int(valid.sum()), "all_in_bounds": bool(valid.all())}
        if not valid.all():
            raise ValueError("Registered coordinates fall outside image")
    nearest = cKDTree(spots[["x", "y"]].to_numpy()).query(spots[["x", "y"]].to_numpy(), k=2)[0][:, 1]
    local_original = Path(r"D:\BMC_Bioinformatics_CADSPOT_REBUILD_20260730\work\GSE203612_extracted") / FILES[0]
    local_check = {"path": str(local_original), "exists": local_original.exists()}
    if local_original.exists():
        local_check["sha256"] = sha256(local_original)
        local_check["same_as_fresh_official_download"] = local_check["sha256"] == sha256(h5path)
    manifest = {
        "sample_id": SAMPLE, "repository": "NCBI Gene Expression Omnibus (GEO)",
        "series_accession": "GSE203612", "sample_accession": "GSM6177618",
        "sample_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6177618",
        "series_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203612",
        "biosample_url": "https://www.ncbi.nlm.nih.gov/biosample/SAMN28600492",
        "original_paper_doi": "10.1038/s41588-022-01141-9",
        "identity_status": "CONFLICT: title/source say PDAC; characteristics say Ovary/Ovarian carcinoma; author-supplement check pending",
        "official_sample_metadata_lines": sample_lines,
        "independence": "Different study from GSE274103. This same external count matrix was already used for a technical-only exploratory audit; not an untouched prospective holdout.",
        "sources": sources, "local_original_count_check": local_check,
        "matrix": {"format": "scipy_CSR_npz", "orientation": "spots_by_genes", "value_type": "nonnegative_integer_deposited_Gene_Expression_counts", "n_spots": counts.shape[0], "n_unique_gene_symbols": counts.shape[1], "source_gene_expression_features": int(keep.sum()), "source_other_features_excluded": int((~keep).sum()), "duplicate_symbol_aggregation": "sum_counts_per_spot", "nnz": counts.nnz, "row_order": "identical_to_positions.tsv_and_barcodes.tsv", "column_order": "identical_to_genes.tsv"},
        "coordinates": {"units": "full_resolution_image_pixels", "x": "pxl_col_in_fullres", "y": "pxl_row_in_fullres", "origin": "top_left_image", "array_coordinates": "dimensionless_10x_hexagonal_array_indices_not_micrometers", "image_conversion": "x_image=x*tissue_image_scalef; y_image=y*tissue_image_scalef", "missing_coordinate_handling": "fail_if_any_expression_barcode_missing_not_imputed", "n_grid_positions": len(positions), "n_expression_positions": len(spots), "n_expression_out_of_tissue": int((spots.in_tissue == 0).sum()), "median_nearest_other_spot_distance_pixels": float(np.median(nearest)), "registered_image_audit": image_audit, "coordinate_documentation_url": COORDINATE_DOC},
        "qc": {"input": "deposited_filtered_matrix", "additional_filter_applied": False, "locked_rule_flag_only": "in_tissue==1 and total_umi>=100 and n_genes>=100", "n_passing_locked_rule": int(spots.passes_locked_qc.sum()), "total_umi_min_median_max": [int(total.min()), float(np.median(total)), int(total.max())], "n_genes_min_median_max": [int(spots.n_genes.min()), float(spots.n_genes.median()), int(spots.n_genes.max())]},
        "scope": "INPUT_PREPARATION_ONLY: no scores, composition proxies, candidate ranking, bootstrap, or random controls generated. No score-dependent sample or spot selection.",
        "composition_note": "Any later gene-module proxies are transcript-derived covariates, not validated cell fractions or malignant-cell calls.",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "script_sha256": sha256(__file__),
        "environment": {"python": sys.version, **{name: importlib.metadata.version(name) for name in ("numpy", "pandas", "scipy", "h5py", "pillow")}},
    }
    manifest["prepared_files"] = [{"file": str(path.relative_to(HERE)), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in sorted(prepared.iterdir()) if path.is_file()]
    (HERE / "source_manifest.json").write_text(json.dumps(manifest, ensure_ascii=True, indent=2) + "\n", encoding="utf-8")
    print("INPUT_AUDIT", json.dumps({"n_spots": len(spots), "n_genes": len(genes), "all_coordinates_present": True, "qc_pass": int(spots.passes_locked_qc.sum()), "official_metadata": sample_lines}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-documentation-fetch", action="store_true")
    args = parser.parse_args()
    prepare(skip_documentation_fetch=args.skip_documentation_fetch)
