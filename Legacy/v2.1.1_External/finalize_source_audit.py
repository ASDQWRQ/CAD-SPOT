"""Verify prepared inputs and record unresolved source identity; compute no scores."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import zipfile

import h5py
import numpy as np
import openpyxl
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parent
SAMPLE = "GSM6177618_NYU_PDAC1_Vis"
PAPER_URL = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-022-01141-9/MediaObjects/41588_2022_1141_MOESM6_ESM.xlsx"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(value, destination):
    destination.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def main():
    original_bad_name = ROOT / "metadata/Barkley2022_Supplementary_Tables.xlsx"
    if original_bad_name.exists() and not zipfile.is_zipfile(original_bad_name):
        original_bad_name.rename(ROOT / "metadata/PMC_supplement_blocked_response.html")
    workbook = ROOT / "metadata/Barkley2022_publisher_Supplementary_Tables.xlsx"
    if not zipfile.is_zipfile(workbook):
        raise ValueError("Publisher supplement is not a valid XLSX")
    book = openpyxl.load_workbook(workbook, read_only=True, data_only=True)
    rows = list(book["TableS1"].iter_rows(values_only=True))
    pathology = [(i, row) for i, row in enumerate(rows, 1) if str(row[0]).strip() == "PDAC NYU1"]
    if len(pathology) != 1:
        raise ValueError("Published PDAC NYU1 row missing or ambiguous")
    row_number, row = pathology[0]
    pd.DataFrame(rows).to_csv(ROOT / "metadata/Barkley2022_TableS1.tsv", sep="\t", header=False, index=False)
    evidence = {
        "publisher_supplement_url": PAPER_URL, "file": str(workbook.relative_to(ROOT)),
        "sha256": sha256(workbook), "worksheet": "TableS1", "excel_row": row_number,
        "sample_label": row[0], "pathology_value": row[1],
        "interpretation": "The published authors' sample table lists PDAC NYU1, but pathology is NA. This corroborates the PDAC sample label but does not resolve the inconsistent GEO/BioSample organ and cancer-type fields. No author correction or specimen-level pathology record was obtained.",
        "status": "repository_identity_conflict_unresolved",
    }
    write_json(evidence, ROOT / "metadata/NYU_PDAC1_identity_audit.json")

    prepared = ROOT / "prepared" / SAMPLE
    matrix = sparse.load_npz(prepared / "counts.npz").tocsr()
    genes = pd.read_csv(prepared / "genes.tsv", sep="\t")
    positions = pd.read_csv(prepared / "positions.tsv", sep="\t")
    source_h5 = ROOT / "raw" / SAMPLE / (SAMPLE + "_processed_filtered_feature_bc_matrix.h5")
    with h5py.File(source_h5, "r") as handle:
        source_barcodes = [value.decode() for value in handle["matrix/barcodes"][:]]
        source_dtype = str(handle["matrix/data"].dtype)
    manifest_path = ROOT / "source_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    tests = {
        "matrix_is_CSR": sparse.isspmatrix_csr(matrix),
        "matrix_shape_matches_axes": matrix.shape == (len(positions), len(genes)),
        "genes_unique": bool(genes.gene.is_unique),
        "barcodes_unique": bool(positions.barcode.is_unique),
        "raw_h5_barcode_order_matches_prepared_rows": source_barcodes == positions.barcode.tolist(),
        "sample_id_absent_for_prepare_section_insert": "sample_id" not in positions,
        "coordinates_finite": bool(np.isfinite(positions[["x", "y"]].to_numpy()).all()),
        "coordinates_unique": not bool(positions.duplicated(["x", "y"]).any()),
        "counts_nonnegative_integers": bool(np.issubdtype(matrix.dtype, np.integer) and (matrix.data >= 0).all()),
        "total_umi_matches_matrix": bool(np.array_equal(np.asarray(matrix.sum(axis=1)).ravel(), positions.total_umi.to_numpy())),
        "n_genes_matches_matrix": bool(np.array_equal(np.diff(matrix.indptr), positions.n_genes.to_numpy())),
        "qc_flag_matches_locked_rule": bool(np.array_equal(positions.passes_locked_qc.to_numpy(), ((positions.in_tissue == 1) & (positions.total_umi >= 100) & (positions.n_genes >= 100)).to_numpy())),
        "all_prepared_hashes_match": all(sha256(ROOT / record["file"]) == record["sha256"] for record in manifest["prepared_files"]),
        "no_score_columns_generated": not any("score" in column.lower() for column in positions),
        "published_workbook_valid": True,
    }
    if not all(tests.values()):
        raise ValueError(f"Input audit failed: {tests}")
    audit = {"sample_id": SAMPLE, "verified_utc": datetime.now(timezone.utc).isoformat(),
             "n_spots": len(positions), "n_unique_gene_symbols": len(genes),
             "n_pass_locked_qc": int(positions.passes_locked_qc.sum()),
             "source_h5_count_dtype": source_dtype, "tests": tests,
             "ready_for_technical_workflow": True, "confirmed_PDAC_validation_identity": False,
             "scope": "Input validation only; no scoring, model fitting, ranking or random controls."}
    write_json(audit, ROOT / "Input_Verification.json")
    manifest["identity_status"] = "Published TableS1 corroborates PDAC NYU1 label; pathology is NA; GEO/BioSample Ovary/Ovarian-carcinoma conflict remains unresolved. Technical-transfer use only, not confirmed PDAC validation."
    manifest["published_identity_evidence"] = evidence
    manifest["ready_for_technical_workflow"] = True
    manifest["confirmed_PDAC_validation_identity"] = False
    manifest["untouched_holdout"] = False
    manifest["permitted_external_claim"] = "Independent-study full-workflow technical transfer example with unresolved tumor identity; not a confirmed PDAC validation case or an untouched holdout."
    blocked_response = ROOT / "metadata/PMC_supplement_blocked_response.html"
    manifest["packaging_exclusions"] = {
        "patterns": ["**/*.partial", "raw/Visium_HD_3prime_Human_Pancreatic_Cancer/**", "runtime_deps/**"],
        "files": [{"file": str(blocked_response.relative_to(ROOT)), "bytes": blocked_response.stat().st_size,
                   "sha256": sha256(blocked_response), "reason": "HTML response masqueraded as XLSX; not supplementary evidence"}] if blocked_response.exists() else [],
        "new_archive_created": False,
    }
    manifest["analysis_plan"] = {"path": "../analysis_plan.md", "sha256": sha256(ROOT.parent / "analysis_plan.md")}
    manifest["technical_input_verification"] = "Input_Verification.json"
    write_json(manifest, manifest_path)
    hd_manifest_path = ROOT / "official_10x_download_manifest.json"
    if hd_manifest_path.exists():
        hd = json.loads(hd_manifest_path.read_text(encoding="utf-8"))
        hd["status"] = "STOPPED_AT_USER_REQUEST_NOT_INCLUDED_IN_CURRENT_ANALYSIS"
        hd["stopping_reason"] = "Explicit user instruction to avoid multi-GB HD download and extra tuning; not a demonstrated network blockade. Spatial download had a recoverable read timeout; binned download was actively terminated."
        hd["incomplete_files_do_not_package"] = [
            {"file": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256_of_partial_file_only": sha256(path), "is_complete_archive": False}
            for path in (ROOT / "raw/Visium_HD_3prime_Human_Pancreatic_Cancer").glob("*.partial")]
        hd["HD_inputs_prepared"] = False
        write_json(hd, hd_manifest_path)
    records = []
    for record in manifest["sources"]:
        if "file" in record:
            records.append({"file": record["file"], "url": record["url"], "bytes": record["bytes"], "sha256": record["sha256"], "scope": "NYU_official_source"})
    records.append({"file": evidence["file"], "url": PAPER_URL, "bytes": workbook.stat().st_size, "sha256": evidence["sha256"], "scope": "published_identity_audit"})
    pd.DataFrame(records).to_csv(ROOT / "source_checksums.tsv", sep="\t", index=False)
    print(json.dumps(audit, indent=2))
    print(json.dumps(evidence, indent=2))


if __name__ == "__main__":
    main()
