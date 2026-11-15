#!/usr/bin/env python
"""Fresh matched gene-set generation and model-equivalent empirical calibration."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import maximum_bipartite_matching

import cadspot as core


def gene_list(path, gene_col="gene"):
    data = core.read_table(path)
    if gene_col not in data:
        raise ValueError(f"Gene list {path} requires a {gene_col} header")
    genes = data[gene_col]
    if genes.isna().any() or genes.astype(str).str.strip().eq("").any() or genes.duplicated().any():
        raise ValueError(f"Missing/duplicate gene in {path}")
    if not len(genes):
        raise ValueError(f"Empty gene list: {path}")
    return genes.astype(str).tolist()


def _spot_keys(spots, key_cols):
    for col in key_cols:
        if col not in spots or spots[col].isna().any() or spots[col].astype(str).str.strip().eq("").any():
            raise ValueError(f"Missing expression spot key: {col}")
    keys = pd.MultiIndex.from_frame(spots[list(key_cols)].astype(str))
    if keys.has_duplicates:
        raise ValueError("Expression/analysis spot keys must be unique; do not join samples by barcode alone")
    return keys


@dataclass
class ExpressionData:
    """Gene-by-spot matrix and ordered spot metadata; sparse input stays sparse."""
    values: object
    genes: list
    spots: pd.DataFrame

    def aligned(self, spots, key_cols=("sample_id", "barcode")):
        source, target = _spot_keys(self.spots, key_cols), _spot_keys(spots, key_cols)
        index = source.get_indexer(target)
        if (index < 0).any() or len(source) != len(target):
            raise ValueError("Expression and analysis must contain exactly the same keyed spots (order may differ)")
        return ExpressionData(self.values[:, index], self.genes, spots.copy().reset_index(drop=True))


def load_expression(path, spots_path, gene_col="gene", genes_path=None, orientation=None):
    """Read sparse NPZ with explicit orientation or a dense gene-by-spot TSV.

    The spots TSV is an ordered mapping of matrix spot axis to section/spot IDs.
    For dense TSV, optional expression_column verifies the actual header order.
    """
    spots = core.read_table(spots_path)
    if str(path).lower().endswith(".npz"):
        if not genes_path or orientation not in ("genes-by-spots", "spots-by-genes"):
            raise ValueError("NPZ requires --genes and explicit --matrix-orientation")
        genes = gene_list(genes_path, gene_col)
        values = sparse.load_npz(path)
        if orientation == "spots-by-genes":
            values = values.T
        values = values.tocsr()
    else:
        frame = core.read_table(path)
        if gene_col not in frame or frame[gene_col].isna().any() or frame[gene_col].duplicated().any():
            raise ValueError("Dense expression requires unique nonmissing gene identifiers")
        genes = frame.pop(gene_col).astype(str).tolist()
        if "expression_column" in spots and list(frame.columns) != spots.expression_column.astype(str).tolist():
            raise ValueError("Expression column headers do not match ordered expression_column map")
        values = frame.apply(pd.to_numeric, errors="raise").to_numpy(dtype=float)
    if values.shape != (len(genes), len(spots)):
        raise ValueError(f"Expression shape {values.shape} != ({len(genes)} genes, {len(spots)} spots)")
    data = values.data if sparse.issparse(values) else values
    if not np.isfinite(data).all():
        raise ValueError("Nonfinite expression values; no missing-gene/spot zero imputation")
    return ExpressionData(values, genes, spots)


def audit_gene_sets(expression, target, reference, intended_target=None, intended_reference=None):
    available = set(expression.genes)
    result = {}
    for label, used, intended in (("target", target, intended_target), ("reference", reference, intended_reference)):
        if not used or len(used) != len(set(used)):
            raise ValueError(f"{label} set must contain unique genes and be nonempty")
        missing = set(used) - available
        if missing:
            raise ValueError(f"Explicit used {label} genes missing in expression: {sorted(missing)}; provide available list, never zero-impute")
        if intended is not None and not set(used) <= set(intended):
            raise ValueError(f"Used {label} genes are not a subset of intended genes")
        result[label] = dict(used_genes=list(used), n_used=len(used), intended_genes=intended,
            n_intended=len(intended) if intended is not None else None,
            missing_from_expression=sorted(set(intended or [])-available),
            intended_present_but_not_used=sorted((set(intended or []) & available)-set(used)),
            missing_gene_policy="explicit_available_set_no_zero_imputation")
    return result


def score_gene_set(expression, genes, scoring="mean_gene_z", sample_col="sample_id"):
    """Arithmetic mean of expression or gene z-scores, same function for all sets.

    Input values are used as supplied (e.g. log1p CP10K); there is no implicit
    count normalization. Gene z-scores use population SD (ddof=0).
    """
    if not genes or len(genes) != len(set(genes)):
        raise ValueError("A scored gene set must be nonempty and unique")
    lookup = {g: i for i, g in enumerate(expression.genes)}
    missing = set(genes)-set(lookup)
    if missing:
        raise ValueError(f"Missing scored genes: {sorted(missing)}; no zero imputation")
    selected = expression.values[[lookup[g] for g in genes], :]
    values = selected.toarray().astype(float) if sparse.issparse(selected) else np.asarray(selected, float).copy()
    if scoring == "mean_expression":
        return values.mean(axis=0)
    if scoring == "mean_gene_z":
        groups = [np.arange(values.shape[1])]
    elif scoring == "mean_gene_z_by_section":
        if sample_col not in expression.spots:
            raise ValueError(f"Missing section column for gene-z scoring: {sample_col}")
        groups = [np.asarray(v) for v in expression.spots.groupby(sample_col, sort=True).indices.values()]
    else:
        raise ValueError(f"Unknown scoring rule: {scoring}")
    for group in groups:
        x = values[:, group]
        sd = x.std(axis=1, ddof=0)
        if (sd <= 0).any():
            bad = [g for g, s in zip(genes, sd) if s <= 0]
            raise ValueError(f"Constant gene(s) cannot be z-standardized under {scoring}: {bad}")
        values[:, group] = (x - x.mean(axis=1)[:, None]) / sd[:, None]
    return values.mean(axis=0)


def generate_matched_sets(expression, metadata, target_genes, reference_genes, universe,
                           exclusions=(), n_sets=500, seed=20261115, mean_tolerance=.25,
                           detection_tolerance=.05, length_tolerance=.25, gene_col="gene",
                           mean_col="mean_expression", detection_col="detection_rate", length_col="gene_length"):
    """Randomized maximum bipartite matches; unique genes per set, strict calipers.

    Match each target slot on log1p(mean), detection probability, and log(length).
    Target/reference and explicit exclusions are never selected. No relaxed or
    nearest-available fallback is permitted. This is not uniform set sampling.
    """
    audit_gene_sets(expression, target_genes, reference_genes)
    tolerances = np.array([mean_tolerance, detection_tolerance, length_tolerance], float)
    if not np.isfinite(tolerances).all() or (tolerances < 0).any() or n_sets < 1:
        raise ValueError("Require finite nonnegative tolerances and n_sets >= 1")
    if len(universe) != len(set(universe)) or not universe:
        raise ValueError("Universe must be explicit, unique, and nonempty")
    if gene_col not in metadata or metadata[gene_col].isna().any() or metadata[gene_col].duplicated().any():
        raise ValueError("Metadata must have unique nonmissing gene IDs")
    meta = metadata.set_index(gene_col).copy()
    required = set(universe) | set(target_genes) | set(reference_genes)
    if required-set(expression.genes) or required-set(meta.index):
        raise ValueError(f"Universe/observed genes must exist in expression and metadata; missing expression={sorted(required-set(expression.genes))[:10]}, metadata={sorted(required-set(meta.index))[:10]}")
    cols = [mean_col, detection_col, length_col]
    for col in cols:
        if col not in meta:
            raise ValueError(f"Missing matching metadata column {col}")
        meta[col] = pd.to_numeric(meta[col], errors="raise")
    check = meta.loc[sorted(required), cols].to_numpy(float)
    if not np.isfinite(check).all() or (check[:, 0] < 0).any() or ((check[:, 1] < 0) | (check[:, 1] > 1)).any() or (check[:, 2] <= 0).any():
        raise ValueError("Invalid matching metadata: require mean>=0, detection in [0,1], length>0, all finite")
    excluded = set(exclusions) | set(target_genes) | set(reference_genes)
    pool = sorted(set(universe)-excluded)
    if len(pool) < len(target_genes):
        raise ValueError("Too few eligible genes after explicit exclusions")
    def transformed(genes):
        x = meta.loc[genes, cols].to_numpy(float)
        x[:, 0], x[:, 2] = np.log1p(x[:, 0]), np.log(x[:, 2])
        return x
    targets, candidates = transformed(target_genes), transformed(pool)
    delta = candidates[None, :, :] - targets[:, None, :]
    eligible = np.all(np.abs(delta) <= tolerances[None, None, :] + 1e-12, axis=2)
    no_match = [g for g, row in zip(target_genes, eligible) if not row.any()]
    if no_match:
        raise ValueError(f"No candidate within declared calipers for: {no_match}; tolerances not relaxed")
    rng, rows = np.random.default_rng(seed), []
    for set_id in range(n_sets):
        row_order, col_order = rng.permutation(len(target_genes)), rng.permutation(len(pool))
        graph = sparse.csr_matrix(eligible[row_order][:, col_order].astype(np.int8))
        match = maximum_bipartite_matching(graph, perm_type="column")
        if (match < 0).any():
            raise ValueError("No unique-gene full matching exists under declared calipers; tolerances not relaxed")
        for perm_row, perm_col in enumerate(match):
            slot, chosen = int(row_order[perm_row]), int(col_order[perm_col])
            d = delta[slot, chosen]
            rows.append(dict(set_id=set_id, slot=slot, matched_to=target_genes[slot], gene=pool[chosen], seed=seed,
                             delta_log1p_mean=float(d[0]), delta_detection=float(d[1]), delta_log_length=float(d[2]),
                             candidate_count=int(eligible[slot].sum())))
    return pd.DataFrame(rows).sort_values(["set_id", "slot"]).reset_index(drop=True)


def calibrate_gene_sets(expression, spots, target_genes, reference_genes, matched_sets, technical,
                         composition=(), sample_col="sample_id", id_col="barcode", key_cols=("sample_id", "barcode"),
                         scoring="mean_gene_z", method="score_residual_spearman", spline_knots=4,
                         alternative="greater"):
    """Recompute observed AND every random target with identical scoring/model.

    The reference is fixed. This statistic is separate from an archived score
    table's .199 demonstration unless independently shown equivalent.
    """
    expression = expression.aligned(spots, key_cols)
    audit_gene_sets(expression, target_genes, reference_genes)
    work = spots.copy().reset_index(drop=True)
    target, reference = "fresh_target_score", "fresh_reference_score"
    if target in work or reference in work:
        raise ValueError("Input spot table already contains fresh score output names")
    work[target] = score_gene_set(expression, target_genes, scoring, sample_col)
    work[reference] = score_gene_set(expression, reference_genes, scoring, sample_col)
    work = core.validate_table(work, target, reference, technical, composition, sample_col, id_col)
    covs = list(technical) + list(composition)
    observed = core.fit_pair(work, target, reference, covs, sample_col, method, spline_knots)[0]
    required = {"set_id", "slot", "matched_to", "gene", "seed"}
    if not required <= set(matched_sets.columns):
        raise ValueError(f"Matched sets require {sorted(required)}")
    controls = []
    for set_id, group in matched_sets.groupby("set_id", sort=True):
        genes = group.gene.astype(str).tolist()
        if (len(genes) != len(target_genes) or len(set(genes)) != len(genes) or
            set(group.matched_to) != set(target_genes) or group.slot.duplicated().any() or
            set(genes) & (set(target_genes) | set(reference_genes))):
            raise ValueError(f"Invalid matched set {set_id}: size, slots, uniqueness, observed-set exclusion")
        trial = work[list(dict.fromkeys([sample_col, target, reference] + covs))].copy()
        trial[target] = score_gene_set(expression, genes, scoring, sample_col)
        rho = core.fit_pair(trial, target, reference, covs, sample_col, method, spline_knots)[0]
        controls.append(dict(set_id=set_id, rho=rho, n_genes=len(genes), scoring=scoring, method=method))
    null = pd.DataFrame(controls)
    if null.empty:
        raise ValueError("No random sets provided")
    summary = core.empirical_calibration(observed, null.rho.to_numpy(), alternative)
    summary.update(n_target_genes=len(target_genes), n_reference_genes=len(reference_genes),
                   scoring=scoring, method=method, n_spots=len(work),
                   adjustment="technical_plus_composition" if composition else "technical",
                   observed_random_model_equivalent=True,
                   score_provenance="fresh_expression_available_genes_not_archived_score_reproduction")
    return pd.DataFrame([summary]), null, work


def _expression_args(ap):
    ap.add_argument("--expression", required=True, help="Sparse scipy NPZ or gene-by-spot TSV")
    ap.add_argument("--expression-spots", required=True, help="Ordered matrix spot-axis TSV")
    ap.add_argument("--genes", help="Ordered NPZ gene-axis TSV")
    ap.add_argument("--matrix-orientation", choices=("genes-by-spots", "spots-by-genes"))
    ap.add_argument("--gene-col", default="gene")
    ap.add_argument("--spot-key-cols", nargs="+", default=["sample_id", "barcode"])
    ap.add_argument("--expression-scale", required=True, help="Audit label, e.g. log1p_CP10K; no implicit normalization")
    ap.add_argument("--target-genes", required=True, help="Explicit available gene list; no missing-gene filling")
    ap.add_argument("--reference-genes", required=True)
    ap.add_argument("--intended-target-genes")
    ap.add_argument("--intended-reference-genes")


def parser():
    ap = argparse.ArgumentParser(description=__doc__)
    commands = ap.add_subparsers(dest="command", required=True)
    gen = commands.add_parser("generate", help="Strict mean/detection/length matching; no calibration yet")
    _expression_args(gen)
    gen.add_argument("--metadata", required=True)
    gen.add_argument("--universe", required=True, help="Explicit gene TSV, not an implicit intersection")
    gen.add_argument("--exclusions", help="Additional excluded gene TSV; observed sets always excluded")
    gen.add_argument("--mean-col", default="mean_expression")
    gen.add_argument("--detection-col", default="detection_rate")
    gen.add_argument("--length-col", default="gene_length")
    gen.add_argument("--mean-tolerance", type=float, required=True, help="Absolute caliper on log1p(mean)")
    gen.add_argument("--detection-tolerance", type=float, required=True)
    gen.add_argument("--length-tolerance", type=float, required=True, help="Absolute caliper on natural log(length)")
    gen.add_argument("--n-sets", type=int, default=500)
    gen.add_argument("--seed", type=int, required=True)
    gen.add_argument("--output", required=True)
    cal = commands.add_parser("calibrate", help="Same fresh score/residual model for observed and random target")
    _expression_args(cal)
    cal.add_argument("--spot-table", required=True)
    cal.add_argument("--matched-sets", required=True)
    cal.add_argument("--matching-manifest", help="Default: MATCHED_SETS.manifest.json; required provenance")
    cal.add_argument("--technical", nargs="+", required=True)
    cal.add_argument("--composition", nargs="*", default=[])
    cal.add_argument("--sample-col", default="sample_id")
    cal.add_argument("--id-col", default="barcode")
    cal.add_argument("--scoring", choices=("mean_gene_z", "mean_gene_z_by_section", "mean_expression"), default="mean_gene_z")
    cal.add_argument("--method", choices=core.METHODS, default="score_residual_spearman")
    cal.add_argument("--spline-knots", type=int, default=4)
    cal.add_argument("--alternative", choices=("greater", "two-sided"), default="greater")
    cal.add_argument("--null-output", required=True)
    cal.add_argument("--score-output")
    cal.add_argument("--output", required=True)
    return ap


def main(argv=None):
    ap = parser()
    args = ap.parse_args(argv)
    try:
        expression = load_expression(args.expression, args.expression_spots, args.gene_col, args.genes, args.matrix_orientation)
        target, reference = gene_list(args.target_genes, args.gene_col), gene_list(args.reference_genes, args.gene_col)
        intended_t = gene_list(args.intended_target_genes, args.gene_col) if args.intended_target_genes else None
        intended_r = gene_list(args.intended_reference_genes, args.gene_col) if args.intended_reference_genes else None
        audit = audit_gene_sets(expression, target, reference, intended_t, intended_r)
        _spot_keys(expression.spots, args.spot_key_cols)
        paths = [args.expression, args.expression_spots, args.target_genes, args.reference_genes]
        paths += [p for p in (args.genes, args.intended_target_genes, args.intended_reference_genes) if p]
        if args.command == "generate":
            universe = gene_list(args.universe, args.gene_col)
            exclusions = gene_list(args.exclusions, args.gene_col) if args.exclusions else []
            matched = generate_matched_sets(expression, core.read_table(args.metadata), target, reference,
                universe, exclusions, args.n_sets, args.seed, args.mean_tolerance, args.detection_tolerance,
                args.length_tolerance, args.gene_col, args.mean_col, args.detection_col, args.length_col)
            core.write_tsv(matched, args.output)
            paths += [args.metadata, args.universe] + ([args.exclusions] if args.exclusions else [])
            payload = dict(arguments=vars(args), observed_gene_audit=audit,
                input_sha256={str(p): core.file_sha256(p) for p in paths},
                expression_sha256=core.file_sha256(args.expression), matched_sets_sha256=core.file_sha256(args.output),
                universe_genes=universe, additional_exclusions=exclusions,
                matching="randomized_maximum_bipartite_matching_not_uniform_over_all_feasible_sets",
                fixed_reference_genes=reference, seed=args.seed, tolerance_relaxation=False,
                unique_genes_per_set=True, repeated_sets_across_draws_permitted=True,
                null_target_set_size=len(target), reference_set_size=len(reference))
        else:
            provenance_path = args.matching_manifest or str(args.matched_sets)+".manifest.json"
            manifest = json.loads(Path(provenance_path).read_text(encoding="utf-8"))
            if manifest.get("expression_sha256") != core.file_sha256(args.expression):
                raise ValueError("Expression differs from matching manifest")
            for key in ("genes", "expression_spots"):
                original = manifest["arguments"].get(key)
                current = getattr(args, key)
                if bool(original) != bool(current) or (original and manifest["input_sha256"][original] != core.file_sha256(current)):
                    raise ValueError(f"{key} axis mapping differs from generation manifest")
            for key in ("matrix_orientation", "expression_scale", "gene_col", "spot_key_cols"):
                prior = manifest["arguments"].get(key, ["sample_id", "barcode"] if key == "spot_key_cols" else None)
                if prior != getattr(args, key):
                    raise ValueError(f"{key} differs from generation manifest")
            if manifest.get("matched_sets_sha256") != core.file_sha256(args.matched_sets):
                raise ValueError("Matched sets differ from generation manifest")
            if (manifest["observed_gene_audit"]["target"]["used_genes"] != target or
                manifest["observed_gene_audit"]["reference"]["used_genes"] != reference):
                raise ValueError("Observed gene sets differ from matching manifest")
            matched = core.read_table(args.matched_sets)
            if not set(matched.gene) <= (set(manifest["universe_genes"])-set(manifest["additional_exclusions"])-set(target)-set(reference)):
                raise ValueError("Random sets violate explicit universe/exclusions")
            summary, null, scores = calibrate_gene_sets(expression, core.read_table(args.spot_table), target,
                reference, matched, args.technical, args.composition, args.sample_col, args.id_col,
                args.spot_key_cols, args.scoring, args.method, args.spline_knots, args.alternative)
            core.write_tsv(summary, args.output)
            core.write_tsv(null, args.null_output)
            if args.score_output:
                core.write_tsv(scores, args.score_output)
            paths += [args.spot_table, args.matched_sets, provenance_path]
            payload = dict(arguments=vars(args), observed_gene_audit=audit,
                input_sha256={str(p): core.file_sha256(p) for p in paths}, generation_manifest=manifest,
                observed_random_model_equivalent=True, spatial_inference=False,
                generation_spot_key_default_inferred="spot_key_cols" not in manifest["arguments"],
                archived_score_table_equivalence="not_assumed_or_claimed")
        core.write_manifest(str(args.output)+".manifest.json", payload)
        return 0
    except (ValueError, KeyError, FileNotFoundError) as exc:
        ap.error(str(exc))


if __name__ == "__main__":
    sys.exit(main())
