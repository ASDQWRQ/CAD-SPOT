#!/usr/bin/env python
"""Score-table CAD-SPOT: associations, spatial rankings, and block uncertainty."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import SplineTransformer, StandardScaler

METHODS = ("score_residual_spearman", "partial_spearman", "spline_score_residual_spearman")
DERIVED = ("target_resid", "reference_resid", "target_resid_z", "reference_resid_z",
           "neighbor_reference_resid_z", "local_support", "priority_score",
           "priority_robust_z", "cadspot_candidate", "priority_rank", "spatial_status")


def read_table(path):
    """Read a TSV (optionally gzip compressed) without coercing identifiers."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=True)


def _numeric(df, columns, allow_missing=False):
    for col in dict.fromkeys(columns):
        if col not in df:
            raise ValueError(f"Missing column: {col}")
        v = pd.to_numeric(df[col], errors="coerce")
        bad = df[col].notna() & v.isna()
        if bad.any():
            raise ValueError(f"Non-numeric values in {col}: rows {np.flatnonzero(bad)[:5].tolist()}")
        if np.isinf(v.to_numpy(dtype=float)).any():
            raise ValueError(f"Infinite values in {col}")
        if not allow_missing and v.isna().any():
            raise ValueError(f"Missing numeric values in {col}; supply an explicit complete analysis table")
        df[col] = v.astype(float)


def validate_table(df, target, reference, technical, composition=(), sample_col="sample_id",
                   id_col="barcode", x_col="x", y_col="y", spatial=False):
    """Validate one row per (section, spot); missing coordinates remain explicit."""
    out = df.copy().reset_index(drop=True)
    covs = list(technical) + list(composition)
    if target == reference or target in covs or reference in covs or len(covs) != len(set(covs)):
        raise ValueError("Target/reference must differ, not be covariates; covariates must be unique")
    if out.empty:
        raise ValueError("Empty analysis table")
    for col in (sample_col, id_col):
        if col not in out or out[col].isna().any() or out[col].astype(str).str.strip().eq("").any():
            raise ValueError(f"Missing/empty identifier: {col}")
        out[col] = out[col].astype(str)
    if out.duplicated([sample_col, id_col]).any():
        raise ValueError("Duplicate spot IDs within section")
    if len(out) < 3:
        raise ValueError("At least three spots required")
    _numeric(out, [target, reference] + covs)
    if spatial:
        _numeric(out, [x_col, y_col], allow_missing=True)
    return out


def design(df, covariates, sample_col="sample_id", method="score_residual_spearman", spline_knots=4):
    """Fit covariate transformations on this table, including section fixed effects.

    The linear score design matches the original standardized OLS estimator on
    complete data. Partial Spearman rank-transforms numeric covariates only.
    Splines are additive cubic B-splines with quantile knots, not interactions.
    """
    x = df[list(covariates)].to_numpy(dtype=float)
    if method == "partial_spearman" and x.shape[1]:
        x = np.column_stack([stats.rankdata(v) for v in x.T])
    if method == "spline_score_residual_spearman" and x.shape[1]:
        # Constant covariates do not admit distinct quantile knots.
        pieces = []
        for col in x.T:
            nk = min(spline_knots, len(np.unique(col)))
            if nk >= 2:
                knots = np.unique(np.quantile(col, np.linspace(0, 1, nk)))[:, None]
                pieces.append(SplineTransformer(degree=3, knots=knots,
                                               include_bias=False).fit_transform(col[:, None]))
        x = np.column_stack(pieces) if pieces else np.empty((len(df), 0))
    if sample_col in df:
        dummies = pd.get_dummies(df[sample_col].astype(str), drop_first=True, dtype=float).to_numpy()
        x = np.column_stack([x, dummies])
    return StandardScaler().fit_transform(x) if x.shape[1] else np.empty((len(df), 0))


def resid(y, X):
    """OLS residuals with an intercept; no target/covariate imputation."""
    y = np.asarray(y, dtype=float)
    if not np.isfinite(y).all() or not np.isfinite(X).all():
        raise ValueError("Residualization requires finite scores and design")
    if not X.shape[1]:
        return y - y.mean(axis=0)
    return y - LinearRegression().fit(X, y).predict(X)


def correlation(a, b, rank=True):
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    if len(a) < 3 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return np.nan
    if rank:
        a, b = stats.rankdata(a), stats.rankdata(b)
    return float(np.corrcoef(a, b)[0, 1])


def fit_pair(df, target, reference, covariates, sample_col="sample_id",
             method="score_residual_spearman", spline_knots=4):
    """Single shared observed/random model entry point. Returns rho and residuals."""
    y = df[[target, reference]].to_numpy(dtype=float)
    if method == "partial_spearman":
        y = np.column_stack([stats.rankdata(y[:, i]) for i in range(2)])
    X = design(df, covariates, sample_col, method, spline_knots)
    r = resid(y, X)
    return correlation(r[:, 0], r[:, 1], rank=method != "partial_spearman"), r[:, 0], r[:, 1]


def _layers(technical, composition):
    result = [("raw", []) , ("technical_adjusted", list(technical))]
    if composition:
        result.append(("composition_adjusted", list(technical) + list(composition)))
    return result


def _association_rows(df, target, reference, technical, composition, sample_col,
                      methods, spline_knots, scope, section):
    rows = []
    for layer, covs in _layers(technical, composition):
        layer_methods = ("score_residual_spearman",) if layer == "raw" else methods
        for method in layer_methods:
            rho = (correlation(df[target], df[reference]) if layer == "raw" else
                   fit_pair(df, target, reference, covs, sample_col, method, spline_knots)[0])
            rows.append(dict(scope=scope, section=section, layer=layer, method=method,
                             spearman_rho=rho, n_spots=len(df),
                             status="ok" if np.isfinite(rho) else "undefined_constant_or_too_few"))
    return rows


def associations(df, target, reference, technical, composition=(), sample_col="sample_id",
                 methods=METHODS, spline_knots=4, per_section=True):
    """Raw Spearman; score-residual Spearman, partial Spearman, spline sensitivity."""
    rows = _association_rows(df, target, reference, technical, composition, sample_col,
                             methods, spline_knots, "pooled", "__pooled__")
    if per_section:
        for section, sub in df.groupby(sample_col, sort=True):
            rows += _association_rows(sub, target, reference, technical, composition, sample_col,
                                      methods, spline_knots, "section", section)
    return pd.DataFrame(rows)


def z_by_sample(df, col, sample_col):
    return df.groupby(sample_col)[col].transform(
        lambda s: (s - s.mean()) / (s.std(ddof=0) if s.std(ddof=0) > 0 else 1.0))


def positive_product(a, b):
    """Positive co-enrichment baseline: a negative-negative pair scores zero."""
    return np.maximum(a, 0) * np.maximum(b, 0)


def spatial_ranking(df, target_resid, reference_resid, sample_col="sample_id", id_col="barcode",
                    x_col="x", y_col="y", robust_z_cutoff=3.0, ranking_method="legacy",
                    target_scores=None, reference_scores=None):
    """Recompute exact six-neighbor within-section rankings, never reuse archives."""
    out = df.copy().reset_index(drop=True)
    for col in DERIVED:
        if col in out:
            archive = f"archived__{col}"
            if archive in out:
                raise ValueError(f"Archive-name collision: {archive}")
            out = out.rename(columns={col: archive})
    out["target_resid"], out["reference_resid"] = target_resid, reference_resid
    out["target_resid_z"] = z_by_sample(out, "target_resid", sample_col)
    out["reference_resid_z"] = z_by_sample(out, "reference_resid", sample_col)
    out["neighbor_reference_mean_z"] = np.nan
    out["neighbor_reference_resid_z"] = np.nan
    out["neighbor_ids"] = ""
    out["neighbor_count"] = 0
    valid = out[[x_col, y_col]].notna().all(axis=1)
    out["spatial_status"] = np.where(valid, "pending", "excluded_missing_coordinates")
    for section, group in out[valid].groupby(sample_col, sort=True):
        # Identifier sorting fixes distance ties independently of input row order.
        sub = group.sort_values(id_col, kind="stable")
        idx = sub.index.to_numpy()
        if len(sub) < 7:
            out.loc[idx, "spatial_status"] = "excluded_insufficient_neighbors"
            continue
        xy = sub[[x_col, y_col]].to_numpy(float)
        nn = NearestNeighbors(n_neighbors=7, metric="euclidean", algorithm="kd_tree").fit(xy)
        dist, near = nn.kneighbors(xy)
        rz = sub["reference_resid_z"].to_numpy()
        ids = sub[id_col].astype(str).to_numpy()
        for i, row_idx in enumerate(idx):
            # Include all ties at the seventh distance, then order by ID.
            candidates = nn.radius_neighbors(xy[i:i+1], radius=np.nextafter(dist[i, -1], np.inf),
                                            return_distance=False)[0]
            candidates = candidates[candidates != i]
            d2 = np.sum((xy[candidates] - xy[i]) ** 2, axis=1)
            order = np.lexsort((ids[candidates], d2))
            chosen = candidates[order[:6]]
            if len(chosen) != 6:
                raise RuntimeError("Internal error: failed exact six-neighbor construction")
            out.loc[row_idx, "neighbor_reference_mean_z"] = float(rz[chosen].mean())
            out.loc[row_idx, "neighbor_ids"] = json.dumps(ids[chosen].tolist())
            out.loc[row_idx, "neighbor_count"] = 6
        mean = out.loc[idx, "neighbor_reference_mean_z"]
        sd = mean.std(ddof=0)
        out.loc[idx, "neighbor_reference_resid_z"] = (mean - mean.mean()) / (sd if sd > 0 else 1.0)
        out.loc[idx, "spatial_status"] = "ok"
    out["local_support"] = out.target_resid_z * out.neighbor_reference_resid_z
    out["baseline_signed_product"] = out.target_resid_z * out.reference_resid_z
    out["baseline_positive_product"] = positive_product(out.target_resid_z, out.reference_resid_z)
    if target_scores is not None and reference_scores is not None:
        out["raw_target_score"] = np.asarray(target_scores, float)
        out["raw_reference_score"] = np.asarray(reference_scores, float)
        out["raw_target_z"] = z_by_sample(out, "raw_target_score", sample_col)
        out["raw_reference_z"] = z_by_sample(out, "raw_reference_score", sample_col)
        out["raw_signed_product"] = out.raw_target_z * out.raw_reference_z
        out["raw_positive_product"] = positive_product(out.raw_target_z, out.raw_reference_z)
    out["legacy_priority_score"] = out[["target_resid_z", "reference_resid_z",
                                            "neighbor_reference_resid_z", "local_support"]].mean(axis=1, skipna=False)
    out["priority_score"] = out.legacy_priority_score if ranking_method == "legacy" else (
        out.baseline_positive_product * np.maximum(out.neighbor_reference_resid_z, 0))
    ok = out.spatial_status.eq("ok")
    out.loc[~ok, "priority_score"] = np.nan
    out["priority_robust_z"] = np.nan
    out["priority_rank"] = np.nan
    out["baseline_positive_rank"] = np.nan
    out["positive_enrichment_eligible"] = ok & (out.target_resid_z > 0) & (out.reference_resid_z > 0) & (out.neighbor_reference_resid_z > 0)
    out["cadspot_candidate"] = False
    for _, sub in out[ok].groupby(sample_col, sort=True):
        s = sub.priority_score
        mad = float(np.median(np.abs(s - s.median()))) * 1.4826
        scale = mad if mad > 0 else (s.std(ddof=0) or 1.0)
        out.loc[sub.index, "priority_robust_z"] = (s - s.median()) / scale
        eligible = sub.index[out.loc[sub.index, "positive_enrichment_eligible"]]
        out.loc[eligible, "priority_rank"] = out.loc[eligible, "priority_score"].rank(ascending=False, method="min")
        baseline = sub.index[(sub.target_resid_z > 0) & (sub.reference_resid_z > 0)]
        out.loc[baseline, "baseline_positive_rank"] = out.loc[baseline, "baseline_positive_product"].rank(ascending=False, method="min")
    out["cadspot_candidate"] = out.positive_enrichment_eligible & (out.priority_robust_z >= robust_z_cutoff)
    out["ranking_method"] = ranking_method
    out["ranking_provenance"] = "recomputed_from_scores_not_archived_overrides"
    return out


def ranking_comparison(ranked, sample_col="sample_id", id_col="barcode", top_k=50,
                       top_fractions=(.01, .05, .1)):
    """Requested raw signed product AND positive-quadrant sensitivity, not accuracy.

    All comparisons use a common spatial-valid universe for both rankings.
    Positive-quadrant and positive-support restrictions apply to BOTH scores.
    The signed product itself remains exported for spots missing coordinates.
    """
    rows = []
    if any(not 0 < f <= 1 for f in top_fractions):
        raise ValueError("Top fractions must be in (0,1]")
    groups = [("pooled", "__pooled__", ranked)]
    groups += [("section", sec, sub) for sec, sub in ranked.groupby(sample_col, sort=True)]
    for scope, section, group in groups:
        sub = group.sort_values([sample_col, id_col], kind="stable")
        valid = sub[sub.spatial_status.eq("ok")]
        comparators = [("signed_residual_product", "baseline_signed_product", restriction,
                        sub.target_resid_z, sub.reference_resid_z)
                       for restriction in ("all_spots", "positive_quadrant", "cadspot_positive_support")]
        comparators += [("positive_part_residual_product", "baseline_positive_product", "positive_quadrant",
                         sub.target_resid_z, sub.reference_resid_z)]
        if "raw_signed_product" in sub:
            comparators += [("raw_score_z_signed_product", "raw_signed_product", "all_spots", sub.raw_target_z, sub.raw_reference_z),
                            ("raw_score_z_signed_product", "raw_signed_product", "positive_quadrant", sub.raw_target_z, sub.raw_reference_z)]
        for comparator, col, restriction, tz, rz in comparators:
            universe = valid
            if restriction != "all_spots":
                universe = universe[(tz.loc[universe.index] > 0) & (rz.loc[universe.index] > 0)]
            if restriction == "cadspot_positive_support":
                universe = universe[universe.positive_enrichment_eligible]
            selections = [("top_k", np.nan, top_k)] + [("top_fraction", f, max(1, int(np.ceil(f*len(universe))))) for f in top_fractions]
            for selection, fraction, k in selections:
                a = universe.sort_values("priority_score", ascending=False, kind="stable").head(k)
                b = universe.sort_values(col, ascending=False, kind="stable").head(k)
                aa = set(zip(a[sample_col], a[id_col]))
                bb = set(zip(b[sample_col], b[id_col]))
                overlap, union = len(aa & bb), aa | bb
                rows.append(dict(scope=scope, section=section, selection=selection, top_fraction=fraction,
                    top_k_requested=k, n_scope_spots=len(sub), n_spatial_valid=len(valid),
                    n_baseline_universe=len(universe), n_comparison_universe=len(universe),
                    n_missing_or_insufficient_coordinates=len(sub)-len(valid),
                    n_cadspot_top=len(aa), n_baseline_top=len(bb), overlap=overlap,
                    overlap_fraction=overlap/len(aa) if aa else np.nan,
                    overlap_fraction_of_cadspot=overlap/len(aa) if aa else np.nan,
                    overlap_fraction_of_baseline=overlap/len(bb) if bb else np.nan,
                    overlap_fraction_of_requested_k=overlap/k,
                    jaccard=overlap/len(union) if union else np.nan,
                    score_rank_spearman=correlation(universe.priority_score, universe[col]),
                    signed_product_negative_negative_count=int(((tz.loc[universe.index] < 0) & (rz.loc[universe.index] < 0)).sum()),
                    negative_negative_top_count=int(((tz.loc[b.index] < 0) & (rz.loc[b.index] < 0)).sum()),
                    comparator=comparator, restriction=restriction,
                    baseline_universe="common_spatial_valid_spots_same_restriction_for_both_scores",
                    top_fraction_denominator="n_comparison_universe",
                    cadspot_score_ranking="continuous_priority_score_no_hidden_candidate_filter",
                    inference="descriptive_not_validation"))
    return pd.DataFrame(rows)


def archived_candidate_comparison(ranked, sample_col="sample_id", id_col="barcode"):
    """Descriptive archived/fresh overlap; the archive never supplies fresh scores."""
    col = "archived__cadspot_candidate"
    if col not in ranked:
        return pd.DataFrame([dict(scope="pooled", section="__pooled__", status="no_archived_candidate_column")])
    text = ranked[col].astype(str).str.lower()
    if not text.isin(["true", "false", "1", "0"]).all():
        raise ValueError("Archived candidate column must contain explicit true/false or 1/0")
    old = text.isin(["true", "1"])
    rows = []
    groups = [("pooled", "__pooled__", ranked)]
    groups += [("section", sec, sub) for sec, sub in ranked.groupby(sample_col, sort=True)]
    for scope, section, sub in groups:
        a = set(zip(sub.loc[old.loc[sub.index], sample_col], sub.loc[old.loc[sub.index], id_col]))
        b = set(zip(sub.loc[sub.cadspot_candidate, sample_col], sub.loc[sub.cadspot_candidate, id_col]))
        rows.append(dict(scope=scope, section=section, n_archived_candidates=len(a), n_fresh_candidates=len(b),
                         overlap=len(a & b), archived_only=len(a-b), fresh_only=len(b-a),
                         jaccard=len(a & b)/len(a | b) if a | b else np.nan,
                         status="descriptive_archive_comparison_not_de_novo_reproduction",
                         archived_basis="saved_values_may_use_rank_residuals",
                         fresh_basis="linear_score_residuals_and_recomputed_six_neighbors"))
    return pd.DataFrame(rows)


def make_spatial_blocks(df, block_size, sample_col="sample_id", x_col="x", y_col="y"):
    """Non-overlapping square blocks, section-specific origin at observed min x/y."""
    if not np.isfinite(block_size) or block_size <= 0:
        raise ValueError("block_size must be finite and positive in declared coordinate units")
    out = df[df[[x_col, y_col]].notna().all(axis=1)].copy().reset_index(drop=True)
    if set(out[sample_col]) != set(df[sample_col]):
        raise ValueError("A section has no finite coordinates; fixed-section bootstrap cannot include it")
    out["block_x"], out["block_y"] = 0, 0
    for section, sub in out.groupby(sample_col, sort=True):
        xy = sub[[x_col, y_col]].to_numpy(float)
        blocks = np.floor((xy - xy.min(axis=0)) / block_size).astype(np.int64)
        out.loc[sub.index, ["block_x", "block_y"]] = blocks
        if len(np.unique(blocks, axis=0)) < 2:
            raise ValueError(f"Section {section} has fewer than two occupied blocks; reduce block size")
    return out


def sample_blocks(blocked, rng, sample_col="sample_id"):
    """Draw the original number of occupied blocks per fixed section with replacement."""
    parts, audit = [], []
    for section, sub in blocked.groupby(sample_col, sort=True):
        groups = [g for _, g in sub.groupby(["block_x", "block_y"], sort=True)]
        draws = rng.integers(0, len(groups), size=len(groups))
        part = pd.concat([groups[i] for i in draws], ignore_index=True)
        parts.append(part)
        audit.append(dict(section=section, n_blocks=len(groups), n_drawn_blocks=len(draws),
                          n_unique_drawn_blocks=len(np.unique(draws)), n_spots=len(part)))
    return pd.concat(parts, ignore_index=True), audit


def block_bootstrap(df, target, reference, technical, composition=(), sample_col="sample_id",
                    x_col="x", y_col="y", block_size=1000.0, n_boot=200, seed=20261115,
                    confidence=0.95, methods=("score_residual_spearman",), spline_knots=4):
    """Refit every model per draw; percentile CI conditional on observed sections.

    Missing-coordinate spots are excluded only here, not from score associations.
    Sampled block multiplicities are retained; spot counts may vary by replicate.
    """
    if n_boot < 2 or not 0 < confidence < 1:
        raise ValueError("n_boot >= 2 and 0 < confidence < 1 required")
    blocked = make_spatial_blocks(df, block_size, sample_col, x_col, y_col)
    point = associations(blocked, target, reference, technical, composition, sample_col, methods, spline_knots)
    rng = np.random.default_rng(seed)
    draws, audits = [], []
    for rep in range(n_boot):
        sampled, audit = sample_blocks(blocked, rng, sample_col)
        result = associations(sampled, target, reference, technical, composition, sample_col, methods, spline_knots)
        result["replicate"] = rep
        draws.append(result)
        audits.extend(dict(replicate=rep, **r) for r in audit)
    draws = pd.concat(draws, ignore_index=True)
    alpha = (1-confidence)/2
    keys = ["scope", "section", "layer", "method"]
    rows = []
    for _, estimate in point.iterrows():
        match = np.ones(len(draws), dtype=bool)
        for key in keys:
            match &= draws[key].eq(estimate[key]).to_numpy()
        vals = draws.loc[match, "spearman_rho"].dropna().to_numpy()
        enough = len(vals) >= max(2, int(np.ceil(0.8*n_boot)))
        row = estimate.to_dict()
        row.update(ci_low=float(np.quantile(vals, alpha)) if enough else np.nan,
                   ci_high=float(np.quantile(vals, 1-alpha)) if enough else np.nan,
                   confidence=confidence, n_boot_requested=n_boot, n_boot_valid=len(vals),
                   block_size=block_size, seed=seed, n_missing_coordinate_excluded=len(df)-len(blocked),
                   ci_status="ok" if enough else "insufficient_valid_resamples",
                   inference="conditional_on_observed_sections_not_population_of_sections",
                   residualization="refitted_in_each_resample")
        rows.append(row)
    return pd.DataFrame(rows), draws, pd.DataFrame(audits), blocked


def empirical_calibration(observed, null_rhos, alternative="greater"):
    """Plus-one gene-set empirical p. Never a spatial-dependence test."""
    values = np.asarray(null_rhos, dtype=float)
    if not np.isfinite(observed) or not len(values) or not np.isfinite(values).all():
        raise ValueError("Observed and all null statistics must be finite; no silent dropped controls")
    if alternative == "greater":
        count = np.sum(values >= observed)
    elif alternative == "two-sided":
        count = np.sum(np.abs(values) >= abs(observed))
    else:
        raise ValueError("alternative must be greater or two-sided")
    return dict(observed_rho=float(observed), n_controls=len(values), null_median=float(np.median(values)),
                null_95th_percentile=float(np.quantile(values, .95)), null_max=float(values.max()),
                plus_one_p=float((1+count)/(1+len(values))), alternative=alternative,
                inference="gene_set_empirical_calibration_not_spatial_inference")


def write_tsv(df, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def write_manifest(path, payload):
    import scipy
    import sklearn
    data = dict(cadspot_version="2.0.0", python=platform.python_version(), numpy=np.__version__,
                pandas=pd.__version__, scipy=scipy.__version__, sklearn=sklearn.__version__, **payload)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, indent=2, default=str) + "\n", encoding="utf-8")


def parser():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--technical", nargs="+", required=True)
    ap.add_argument("--composition", nargs="*", default=[])
    ap.add_argument("--sample-col", "--section-col", default="sample_id")
    ap.add_argument("--id-col", default="barcode")
    ap.add_argument("--x-col", default="x")
    ap.add_argument("--y-col", default="y")
    ap.add_argument("--coordinate-units", help="Shared x/y Euclidean units within each section; no axis mixing")
    ap.add_argument("--ranking-method", choices=("legacy", "positive_support"), default="legacy")
    ap.add_argument("--robust-z-cutoff", type=float, default=3.0)
    ap.add_argument("--top-k", type=int, default=50)
    ap.add_argument("--top-fractions", nargs="+", type=float, default=[.01, .05, .1])
    ap.add_argument("--spline-knots", type=int, default=4)
    ap.add_argument("--methods", nargs="+", choices=METHODS, default=list(METHODS))
    ap.add_argument("--coordinate-output")
    ap.add_argument("--ranking-comparison-output")
    ap.add_argument("--archived-comparison-output")
    ap.add_argument("--bootstrap-output")
    ap.add_argument("--bootstrap-draws-output")
    ap.add_argument("--bootstrap-audit-output")
    ap.add_argument("--block-map-output")
    ap.add_argument("--block-size", type=float)
    ap.add_argument("--bootstrap-reps", type=int, default=200)
    ap.add_argument("--bootstrap-methods", nargs="+", choices=METHODS, default=[METHODS[0]])
    ap.add_argument("--confidence", type=float, default=.95)
    ap.add_argument("--seed", type=int, default=20261115)
    ap.add_argument("--null-rho-file", help="Archived/external calibration only; no model-equivalence claim")
    ap.add_argument("--null-output")
    ap.add_argument("--manifest-output")
    ap.add_argument("--output", required=True)
    return ap


def main(argv=None):
    ap = parser()
    args = ap.parse_args(argv)
    spatial = bool(args.coordinate_output or args.ranking_comparison_output or args.archived_comparison_output or args.bootstrap_output)
    try:
        if spatial and not args.coordinate_units:
            raise ValueError("Spatial output requires --coordinate-units (same x/y units within section)")
        if args.bootstrap_output and args.block_size is None:
            raise ValueError("Bootstrap requires explicit --block-size in declared coordinate units")
        if any([args.bootstrap_draws_output, args.bootstrap_audit_output, args.block_map_output]) and not args.bootstrap_output:
            raise ValueError("Bootstrap auxiliary outputs require --bootstrap-output")
        if args.null_output and not args.null_rho_file:
            raise ValueError("--null-output requires --null-rho-file")
        if args.spline_knots < 2 or args.top_k < 1 or not np.isfinite(args.robust_z_cutoff):
            raise ValueError("Require spline-knots >=2, top-k >=1, finite robust-z-cutoff")
        df = validate_table(read_table(args.input), args.target, args.reference, args.technical,
                            args.composition, args.sample_col, args.id_col, args.x_col, args.y_col, spatial)
        summary = associations(df, args.target, args.reference, args.technical, args.composition,
                               args.sample_col, args.methods, args.spline_knots)
        write_tsv(summary, args.output)
        if args.coordinate_output or args.ranking_comparison_output or args.archived_comparison_output:
            covs = args.technical + args.composition
            _, a, b = fit_pair(df, args.target, args.reference, covs, args.sample_col)
            ranked = spatial_ranking(df, a, b, args.sample_col, args.id_col, args.x_col,
                                     args.y_col, args.robust_z_cutoff, args.ranking_method,
                                     df[args.target].to_numpy(), df[args.reference].to_numpy())
            ranked["coordinate_units"] = args.coordinate_units
            if args.coordinate_output:
                write_tsv(ranked, args.coordinate_output)
            if args.ranking_comparison_output:
                write_tsv(ranking_comparison(ranked, args.sample_col, args.id_col, args.top_k, args.top_fractions), args.ranking_comparison_output)
            if args.archived_comparison_output:
                write_tsv(archived_candidate_comparison(ranked, args.sample_col, args.id_col), args.archived_comparison_output)
        if args.bootstrap_output:
            ci, draws, audit, blocks = block_bootstrap(df, args.target, args.reference, args.technical,
                args.composition, args.sample_col, args.x_col, args.y_col, args.block_size,
                args.bootstrap_reps, args.seed, args.confidence, args.bootstrap_methods, args.spline_knots)
            ci["coordinate_units"] = args.coordinate_units
            write_tsv(ci, args.bootstrap_output)
            for frame, path in ((draws, args.bootstrap_draws_output), (audit, args.bootstrap_audit_output),
                                (blocks[[args.sample_col, args.id_col, "block_x", "block_y"]], args.block_map_output)):
                if path:
                    write_tsv(frame, path)
        if args.null_rho_file:
            null = read_table(args.null_rho_file)
            if "rho" not in null:
                raise ValueError("External null table requires named rho column; use cadspot_null.py for fresh generation")
            layer = "composition_adjusted" if args.composition else "technical_adjusted"
            rho = fit_pair(df, args.target, args.reference, args.technical + args.composition, args.sample_col)[0]
            cal = empirical_calibration(rho, pd.to_numeric(null.rho, errors="raise"))
            cal.update(layer=layer, provenance="external_archived_model_equivalence_unverified_not_fresh_reproduction")
            write_tsv(pd.DataFrame([cal]), args.null_output or str(args.output)+".external_null.tsv")
        write_manifest(args.manifest_output or str(args.output)+".manifest.json",
                       dict(arguments=vars(args), input_sha256=file_sha256(args.input), n_spots=len(df),
                            sections=sorted(df[args.sample_col].unique()),
                            primary_estimator="Spearman_of_linear_score_residuals_with_section_fixed_effects",
                            archived_overrides_used=False, spatial_neighbors=6,
                            coordinate_unit_assumption="User-declared common x/y Euclidean units; not inferable from numbers"))
        return 0
    except (ValueError, KeyError) as exc:
        ap.error(str(exc))


if __name__ == "__main__":
    sys.exit(main())
