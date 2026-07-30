#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

def design(df, covariates, sample_col):
    x = df[covariates].copy().apply(pd.to_numeric, errors="coerce")
    x = x.fillna(x.median(numeric_only=True))
    if sample_col in df.columns:
        x = pd.concat([x, pd.get_dummies(df[sample_col].astype(str), drop_first=True)], axis=1)
    return StandardScaler().fit_transform(x.values)

def resid(y, X):
    y = np.asarray(y, float)
    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    out = np.full_like(y, np.nan, dtype=float)
    model = LinearRegression().fit(X[m], y[m])
    out[m] = y[m] - model.predict(X[m])
    return out

def z_by_sample(df, col, sample_col):
    return df.groupby(sample_col)[col].transform(lambda s: (s - s.mean()) / (s.std(ddof=0) or 1))

def add_coordinates(df, target_resid, reference_resid, sample_col, x_col, y_col, robust_z_cutoff):
    out = df.copy()
    out["target_resid"] = target_resid
    out["reference_resid"] = reference_resid
    out["target_resid_z"] = z_by_sample(out, "target_resid", sample_col)
    out["reference_resid_z"] = z_by_sample(out, "reference_resid", sample_col)
    out["neighbor_reference_resid_z"] = 0.0
    out["local_support"] = 0.0
    if x_col in out.columns and y_col in out.columns:
        for sample, idx in out.groupby(sample_col).groups.items():
            sub = out.loc[idx]
            xy = sub[[x_col, y_col]].apply(pd.to_numeric, errors="coerce").values
            if len(sub) < 8 or np.isnan(xy).any():
                continue
            nn = NearestNeighbors(n_neighbors=7).fit(xy)
            neigh = nn.kneighbors(xy, return_distance=False)[:, 1:]
            rz = sub["reference_resid_z"].values
            tz = sub["target_resid_z"].values
            out.loc[idx, "neighbor_reference_resid_z"] = np.nanmean(rz[neigh], axis=1)
            out.loc[idx, "local_support"] = np.nanmean(np.minimum(tz[neigh], rz[neigh]), axis=1)
    for c in ["neighbor_reference_resid_z", "local_support"]:
        out[c] = z_by_sample(out, c, sample_col)
    out["priority_score"] = out[["target_resid_z", "reference_resid_z", "neighbor_reference_resid_z", "local_support"]].mean(axis=1)
    out["priority_robust_z"] = np.nan
    out["cadspot_candidate"] = False
    for sample, idx in out.groupby(sample_col).groups.items():
        s = out.loc[idx, "priority_score"]
        med = s.median()
        mad = np.nanmedian(np.abs(s - med)) * 1.4826
        if not np.isfinite(mad) or mad == 0:
            mad = s.std(ddof=0) or 1.0
        rz = (s - med) / mad
        out.loc[idx, "priority_robust_z"] = rz
        out.loc[idx, "cadspot_candidate"] = (rz >= robust_z_cutoff) & (out.loc[idx, "target_resid_z"] > 0) & (out.loc[idx, "reference_resid_z"] > 0) & (out.loc[idx, "neighbor_reference_resid_z"] > 0) & (out.loc[idx, "local_support"] > 0)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--technical", nargs="+", required=True)
    ap.add_argument("--composition", nargs="*", default=[])
    ap.add_argument("--sample-col", default="sample_id")
    ap.add_argument("--x-col", default="x")
    ap.add_argument("--y-col", default="y")
    ap.add_argument("--robust-z-cutoff", type=float, default=3.0)
    ap.add_argument("--coordinate-output", default=None)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    df = pd.read_csv(args.input, sep="\t")
    layers = []
    layers.append(("raw", df[args.target].values, df[args.reference].values))
    Xt = design(df, args.technical, args.sample_col)
    layers.append(("technical_adjusted", resid(df[args.target], Xt), resid(df[args.reference], Xt)))
    if args.composition:
        Xc = design(df, args.technical + args.composition, args.sample_col)
        comp_target = resid(df[args.target], Xc)
        comp_ref = resid(df[args.reference], Xc)
        layers.append(("composition_adjusted", comp_target, comp_ref))
    rows = []
    for name, a, b in layers:
        m = np.isfinite(a) & np.isfinite(b)
        rows.append({"layer": name, "spearman_rho": stats.spearmanr(a[m], b[m]).statistic, "n_spots": int(m.sum())})
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    if args.coordinate_output and args.composition:
        coord = add_coordinates(df, comp_target, comp_ref, args.sample_col, args.x_col, args.y_col, args.robust_z_cutoff)
        coord.to_csv(args.coordinate_output, sep="\t", index=False)

if __name__ == "__main__":
    main()
