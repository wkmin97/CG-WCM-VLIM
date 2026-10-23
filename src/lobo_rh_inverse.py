#!/usr/bin/env python3
"""
Leave-one-UAV-block-out (LOBO) RH correction, one step only.

Fits the original four-coefficient v1/v2 *training regression*:
    G_p - H_p = d_p + a_p*S + b_p*H_p + c_p*Z
    S = 12.5 * tan(slope in degrees), Z = ALS_DEM
where G_p=Gs_RHp and H_p=ALS_RHp.

IMPORTANT: The original stored Pred_RH was fitted/reference-assisted.
For an independently held-out test, the test block's H_p cannot be an input.
This script uses the explicit algebraic inverse:
    Hhat_p = (G_p - d_p - a_p*S - c_p*Z) / (1 + b_p).
This prediction rule is NEW relative to the stored Pred_RH and should be
reported as "inverse-calibrated four-coefficient GC correction", not as
an unchanged implementation of the manuscript's two-term GC-WCM equation.

Outputs: held-out RH predictions, 12-block scores, pooled scores, coefficients.
Does NOT fit VLIM/AGBD, make figures, or modify existing prediction columns.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

PERCENTILES = list(range(5, 101, 5))


def score(ref, raw, corrected):
    """Common rows for a fair, paired raw-vs-corrected comparison."""
    y = np.asarray(ref, dtype=float)
    g = np.asarray(raw, dtype=float)
    p = np.asarray(corrected, dtype=float)
    keep = np.isfinite(y) & np.isfinite(g) & np.isfinite(p)
    y, g, p = y[keep], g[keep], p[keep]

    def calc(v):
        residual = v - y                  # prediction minus UAV reference
        rmse = np.sqrt(np.mean(residual**2))
        bias = np.mean(residual)
        r = (np.corrcoef(y, v)[0, 1]
             if len(y) > 1 and np.std(y) > 0 and np.std(v) > 0
             else np.nan)
        return r, rmse, bias

    if not len(y):
        return dict(n=0, raw_r=np.nan, raw_RMSE=np.nan, raw_Bias=np.nan,
                    corrected_r=np.nan, corrected_RMSE=np.nan,
                    corrected_Bias=np.nan, RMSE_reduction=np.nan)
    rr, re, rb = calc(g)
    cr, ce, cb = calc(p)
    return dict(n=len(y), raw_r=rr, raw_RMSE=re, raw_Bias=rb,
                corrected_r=cr, corrected_RMSE=ce, corrected_Bias=cb,
                RMSE_reduction=re - ce)


def run(input_csv, outdir, site_col="site", id_col="ID",
        slope_col="ALS_Slope", dem_col="ALS_DEM"):
    df = pd.read_csv(input_csv)
    need = [site_col, slope_col, dem_col]
    need += [f"{prefix}{p}" for p in PERCENTILES
             for prefix in ("Gs_RH", "ALS_RH")]
    absent = sorted(set(need).difference(df.columns))
    if absent:
        raise ValueError(f"Missing input columns: {absent}")

    if df[site_col].isna().any():
        raise ValueError("Some records lack a UAV block/site ID.")
    df[site_col] = df[site_col].astype(str)
    if id_col in df.columns:
        if df[[site_col, id_col]].isna().any().any():
            raise ValueError("Missing site+ID key detected.")
        if df.duplicated([site_col, id_col]).any():
            raise ValueError("Duplicate site+ID key detected.")

    sites = sorted(df[site_col].unique())
    if len(sites) != 12:
        raise ValueError(f"Expected exactly 12 UAV blocks; found {len(sites)}.")

    for name in set(need) - {site_col}:
        df[name] = pd.to_numeric(df[name], errors="raise")
        df[name] = df[name].replace([np.inf, -np.inf], np.nan)

    # From original R fitting definition; ALS_Slope MUST be in DEGREES.
    slope = df[slope_col].to_numpy(dtype=float)
    S = 12.5 * np.tan(np.deg2rad(slope))
    Z = df[dem_col].to_numpy(dtype=float)
    if np.any(np.isfinite(slope) & (np.abs(slope) >= 90)):
        raise ValueError("Slope must be in degrees and lie inside (-90,90).")

    site = df[site_col].to_numpy()
    fold_rows, coefficient_rows = [], []

    for p in PERCENTILES:
        g = df[f"Gs_RH{p}"].to_numpy(dtype=float)
        h = df[f"ALS_RH{p}"].to_numpy(dtype=float)
        pred = np.full(len(df), np.nan, dtype=float)

        train_valid = np.isfinite(g) & np.isfinite(h) & np.isfinite(S) & np.isfinite(Z)
        test_valid = np.isfinite(g) & np.isfinite(S) & np.isfinite(Z)

        for held_out_site in sites:
            train = (site != held_out_site) & train_valid
            test = (site == held_out_site) & test_valid
            if train.sum() < 5 or not test.any():
                raise ValueError(f"Too few valid rows: site={held_out_site}, RH{p}")

            # OLS equivalent to lm(G-H ~ S + ALS_RH + DEM) on 11 blocks.
            X = np.column_stack([np.ones(train.sum()), S[train], h[train], Z[train]])
            y = (g - h)[train]
            beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
            if rank != 4:
                raise ValueError(f"Rank-deficient training design: site={held_out_site}, RH{p}")

            d, a, b, c = beta
            if abs(1.0 + b) < 1e-8:
                raise ValueError(f"Cannot invert correction: 1+b≈0 for site={held_out_site}, RH{p}")

            # Never access h[test] in the prediction formula:
            pred[test] = (g[test] - d - a*S[test] - c*Z[test]) / (1.0 + b)

            metric = score(h[site == held_out_site],
                           g[site == held_out_site],
                           pred[site == held_out_site])
            fold_rows.append({"held_out_site": held_out_site, "RH": p,
                              "n_train": int(train.sum()), **metric})
            coefficient_rows.append({
                "held_out_site": held_out_site, "RH": p,
                "intercept_d": d, "slope_a": a,
                "reference_RH_b": b, "elevation_c": c,
                "inverse_denominator": 1 + b
            })

        if np.any(test_valid & ~np.isfinite(pred)):
            raise AssertionError(f"RH{p}: some eligible samples lacked OOF predictions.")
        out_col = f"LOBO_RH{p}"
        if out_col in df.columns:
            raise ValueError(f"Refusing to overwrite existing column {out_col}")
        df[out_col] = pred

    pooled_rows = []
    for p in PERCENTILES:
        pooled_rows.append({
            "RH": p,
            **score(df[f"ALS_RH{p}"], df[f"Gs_RH{p}"], df[f"LOBO_RH{p}"])
        })

    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    df.to_csv(output / "lobo_rh_predictions.csv", index=False)
    pd.DataFrame(fold_rows).to_csv(output / "lobo_rh_metrics_by_block.csv", index=False)
    pd.DataFrame(pooled_rows).to_csv(output / "lobo_rh_metrics_pooled.csv", index=False)
    pd.DataFrame(coefficient_rows).to_csv(output / "lobo_rh_coefficients.csv", index=False)

    print(f"LOBO done: {len(df)} footprints, {len(sites)} held-out blocks.")
    print(f"Results saved in: {output.resolve()}")
    print(pd.DataFrame(pooled_rows)[
        ["RH", "n", "raw_RMSE", "corrected_RMSE", "RMSE_reduction"]
    ].to_string(index=False, float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="data_pred_all_data_v2.csv")
    ap.add_argument("--outdir", default="lobo_rh_results")
    ap.add_argument("--site-col", default="site")
    ap.add_argument("--id-col", default="ID")
    ap.add_argument("--slope-col", default="ALS_Slope")
    ap.add_argument("--dem-col", default="ALS_DEM")
    args = ap.parse_args()
    run(args.input, args.outdir, args.site_col, args.id_col,
        args.slope_col, args.dem_col)
