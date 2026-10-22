#!/usr/bin/env python3
"""Reproduce the original four-coefficient reference-assisted RH calibration.

This is a full-cohort calibration, not independent validation. It reproduces
the original R v2 regression, not the manuscript's written two-term equation.
The independent inverse-calibration experiment has a separate entry point.
"""
import numpy as np
import pandas as pd

from src.common import PERCENTILES, RESULTS, load_cohort, score, write_json


def main():
    d = load_cohort()
    out = RESULTS/"rh_calibration"
    out.mkdir(parents=True, exist_ok=True)
    slope = 12.5*np.tan(np.deg2rad(d.ALS_Slope.to_numpy(float)))
    estimates = d[["site", "ID"]].copy()
    coefficients, rows = [], []
    max_error = 0.0
    for p in PERCENTILES:
        g = d[f"Gs_RH{p}"].to_numpy(float)
        h = d[f"ALS_RH{p}"].to_numpy(float)
        X = np.column_stack([np.ones(len(d)), slope, h, d.ALS_DEM])
        beta, _, rank, _ = np.linalg.lstsq(X, g-h, rcond=None)
        if rank != 4:
            raise ValueError(f"Rank-deficient calibration design for RH{p}.")
        corrected = g-X@beta
        estimates[f"recalibrated_RH{p}"] = corrected
        error = float(np.max(np.abs(corrected-d[f"Pred_RH{p}"].to_numpy(float))))
        max_error = max(max_error, error)
        coefficients.append({"RH": p, "intercept_d": beta[0], "slope_a": beta[1],
                             "reference_RH_b": beta[2], "elevation_c": beta[3]})
        for label, pred in [("raw", g), ("reference_assisted_calibration", corrected)]:
            rows.append({"RH": p, "input": label, **score(h, pred),
                         "max_absolute_difference_from_stored_Pred_RH": error})
    pd.DataFrame(coefficients).to_csv(out/"coefficients.csv", index=False)
    pd.DataFrame(rows).to_csv(out/"calibration_metrics.csv", index=False)
    estimates.to_csv(out/"recalibrated_RH.csv", index=False)
    write_json(out/"reproduction_check.json", {
        "equation": "G_p-H_p = d_p+a_p*(12.5*tan(slope))+b_p*H_p+c_p*DEM",
        "calibration_rows": len(d), "max_absolute_difference_from_stored_Pred_RH_m": max_error,
        "matches_stored_RH_within_1e-7_m": bool(max_error < 1e-7),
        "reference_RH_required_at_prediction_time": True,
        "independent_validation": False,
        "matches_manuscript_written_two_term_equation": False})
    if max_error >= 1e-7:
        raise AssertionError("The recalibrated RH columns do not reproduce the supplied Pred_RH columns.")
    print(f"Stored RH calibration reproduced; maximum difference = {max_error:.3g} m.")
    print("This calculation uses reference RH and is not independent validation.")


if __name__ == "__main__":
    main()
