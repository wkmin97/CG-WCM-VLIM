#!/usr/bin/env python3
"""Independently recompute exported scores and check prediction/data alignment."""
import json
import numpy as np
import pandas as pd

from src.common import ROOT, RESULTS, load_cohort, score, write_json


def main():
    source = load_cohort()
    out = RESULTS/"agbd_conditional_LOBO"
    pred = pd.read_csv(out/"LOBO_AGBD_all_predictions.csv", dtype={"site": str, "ID": str})
    if len(pred) != len(source) or pred.duplicated(["site", "ID"]).any():
        raise AssertionError("Predictions must contain exactly one row per source sample.")
    paired = source[["site", "ID", "ALS_AGB"]].merge(pred, on=["site", "ID"], suffixes=("_source", ""), validate="one_to_one")
    np.testing.assert_array_equal(paired.ALS_AGB_source, paired.ALS_AGB)
    metrics = pd.read_csv(out/"LOBO_AGBD_all_metrics_pooled.csv")
    supplied = pd.read_csv(ROOT/"data/reported/LOBO_AGBD_all_metrics_pooled.csv")
    source_comparison = []
    for prefix in ["Gs_RH", "Pred_RH"]:
        col = f"LOBO_layer_{prefix}"
        if not np.isfinite(paired[col]).all():
            raise AssertionError(f"Nonfinite held-out predictions: {col}")
        actual = score(paired.ALS_AGB_source, paired[col])
        row = metrics[metrics.model.eq("layer") & metrics.RH_input.eq(prefix)].iloc[0]
        for name in ["r", "RMSE", "Bias", "R2"]:
            np.testing.assert_allclose(actual[name], row[name], rtol=1e-9, atol=1e-9)
        expected = supplied[supplied.model.eq("layer") & supplied.RH_input.eq(prefix)].iloc[0]
        source_comparison.append({"RH_input": prefix, "recomputed_RMSE": actual["RMSE"],
                                  "supplied_RMSE": float(expected.RMSE),
                                  "difference": actual["RMSE"]-float(expected.RMSE)})
    fold_scores = pd.read_csv(out/"LOBO_AGBD_all_metrics_by_block.csv", dtype={"held_out_site": str})
    layer_scores = fold_scores[fold_scores.model.eq("layer")]
    for prefix in ["Gs_RH", "Pred_RH"]:
        scores = layer_scores[layer_scores.RH_input.eq(prefix)]
        if len(scores) != 12 or scores.held_out_site.nunique() != 12:
            raise AssertionError("Expected one layer-model score per held-out block and RH input.")
        for _, row in scores.iterrows():
            n_test = int(source.site.eq(row.held_out_site).sum())
            if int(row.n) != n_test or int(row.n_train) != len(source)-n_test:
                raise AssertionError("Fold counts do not agree with source block membership.")
    comparison = pd.read_csv(out/"Table3_reproduction_comparison.csv")
    diff_columns = [c for c in comparison if c.endswith("_difference")]
    within_rounding = bool((comparison[diff_columns].abs().to_numpy() <= .00500001).all())
    calibration = json.loads((RESULTS/"rh_calibration/reproduction_check.json").read_text())
    if not calibration["matches_stored_RH_within_1e-7_m"]:
        raise AssertionError("Stored RH calibration could not be reproduced.")
    result = {"technical_checks": "passed", "prediction_rows": len(pred), "blocks": 12,
              "checks": ["source keys and targets preserved", "one prediction per analysis row",
                         "finite layer predictions", "metrics independently recomputed",
                         "12-block score counts agree with input membership", "stored RH reproduced"],
              "matches_manuscript_Table3_within_reported_rounding": within_rounding,
              "matches_supplied_layer_RMSE_within_1e-3": bool(all(abs(x["difference"]) < .001 for x in source_comparison)),
              "comparison_with_supplied_layer_scores": source_comparison,
              "whole_manuscript_reproduced": False,
              "release_status": "draft; see docs/RELEASE_STATUS.md"}
    write_json(RESULTS/"verification.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
