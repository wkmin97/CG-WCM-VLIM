#!/usr/bin/env python3
"""Check supplied inputs and summarize stored predictions without refitting."""
import json
import platform
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scipy
import sklearn

from src.common import DATA, RAW_METRICS, ROOT, RESULTS, load_cohort, score, sha256, write_json


def main():
    out = RESULTS/"data_audit"
    out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((ROOT/"data/manifest.json").read_text())
    for item in manifest:
        path = ROOT/item["path"]
        if not path.is_file() or sha256(path) != item["sha256"]:
            raise ValueError(f"Input missing or checksum mismatch: {item['path']}")
    d = load_cohort()
    mapping = pd.read_csv(ROOT/"config/block_mapping.csv", dtype={"site": str})
    if set(mapping.site) != set(d.site) or mapping.site.duplicated().any():
        raise ValueError("Block mapping does not uniquely cover the analysis cohort.")
    raw = pd.read_csv(RAW_METRICS, dtype={"site": str, "ID": str})
    keys = ["site", "ID"]
    missing = d[keys].merge(raw[keys].drop_duplicates(), on=keys, how="left", indicator=True)
    missing = missing[missing._merge.eq("left_only")].drop(columns="_merge")
    missing.to_csv(out/"cohort_keys_missing_from_broader_table.csv", index=False)
    summary = d.groupby("site").agg(n=("ID", "size"), DEM_min=("ALS_DEM", "min"),
                    DEM_max=("ALS_DEM", "max"), slope_min=("ALS_Slope", "min"),
                    slope_max=("ALS_Slope", "max"), AGBD_max=("ALS_AGB", "max"))
    summary.reset_index().merge(mapping, on="site").sort_values("block").to_csv(out/"block_summary.csv", index=False)
    columns = [c for c in d if c.startswith("Pred_") and not c.startswith("Pred_RH")]
    columns += [c for c in ["SIM_Gs_AGBD"] if c in d]
    pd.DataFrame([{"stored_column": c, "evaluation": "stored_predictions_only; CV_provenance_not_established",
                   **score(d.ALS_AGB, d[c])} for c in columns]).to_csv(out/"stored_prediction_metrics.csv", index=False)
    report = {"analysis_rows": len(d), "analysis_columns": len(d.columns), "blocks": d.site.nunique(),
              "broader_metrics_rows": len(raw), "cohort_keys_absent_from_broader_table": len(missing),
              "AGBD_above_100_count": int((d.ALS_AGB > 100).sum()),
              "data_checksums_verified": len(manifest), "input_sha256": sha256(DATA),
              "row_filtering": "none", "status": "processed inputs checked; whole-study release gaps remain"}
    write_json(out/"data_audit.json", report)
    versions = {"python": sys.version, "platform": platform.platform(),
                **{m.__name__: m.__version__ for m in [np, pd, scipy, sklearn, matplotlib]}}
    write_json(RESULTS/"environment.json", versions)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
