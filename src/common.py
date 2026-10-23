"""Repository paths, scoring conventions, and figure export."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/data_pred_all_data_v2.csv"
RAW_METRICS = ROOT / "data/processed/GEDI_SIM_metrics.csv"
REGIONAL = ROOT / "data/regional/GEDI_2019_with_Predicted_AGBD.csv"
PERCENTILES = tuple(range(5, 101, 5))
RESULTS = ROOT / "results"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_cohort(path=DATA):
    df = pd.read_csv(path, dtype={"site": str, "ID": str})
    required = ["ID", "site", "ALS_AGB", "ALS_Slope", "ALS_DEM"]
    required += [f"{prefix}{p}" for prefix in ("Gs_RH", "ALS_RH", "Pred_RH")
                 for p in PERCENTILES]
    missing = sorted(set(required) - set(df))
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    if df[["site", "ID"]].isna().any().any() or df.duplicated(["site", "ID"]).any():
        raise ValueError("site + ID must be complete and unique.")
    values = df[[c for c in required if c not in ("site", "ID")]].to_numpy(float)
    if not np.isfinite(values).all():
        raise ValueError("Required numerical inputs contain missing/nonfinite values; no rows were dropped.")
    if df.site.nunique() != 12:
        raise ValueError("This reproduction expects the 12-block study dataset.")
    if not df.ALS_Slope.between(0, 90, inclusive="left").all():
        raise ValueError("ALS_Slope must be in degrees in [0, 90).")
    return df


def score(reference, prediction):
    y, p = np.asarray(reference, float), np.asarray(prediction, float)
    keep = np.isfinite(y) & np.isfinite(p)
    y, p = y[keep], p[keep]
    if not len(y):
        return {"n": 0, "r": None, "RMSE": None, "Bias": None, "R2": None}
    err = p-y
    sst = np.sum((y-y.mean())**2)
    return {"n": int(len(y)),
            "r": float(np.corrcoef(y, p)[0, 1]) if len(y) > 1 and np.std(y) > 0 and np.std(p) > 0 else None,
            "RMSE": float(np.sqrt(np.mean(err**2))),
            "Bias": float(np.mean(err)),
            "R2": float(1-np.sum(err**2)/sst) if sst > 0 else None}


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def save_figure(fig, folder, name):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg", "pdf"):
        fig.savefig(folder/f"{name}.{extension}", dpi=300, bbox_inches="tight")
