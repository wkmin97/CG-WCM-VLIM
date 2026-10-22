#!/usr/bin/env python3
"""Summarize supplied regional columns directly, with no undocumented rescaling.

This is a diagnostic of the supplied CSV, not a verified recreation of
manuscript Figure 14 or a new application of the correction model.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.common import REGIONAL, RESULTS, save_figure, write_json


def main():
    d = pd.read_csv(REGIONAL)
    required = ["Slope1", "rh95", "Cal_RH95", "GLAD_Ht1", "ETH_Ht1",
                "GEDI_AGBD", "Cal_AGBD", "CCI_AGBD1", "CFATD_AGBD"]
    absent = sorted(set(required)-set(d))
    if absent:
        raise ValueError(f"Missing regional columns: {absent}")
    if not d.Slope1.between(0, 90, inclusive="left").all():
        raise ValueError("Slope1 must contain degrees in [0,90).")
    labels = ["0–5", "5–10", "10–15", "15–20", "20–25", "25–30", "30–35", "35–40", "40–45", "≥45"]
    d["slope_bin"] = pd.cut(d.Slope1, [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 90],
                            labels=labels, right=False)
    panels = [(["rh95", "Cal_RH95", "GLAD_Ht1", "ETH_Ht1"], "Height (m)", "Stored height columns"),
              (["GEDI_AGBD", "Cal_AGBD", "CFATD_AGBD", "CCI_AGBD1"], "Stored AGBD (units to confirm)", "Stored biomass columns")]
    colors = ["#D87970", "#6EAB94", "#A0A0A0", "#7D9FC2"]
    summaries = []
    fig, axes = plt.subplots(2, 1, figsize=(12, 8.2), constrained_layout=True)
    for ax, (columns, ytitle, title) in zip(axes, panels):
        for j, (column, color) in enumerate(zip(columns, colors)):
            groups, positions = [], []
            for k, label in enumerate(labels):
                v = d.loc[d.slope_bin.eq(label), column].to_numpy(float)
                v = v[np.isfinite(v)]
                if not len(v):
                    continue
                groups.append(v)
                positions.append(k+(j-1.5)*.19)
                summaries.append({"column": column, "slope_bin": label, "n": len(v),
                                  "minimum": float(v.min()), "q25": float(np.quantile(v,.25)),
                                  "median": float(np.median(v)), "q75": float(np.quantile(v,.75)),
                                  "maximum": float(v.max())})
            boxes = ax.boxplot(groups, positions=positions, widths=.16, patch_artist=True,
                               showfliers=False, manage_ticks=False,
                               medianprops={"color": "#333333", "linewidth": .8})
            for patch in boxes["boxes"]:
                patch.set_facecolor(color)
                patch.set_linewidth(.5)
            ax.plot([], [], color=color, linewidth=6, label=column)
        ax.set_xticks(range(len(labels)), labels)
        ax.set(xlabel="Terrain slope (degrees)", ylabel=ytitle, title=title)
        ax.legend(ncol=4, frameon=False, loc="upper left", fontsize=8)
        ax.grid(axis="y", alpha=.2)
    fig.suptitle("Regional CSV diagnostics — values plotted as stored", fontsize=12)
    save_figure(fig, RESULTS/"figures", "regional_stored_columns_diagnostic")
    plt.close(fig)
    out = RESULTS/"regional_diagnostics"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summaries).to_csv(out/"slope_bin_summaries.csv", index=False)
    write_json(out/"provenance_review.json", {
        "rows": len(d), "extra_filtering": "none", "value_rescaling": "none",
        "degrade_flag_nonzero_count": int(d.degrade_fl.ne(0).sum()) if "degrade_fl" in d else None,
        "official_L4A_provenance_confirmed": False,
        "regional_footprint_region_and_year_confirmed": False,
        "Cal_AGBD_rescaling_in_original_R": "13apply model.R contains CAL_AGB = Cal_AGBD*3+50; no justification was provided; not applied here",
        "reproduces_manuscript_Figure14": False,
        "review": "Resolve column units, product provenance, geographic/year metadata, correction application, and quality filters before release."})
    print(f"Summarized {len(d):,} regional rows without rescaling; Figure 14 provenance requires review.")


if __name__ == "__main__":
    main()
