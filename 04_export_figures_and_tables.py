#!/usr/bin/env python3
"""Export conditional LOBO block scores, residual plots, and a manuscript comparison."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.common import ROOT, RESULTS, score, save_figure


def main():
    folder = RESULTS/"agbd_conditional_LOBO"
    path = folder/"LOBO_AGBD_all_predictions.csv"
    if not path.is_file():
        raise FileNotFoundError("Run 03_validate_agbd_lobo.py first.")
    d = pd.read_csv(path, dtype={"site": str, "ID": str})
    blocks = pd.read_csv(ROOT/"config/block_mapping.csv", dtype={"site": str})
    d = d.merge(blocks[["site", "block"]], on="site", validate="many_to_one")
    variants = [("raw", "LOBO_layer_Gs_RH", "#D98A9B"),
                ("corrected", "LOBO_layer_Pred_RH", "#699FC6")]
    rows = []
    for block, sub in d.groupby("block", sort=True):
        row = {"block": f"Block {block}", "site": sub.site.iloc[0], "n": len(sub)}
        for label, col, _ in variants:
            row.update({f"{label}_{key}": val for key, val in score(sub.ALS_AGB, sub[col]).items() if key != "n"})
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(folder/"Table3_conditional_LOBO_recomputed.csv", index=False)
    reported = pd.read_csv(ROOT/"data/reported/Table3_reported.csv")
    comparison = table.merge(reported, on="block", suffixes=("_recomputed", "_manuscript"), validate="one_to_one")
    for label, _, _ in variants:
        for metric in ["r", "RMSE", "Bias"]:
            name = f"{label}_{metric}"
            comparison[f"{name}_difference"] = comparison[f"{name}_recomputed"]-comparison[f"{name}_manuscript"]
    comparison.to_csv(folder/"Table3_reproduction_comparison.csv", index=False)
    plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans", "axes.spines.top": False,
                         "axes.spines.right": False, "svg.fonttype": "none"})
    fig, axes = plt.subplots(3, 4, figsize=(12.8, 8.1), sharey=True, constrained_layout=True)
    for ax, (block, sub) in zip(axes.flat, d.groupby("block", sort=True)):
        arrays = [sub[col]-sub.ALS_AGB for _, col, _ in variants]
        boxes = ax.boxplot(arrays, positions=[1, 2], widths=.55, patch_artist=True,
                           showfliers=False, medianprops={"color": "#303030"})
        for box, (_, _, color) in zip(boxes["boxes"], variants):
            box.set_facecolor(color)
        ax.axhline(0, color="#666666", linestyle="--", linewidth=.8)
        ax.set_xticks([1, 2], ["Original RH", "Stored corrected RH"], fontsize=8)
        ax.set_title(f"Block {block} (n={len(sub):,})")
        ax.grid(axis="y", alpha=.15)
    for ax in axes[:, 0]:
        ax.set_ylabel("Prediction − reference (Mg ha$^{-1}$)")
    fig.suptitle("Conditional AGBD-model LOBO residuals across UAV blocks", fontsize=12)
    save_figure(fig, RESULTS/"figures", "LOBO_AGBD_residuals_by_block")
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.4), constrained_layout=True)
    for ax, (label, col, color) in zip(axes, variants):
        met = score(d.ALS_AGB, d[col])
        lo = min(0., float(d[col].min()), float(d.ALS_AGB.min()))
        hi = max(float(d[col].max()), float(d.ALS_AGB.max()))*1.04
        ax.scatter(d.ALS_AGB, d[col], s=3, alpha=.2, color=color, rasterized=True)
        ax.plot([lo, hi], [lo, hi], "--", color="#555555", linewidth=.8)
        ax.set(xlim=(lo, hi), ylim=(lo, hi), xlabel="UAV reference AGBD (Mg ha$^{-1}$)",
               ylabel="LOBO prediction (Mg ha$^{-1}$)", title="Original RH" if label == "raw" else "Stored corrected RH")
        ax.set_aspect("equal")
        annotation = f"r = {met['r']:.3f}\nRMSE = {met['RMSE']:.2f} Mg ha$^{{-1}}$\nBias = {met['Bias']:.2f} Mg ha$^{{-1}}$\nn = {met['n']:,}"
        ax.text(.04, .96, annotation, transform=ax.transAxes, va="top", fontsize=9)
    fig.suptitle("Conditional AGBD-model LOBO validation", fontsize=12)
    save_figure(fig, RESULTS/"figures", "LOBO_AGBD_scatter")
    plt.close(fig)
    print("Exported conditional LOBO figures and Table 3 comparison; original manuscript values remain separate.")


if __name__ == "__main__":
    main()
