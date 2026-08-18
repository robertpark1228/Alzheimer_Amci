#!/usr/bin/env python3
"""
Create compact manuscript tables and vector/raster figures from the locked analysis
outputs. This script does not refit models; it only formats analysis results.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


MODEL_ORDER = [
    "SNP-only",
    "Full clinical-genomic",
    "Strict leakage-reduced clinical-genomic",
    "Strict leakage-reduced without rs429358/APOE proxy",
    "Non-APOE SNP-only",
]


def read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def save_figure(fig, stem: Path) -> None:
    for ext in ["pdf", "svg", "png"]:
        fig.savefig(stem.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def performance_table(df: pd.DataFrame) -> pd.DataFrame:
    needed = [
        "Model", "N", "N_features", "AUROC", "AUROC_CI_low", "AUROC_CI_high",
        "AUPRC", "AUPRC_CI_low", "AUPRC_CI_high", "Accuracy",
        "Balanced_accuracy", "Sensitivity", "Specificity", "Brier",
        "TN", "FP", "FN", "TP",
    ]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise KeyError(f"model performance file missing columns: {missing}")

    x = df[needed].copy()
    x["AUROC (95% CI)"] = x.apply(
        lambda r: f'{r["AUROC"]:.3f} ({r["AUROC_CI_low"]:.3f}–{r["AUROC_CI_high"]:.3f})',
        axis=1,
    )
    x["AUPRC (95% CI)"] = x.apply(
        lambda r: f'{r["AUPRC"]:.3f} ({r["AUPRC_CI_low"]:.3f}–{r["AUPRC_CI_high"]:.3f})',
        axis=1,
    )
    x["Confusion matrix (TN, FP, FN, TP)"] = x.apply(
        lambda r: f'{int(r["TN"])}, {int(r["FP"])}, {int(r["FN"])}, {int(r["TP"])}',
        axis=1,
    )
    out = x[
        ["Model", "N", "N_features", "AUROC (95% CI)", "AUPRC (95% CI)",
         "Accuracy", "Balanced_accuracy", "Sensitivity", "Specificity", "Brier",
         "Confusion matrix (TN, FP, FN, TP)"]
    ].rename(columns={
        "N_features": "Features",
        "Balanced_accuracy": "Balanced accuracy",
    })
    return out


def plot_performance(df: pd.DataFrame, outdir: Path) -> None:
    order_map = {m: i for i, m in enumerate(MODEL_ORDER)}
    x = df.copy()
    x["_order"] = x["Model"].map(order_map)
    x = x.sort_values("_order").reset_index(drop=True)

    y = np.arange(len(x))
    auroc = x["AUROC"].to_numpy(float)
    lo = x["AUROC_CI_low"].to_numpy(float)
    hi = x["AUROC_CI_high"].to_numpy(float)
    err = np.vstack([auroc - lo, hi - auroc])

    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    ax.errorbar(auroc, y, xerr=err, fmt="o", capsize=4)
    ax.set_yticks(y, labels=x["Model"])
    ax.invert_yaxis()
    ax.set_xlim(0.45, 1.01)
    ax.set_xlabel("Pooled OOF AUROC (95% bootstrap CI)")
    ax.set_title("Model discrimination across prespecified configurations")
    ax.axvline(0.5, linewidth=1, linestyle="--")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    save_figure(fig, outdir / "Figure_model_performance_ladder")

    focus_names = [
        "Full clinical-genomic",
        "Strict leakage-reduced clinical-genomic",
        "Strict leakage-reduced without rs429358/APOE proxy",
        "Non-APOE SNP-only",
    ]
    f = x[x["Model"].isin(focus_names)].copy()
    f["_order2"] = f["Model"].map({m: i for i, m in enumerate(focus_names)})
    f = f.sort_values("_order2")
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    xx = np.arange(len(f))
    ax.plot(xx, f["AUROC"].to_numpy(float), marker="o")
    ax.set_xticks(xx, labels=[
        "Full clinical–\ngenomic",
        "Strict leakage-\nreduced",
        "Strict without\nrs429358/APOE",
        "Non-APOE\nSNP-only",
    ])
    ax.set_ylim(0.45, 1.02)
    ax.set_ylabel("Pooled OOF AUROC")
    ax.set_title("Reduction in discrimination after removing diagnostic-proximal variables")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    save_figure(fig, outdir / "Figure_performance_reduction")


def permutation_table(df: pd.DataFrame) -> pd.DataFrame:
    mean_col = "Permutation_importance_mean_AUROC_decrease"
    sd_col = "Permutation_importance_SD"
    needed = ["Feature", mean_col, sd_col, "N_permutations"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise KeyError(f"permutation file missing columns: {missing}")
    x = df.copy()
    if "Rank" not in x.columns:
        x = x.sort_values(mean_col, ascending=False).reset_index(drop=True)
        x.insert(0, "Rank", np.arange(1, len(x) + 1))
    return x[["Rank", "Feature", mean_col, sd_col, "N_permutations"]]


def plot_permutation(df: pd.DataFrame, outdir: Path) -> None:
    mean_col = "Permutation_importance_mean_AUROC_decrease"
    sd_col = "Permutation_importance_SD"
    x = df.sort_values(mean_col, ascending=True).copy()

    fig, ax = plt.subplots(figsize=(8.4, 6.6))
    y = np.arange(len(x))
    ax.barh(
        y,
        x[mean_col].to_numpy(float),
        xerr=x[sd_col].to_numpy(float),
        capsize=2,
    )
    ax.set_yticks(y, labels=x["Feature"])
    ax.axvline(0, linewidth=1)
    ax.set_xlabel("Mean validation-fold AUROC decrease after permutation")
    ax.set_title("Strict leakage-reduced model: fold-wise permutation importance")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    save_figure(fig, outdir / "Figure_strict_permutation_importance")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-performance", required=True)
    ap.add_argument("--permutation", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    perf = read_tsv(args.model_performance)
    perm = read_tsv(args.permutation)

    t3 = performance_table(perf)
    t3.to_csv(outdir / "Table3_model_performance.tsv", sep="\t", index=False)
    t3.to_excel(outdir / "Table3_model_performance.xlsx", index=False)

    t4 = permutation_table(perm)
    t4.to_csv(outdir / "Table4_strict_permutation_importance.tsv", sep="\t", index=False)
    t4.to_excel(outdir / "Table4_strict_permutation_importance.xlsx", index=False)

    plot_performance(perf, outdir)
    plot_permutation(t4, outdir)

    print(f"Wrote tables and figures under: {outdir}")


if __name__ == "__main__":
    main()
