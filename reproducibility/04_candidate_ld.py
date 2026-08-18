#!/usr/bin/env python3
"""
Pairwise linkage-disequilibrium summary for the five manuscript candidate variants.

LD is represented as squared Pearson correlation (r²) between additive genotype
dosages, using pairwise complete observations. Both r and r² are retained.
"""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


SNPS = ["rs429358", "rs440446", "rs28469095", "rs7946", "rs25489"]


def read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(p)
    if p.suffix.lower() == ".csv":
        return pd.read_csv(p, low_memory=False)
    return pd.read_csv(p, sep="\t", low_memory=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = read_table(args.data)
    missing = [s for s in SNPS if s not in df.columns]
    if missing:
        raise KeyError(f"Missing candidate SNP columns: {missing}")

    g = df[SNPS].apply(pd.to_numeric, errors="coerce")
    rows = []

    for a, b in combinations(SNPS, 2):
        d = g[[a, b]].dropna()
        if len(d) < 3 or d[a].nunique() < 2 or d[b].nunique() < 2:
            r = np.nan
        else:
            r = float(d[a].corr(d[b], method="pearson"))
        rows.append({
            "SNP1": a,
            "SNP2": b,
            "N_pairwise_complete": int(len(d)),
            "R": r,
            "R2": r * r if np.isfinite(r) else np.nan,
        })

    pairs = pd.DataFrame(rows)
    pairs.to_csv(outdir / "candidate_pairwise_LD.tsv", sep="\t", index=False)
    pairs.to_excel(outdir / "candidate_pairwise_LD.xlsx", index=False)

    r2 = pd.DataFrame(np.eye(len(SNPS)), index=SNPS, columns=SNPS)
    nmat = pd.DataFrame(np.nan, index=SNPS, columns=SNPS)
    for s in SNPS:
        nmat.loc[s, s] = int(g[s].notna().sum())
    for row in rows:
        r2.loc[row["SNP1"], row["SNP2"]] = row["R2"]
        r2.loc[row["SNP2"], row["SNP1"]] = row["R2"]
        nmat.loc[row["SNP1"], row["SNP2"]] = row["N_pairwise_complete"]
        nmat.loc[row["SNP2"], row["SNP1"]] = row["N_pairwise_complete"]

    r2.to_csv(outdir / "candidate_LD_r2_matrix.tsv", sep="\t")
    nmat.to_csv(outdir / "candidate_LD_n_matrix.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    im = ax.imshow(r2.to_numpy(dtype=float), vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(SNPS)), labels=SNPS, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(SNPS)), labels=SNPS)
    ax.set_title("Candidate-variant dosage LD (r²)")

    for i in range(len(SNPS)):
        for j in range(len(SNPS)):
            val = r2.iloc[i, j]
            label = "NA" if pd.isna(val) else f"{val:.3f}"
            ax.text(j, i, label, ha="center", va="center")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("r²")
    fig.tight_layout()
    for ext in ["pdf", "svg", "png"]:
        fig.savefig(outdir / f"candidate_LD_heatmap.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(pairs.to_string(index=False))
    print(f"\nOutputs written under: {outdir}")


if __name__ == "__main__":
    main()
