#!/usr/bin/env python3
"""
Candidate-level WGS logistic regression sensitivity analysis.

Non-APOE candidates:
    AD ~ SNP dosage + Age + Sex + rs429358 dosage

APOE-region candidates (rs429358, rs440446):
    AD ~ SNP dosage + Age + Sex

This analysis is intentionally separate from XGBoost feature attribution.
Five-test Benjamini-Hochberg FDR and Bonferroni correction are reported.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


CANDIDATES = [
    ("rs7946", "PEMT", "Non-APOE candidate"),
    ("rs25489", "XRCC1", "Non-APOE candidate"),
    ("rs28469095", "NKPD1", "Non-APOE candidate"),
    ("rs429358", "APOE", "APOE-region candidate"),
    ("rs440446", "APOE-region", "APOE-region candidate"),
]


def read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(p)
    if p.suffix.lower() == ".csv":
        return pd.read_csv(p, low_memory=False)
    return pd.read_csv(p, sep="\t", low_memory=False)


def binary_encode(s: pd.Series, name: str) -> pd.Series:
    num = pd.to_numeric(s, errors="coerce")
    obs = sorted(pd.unique(num.dropna()))
    if set(obs).issubset({0, 1}):
        return num.astype(float)

    text = s.astype("string").str.strip().str.lower()
    mappings = {
        "control": 0.0, "cn": 0.0, "normal": 0.0,
        "ad": 1.0, "alzheimer": 1.0, "alzheimer's": 1.0, "case": 1.0,
        "female": 0.0, "f": 0.0, "male": 1.0, "m": 1.0,
    }
    mapped = text.map(mappings)
    if mapped.notna().sum() == text.notna().sum():
        return mapped
    if len(obs) == 2:
        return num.map({obs[0]: 0.0, obs[1]: 1.0})
    raise ValueError(f"{name} is not a recognizable binary field.")


def fit_logit(df: pd.DataFrame, outcome: str, snp: str, covariates: list[str]) -> dict:
    cols = [outcome, snp] + covariates
    d = df[cols].copy()
    d[outcome] = binary_encode(d[outcome], outcome)

    for c in [snp] + covariates:
        if c.lower() == "sex":
            try:
                d[c] = binary_encode(d[c], c)
                continue
            except ValueError:
                pass
        d[c] = pd.to_numeric(d[c], errors="coerce")

    d = d.dropna()
    if len(d) == 0:
        raise ValueError(f"No complete observations for {snp}.")
    if d[snp].nunique() < 2:
        raise ValueError(f"No genotype variation for {snp}.")

    X = sm.add_constant(d[[snp] + covariates], has_constant="add")
    y = d[outcome].astype(float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = sm.GLM(y, X, family=sm.families.Binomial()).fit()

    beta = float(fit.params[snp])
    se = float(fit.bse[snp])
    p = float(fit.pvalues[snp])
    z = 1.959963984540054

    return {
        "NMISS": int(len(d)),
        "Beta": beta,
        "SE": se,
        "OR": float(np.exp(beta)),
        "CI95_low": float(np.exp(beta - z * se)),
        "CI95_high": float(np.exp(beta + z * se)),
        "P": p,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--target", default="Diagnose")
    ap.add_argument("--age", default="Age")
    ap.add_argument("--sex", default="Sex")
    ap.add_argument("--apoe-dosage", default="rs429358", dest="apoe_dosage")
    args = ap.parse_args()

    df = read_table(args.data)
    target = args.target
    if target not in df.columns:
        aliases = [c for c in ["Diagnose", "Diagnosis"] if c in df.columns]
        if len(aliases) == 1:
            target = aliases[0]
        else:
            raise KeyError(f"Target {args.target!r} not found.")

    required = {target, args.age, args.sex, args.apoe_dosage}
    required.update(snp for snp, _, _ in CANDIDATES)
    missing = sorted(required.difference(df.columns))
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    rows = []
    for snp, gene, group in CANDIDATES:
        covars = [args.age, args.sex]
        if group == "Non-APOE candidate":
            covars.append(args.apoe_dosage)
        r = fit_logit(df, target, snp, covars)
        r.update({
            "Gene_or_region": gene,
            "SNP": snp,
            "Model_group": group,
            "Inheritance_model": "Additive dosage",
            "Covariates": ", ".join(covars),
        })
        rows.append(r)

    res = pd.DataFrame(rows)
    reject, q, _, _ = multipletests(res["P"], alpha=0.05, method="fdr_bh")
    res["FDR_BH_primary_5"] = q
    res["Bonferroni_primary_5"] = np.minimum(res["P"] * len(res), 1.0)
    res["FDR_reject_0.05"] = reject
    res["CI95"] = res.apply(
        lambda r: f'{r.CI95_low:.3f}–{r.CI95_high:.3f}', axis=1
    )

    def label(row):
        if row.SNP == "rs429358":
            return "Expected APOE-region positive-control signal."
        if row.SNP == "rs440446":
            return "APOE-region marker; not interpreted as APOE-independent."
        if row.FDR_BH_primary_5 < 0.05 or row.Bonferroni_primary_5 < 0.05:
            return "Corrected WGS candidate-level association; requires independent replication."
        if row.P < 0.05:
            return "Nominal WGS candidate-level association; not significant after correction."
        return "No corrected WGS candidate-level association evidence."

    res["Interpretation"] = res.apply(label, axis=1)
    res = res[
        ["Gene_or_region", "SNP", "Model_group", "Inheritance_model", "Covariates",
         "NMISS", "OR", "CI95", "Beta", "SE", "P",
         "FDR_BH_primary_5", "Bonferroni_primary_5", "FDR_reject_0.05",
         "Interpretation"]
    ]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    res.to_csv(outdir / "wgs_candidate_logistic_sensitivity.tsv", sep="\t", index=False)
    res.to_excel(outdir / "wgs_candidate_logistic_sensitivity.xlsx", index=False)

    print(res.to_string(index=False))
    print(f"\nWrote: {outdir / 'wgs_candidate_logistic_sensitivity.tsv'}")


if __name__ == "__main__":
    main()
