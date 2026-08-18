#!/usr/bin/env python3
"""
Minimal WES candidate-variant association analysis aligned to the revised manuscript.

Primary additive models
-----------------------
Non-APOE candidates:
    aMCI_status ~ SNP + age + sex + APOE_e4 [+ optional ancestry PCs]

APOE-region candidates:
    aMCI_status ~ SNP + age + sex [+ optional ancestry PCs]

Multiple-testing correction is applied across the five candidate-variant tests.
Input SNP columns must be additive dosages (0/1/2) for a consistent effect allele.
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
    suffix = p.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(p)
    if suffix == ".csv":
        return pd.read_csv(p, low_memory=False)
    return pd.read_csv(p, sep="\t", low_memory=False)


def encode_binary(series: pd.Series, name: str) -> pd.Series:
    """Accept numeric 0/1 or common binary string labels."""
    x = series.copy()
    numeric = pd.to_numeric(x, errors="coerce")
    observed = sorted(pd.unique(numeric.dropna()))
    if set(observed).issubset({0, 1}):
        return numeric.astype(float)

    text = x.astype("string").str.strip().str.lower()
    mapping_sets = [
        ({"female": 0, "f": 0, "male": 1, "m": 1}, "sex"),
        ({"control": 0, "cn": 0, "amci": 1, "case": 1}, "status"),
        ({"no": 0, "noncarrier": 0, "non-carrier": 0,
          "yes": 1, "carrier": 1}, "carrier"),
    ]
    for mapping, _ in mapping_sets:
        mapped = text.map(mapping)
        if mapped.notna().sum() == text.notna().sum():
            return mapped.astype(float)

    if len(observed) == 2:
        return numeric.map({observed[0]: 0.0, observed[1]: 1.0})

    raise ValueError(
        f"{name!r} must be binary (0/1) or a recognizable two-level field. "
        f"Observed non-missing values: {series.dropna().unique()[:10]}"
    )


def fit_one(
    df: pd.DataFrame,
    outcome: str,
    snp: str,
    covariates: list[str],
) -> dict:
    cols = [outcome, snp] + covariates
    d = df[cols].copy()

    d[outcome] = encode_binary(d[outcome], outcome)
    for c in [snp] + covariates:
        if c in {"sex", "Sex", "APOE_e4"}:
            try:
                d[c] = encode_binary(d[c], c)
                continue
            except ValueError:
                pass
        d[c] = pd.to_numeric(d[c], errors="coerce")

    d = d.dropna()
    if d.empty:
        raise ValueError(f"No complete observations for {snp} model.")
    if d[outcome].nunique() != 2:
        raise ValueError(f"Outcome is not binary after complete-case filtering for {snp}.")
    if d[snp].nunique() < 2:
        raise ValueError(f"{snp} has no usable dosage variation.")

    X = sm.add_constant(d[[snp] + covariates], has_constant="add")
    y = d[outcome].astype(float)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = sm.GLM(y, X, family=sm.families.Binomial())
        fit = model.fit()

    beta = float(fit.params[snp])
    se = float(fit.bse[snp])
    p = float(fit.pvalues[snp])
    lo = beta - 1.959963984540054 * se
    hi = beta + 1.959963984540054 * se

    return {
        "NMISS": int(len(d)),
        "Beta": beta,
        "SE": se,
        "OR": float(np.exp(beta)),
        "CI95_low": float(np.exp(lo)),
        "CI95_high": float(np.exp(hi)),
        "P": p,
    }


def interpretation(snp: str, p: float, fdr: float, bonf: float) -> str:
    if snp == "rs429358":
        return "Expected APOE-region positive-control signal; not a novel candidate locus."
    if snp == "rs440446":
        return "APOE-region marker; not interpreted as independent of APOE-related mechanisms."
    if fdr < 0.05 or bonf < 0.05:
        return "Candidate-panel corrected association; requires independent replication."
    if p < 0.05:
        return "Nominal exploratory candidate; not significant after candidate-panel correction."
    return "Attenuated/non-significant after covariate adjustment; exploratory candidate only."


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--outcome", default="aMCI_status")
    ap.add_argument("--age", default="age")
    ap.add_argument("--sex", default="sex")
    ap.add_argument("--apoe-e4", default="APOE_e4", dest="apoe_e4")
    ap.add_argument(
        "--pcs", nargs="*", default=[],
        help="Optional genotype-derived ancestry PCs, e.g. PC1 PC2 PC3 PC4 PC5."
    )
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    df = read_table(args.data)

    required = {args.outcome, args.age, args.sex}
    required.update(s for s, _, _ in CANDIDATES)
    if any(group == "Non-APOE candidate" for _, _, group in CANDIDATES):
        required.add(args.apoe_e4)
    required.update(args.pcs)
    missing = sorted(required.difference(df.columns))
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    rows = []
    for snp, gene, group in CANDIDATES:
        covariates = [args.age, args.sex] + list(args.pcs)
        if group == "Non-APOE candidate":
            covariates.insert(2, args.apoe_e4)

        r = fit_one(df, args.outcome, snp, covariates)
        r.update({
            "Gene_or_region": gene,
            "SNP": snp,
            "Model_group": group,
            "Inheritance_model": "Additive dosage",
            "Covariates": ", ".join(covariates),
        })
        rows.append(r)

    res = pd.DataFrame(rows)
    reject, q, _, _ = multipletests(res["P"].values, alpha=0.05, method="fdr_bh")
    res["FDR_BH_primary_5"] = q
    res["Bonferroni_primary_5"] = np.minimum(res["P"] * len(res), 1.0)
    res["FDR_reject_0.05"] = reject
    res["Interpretation"] = [
        interpretation(snp, p, fdr, bonf)
        for snp, p, fdr, bonf in zip(
            res["SNP"], res["P"], res["FDR_BH_primary_5"], res["Bonferroni_primary_5"]
        )
    ]
    res["CI95"] = res.apply(
        lambda x: f'{x["CI95_low"]:.3f}–{x["CI95_high"]:.3f}', axis=1
    )

    order = [
        "Gene_or_region", "SNP", "Model_group", "Inheritance_model", "Covariates",
        "NMISS", "OR", "CI95", "Beta", "SE", "P",
        "FDR_BH_primary_5", "Bonferroni_primary_5", "FDR_reject_0.05",
        "Interpretation",
    ]
    res = res[order]

    stem = "wes_candidate_association_pc_adjusted" if args.pcs else "wes_candidate_association"
    res.to_csv(outdir / f"{stem}.tsv", sep="\t", index=False)
    res.to_excel(outdir / f"{stem}.xlsx", index=False)

    print(res.to_string(index=False))
    print(f"\nWrote: {outdir / f'{stem}.tsv'}")


if __name__ == "__main__":
    main()
