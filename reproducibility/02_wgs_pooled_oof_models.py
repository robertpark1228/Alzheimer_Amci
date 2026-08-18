#!/usr/bin/env python3
"""
Five prespecified WGS clinical-genomic models using one common stratified 5-fold
pooled out-of-fold (OOF) evaluation.

Aligned to the revised manuscript V25 specification:
- n=995 final WGS analytic cohort
- global exclusions: ID, FAMILY_ID, PNTTM, Ages
- fixed XGBoost hyperparameters
- fold-local one-hot encoding
- no numeric imputation; XGBoost native missing-value handling
- probability threshold 0.50
- 2,000 paired class-stratified bootstrap resamples of pooled OOF predictions
- strict-model permutation importance only on outer validation folds
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from xgboost import XGBClassifier


SEED = 42
N_SPLITS = 5
N_BOOT = 2000
N_PERM = 50
THRESHOLD = 0.50

CANDIDATE_SNPS = ["rs7946", "rs25489", "rs28469095", "rs429358", "rs440446"]
NON_APOE_SNPS = ["rs7946", "rs25489", "rs28469095"]
GLOBAL_EXCLUSIONS = ["ID", "FAMILY_ID", "PNTTM", "Ages"]

STRICT_FEATURES_V25 = [
    "Test_Age",
    "YearsOfEducation",
    "APOE",
    "Age",
    "Job",
    "Education",
    "Sex",
    "Education_Level",
    "Handgrip",
    "rs7946",
    "rs25489",
    "rs28469095",
    "rs429358",
    "rs440446",
]

MODEL_ORDER = [
    "SNP-only",
    "Full clinical-genomic",
    "Strict leakage-reduced clinical-genomic",
    "Strict leakage-reduced without rs429358/APOE proxy",
    "Non-APOE SNP-only",
]


def read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(p)
    if p.suffix.lower() == ".csv":
        return pd.read_csv(p, low_memory=False)
    return pd.read_csv(p, sep="\t", low_memory=False)


def read_feature_list(path: str | None) -> list[str]:
    if not path:
        return []
    vals = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            vals.append(line)
    return vals


def resolve_target(df: pd.DataFrame, requested: str) -> str:
    if requested in df.columns:
        return requested
    aliases = ["Diagnose", "Diagnosis"]
    found = [x for x in aliases if x in df.columns]
    if len(found) == 1:
        return found[0]
    raise KeyError(f"Target {requested!r} not found. Available aliases found: {found}")


def binary_target(s: pd.Series) -> pd.Series:
    num = pd.to_numeric(s, errors="coerce")
    observed = sorted(pd.unique(num.dropna()))
    if set(observed).issubset({0, 1}):
        return num.astype(float)

    text = s.astype("string").str.strip().str.lower()
    mapped = text.map({
        "control": 0.0, "cn": 0.0, "normal": 0.0,
        "ad": 1.0, "alzheimer": 1.0, "alzheimer's": 1.0,
        "case": 1.0,
    })
    if mapped.notna().sum() == text.notna().sum():
        return mapped

    if len(observed) == 2:
        return num.map({observed[0]: 0.0, observed[1]: 1.0})
    raise ValueError(f"Outcome must be binary. Observed: {s.dropna().unique()[:10]}")


def replace_sentinels(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    sentinels = [66666, 77777, 99999, "66666", "77777", "99999"]
    return out.replace(sentinels, np.nan)


def determine_categoricals(
    df: pd.DataFrame,
    features: list[str],
    forced: set[str],
) -> list[str]:
    cats = []
    for c in features:
        if c in forced:
            cats.append(c)
            continue
        dtype = df[c].dtype
        if (
            pd.api.types.is_object_dtype(dtype)
            or isinstance(dtype, pd.CategoricalDtype)
            or pd.api.types.is_bool_dtype(dtype)
            or pd.api.types.is_string_dtype(dtype)
        ):
            cats.append(c)
    return cats


def prepare_X(
    df: pd.DataFrame,
    features: list[str],
    forced_categoricals: set[str],
) -> tuple[pd.DataFrame, list[str], list[str]]:
    missing = [c for c in features if c not in df.columns]
    if missing:
        raise KeyError(f"Missing model features: {missing}")

    X = df[features].copy()
    categorical = determine_categoricals(X, features, forced_categoricals)
    numeric = [c for c in features if c not in categorical]

    for c in numeric:
        X[c] = pd.to_numeric(X[c], errors="coerce")

    for c in categorical:
        mask = X[c].isna()
        X[c] = X[c].astype(str)
        X.loc[mask, c] = np.nan

    return X, numeric, categorical


def build_pipeline(numeric: list[str], categorical: list[str]) -> Pipeline:
    transformers = []
    if numeric:
        transformers.append(("num", "passthrough", numeric))
    if categorical:
        transformers.append((
            "cat",
            OneHotEncoder(handle_unknown="ignore", sparse_output=False),
            categorical,
        ))
    prep = ColumnTransformer(
        transformers=transformers,
        remainder="drop",
        verbose_feature_names_out=False,
    )
    model = XGBClassifier(
        objective="binary:logistic",
        eval_metric="logloss",
        n_estimators=500,
        max_depth=3,
        learning_rate=0.03,
        subsample=0.8,
        colsample_bytree=0.8,
        min_child_weight=5,
        reg_lambda=1.0,
        n_jobs=-1,
        random_state=SEED,
        missing=np.nan,
        verbosity=0,
    )
    return Pipeline([("preprocess", prep), ("model", model)])


def make_bootstrap_indices(y: np.ndarray, n_boot: int, seed: int) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    neg = np.flatnonzero(y == 0)
    pos = np.flatnonzero(y == 1)
    indices = []
    for _ in range(n_boot):
        b0 = rng.choice(neg, size=len(neg), replace=True)
        b1 = rng.choice(pos, size=len(pos), replace=True)
        idx = np.concatenate([b0, b1])
        rng.shuffle(idx)
        indices.append(idx)
    return indices


def percentile_ci(values: Iterable[float]) -> tuple[float, float]:
    x = np.asarray(list(values), dtype=float)
    return float(np.nanpercentile(x, 2.5)), float(np.nanpercentile(x, 97.5))


def summarize_predictions(
    y: np.ndarray,
    p: np.ndarray,
    bootstrap_indices: list[np.ndarray],
) -> dict:
    pred = (p >= THRESHOLD).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()

    auroc_boot, auprc_boot = [], []
    for idx in bootstrap_indices:
        yy, pp = y[idx], p[idx]
        auroc_boot.append(roc_auc_score(yy, pp))
        auprc_boot.append(average_precision_score(yy, pp))

    auroc = roc_auc_score(y, p)
    auprc = average_precision_score(y, p)
    auroc_lo, auroc_hi = percentile_ci(auroc_boot)
    auprc_lo, auprc_hi = percentile_ci(auprc_boot)

    sensitivity = recall_score(y, pred, pos_label=1, zero_division=0)
    specificity = tn / (tn + fp) if (tn + fp) else np.nan
    ppv = precision_score(y, pred, pos_label=1, zero_division=0)
    npv = tn / (tn + fn) if (tn + fn) else np.nan

    return {
        "N": len(y),
        "AUROC": auroc,
        "AUROC_CI_low": auroc_lo,
        "AUROC_CI_high": auroc_hi,
        "AUPRC": auprc,
        "AUPRC_CI_low": auprc_lo,
        "AUPRC_CI_high": auprc_hi,
        "Accuracy": accuracy_score(y, pred),
        "Balanced_accuracy": balanced_accuracy_score(y, pred),
        "Sensitivity": sensitivity,
        "Specificity": specificity,
        "PPV": ppv,
        "NPV": npv,
        "Brier": brier_score_loss(y, p),
        "TN": int(tn), "FP": int(fp), "FN": int(fn), "TP": int(tp),
    }


def permute_outer_validation(
    pipe: Pipeline,
    X_val: pd.DataFrame,
    y_val: np.ndarray,
    fold: int,
    n_repeats: int,
) -> list[dict]:
    base_prob = pipe.predict_proba(X_val)[:, 1]
    base_auc = roc_auc_score(y_val, base_prob)
    rows = []

    for feature_idx, feature in enumerate(X_val.columns):
        rng = np.random.default_rng(SEED + fold * 10000 + feature_idx * 100)
        values = X_val[feature].to_numpy(copy=True)
        for repeat in range(n_repeats):
            Xp = X_val.copy()
            perm = rng.permutation(len(Xp))
            Xp[feature] = values[perm]
            auc_perm = roc_auc_score(y_val, pipe.predict_proba(Xp)[:, 1])
            rows.append({
                "Fold": fold,
                "Repeat": repeat + 1,
                "Feature": feature,
                "Base_AUROC": base_auc,
                "Permuted_AUROC": auc_perm,
                "AUROC_decrease": base_auc - auc_perm,
            })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--target", default="Diagnose")
    ap.add_argument("--outdir", required=True)
    ap.add_argument(
        "--categorical-features",
        help="Optional one-column text file forcing named integer-coded fields to categorical.",
    )
    ap.add_argument(
        "--strict-features",
        help="Optional one-column text file overriding the locked V25 14-feature strict list.",
    )
    ap.add_argument("--bootstrap", type=int, default=N_BOOT)
    ap.add_argument("--permutation-repeats", type=int, default=N_PERM)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = replace_sentinels(read_table(args.data))
    target = resolve_target(df, args.target)
    y_series = binary_target(df[target])

    if y_series.isna().any():
        keep = y_series.notna()
        df = df.loc[keep].reset_index(drop=True)
        y_series = y_series.loc[keep].reset_index(drop=True)
    y = y_series.astype(int).to_numpy()

    for c in CANDIDATE_SNPS:
        if c not in df.columns:
            raise KeyError(f"Missing candidate SNP column: {c}")

    strict_features = read_feature_list(args.strict_features) or STRICT_FEATURES_V25
    forced_cats = set(read_feature_list(args.categorical_features))

    forbidden = set(GLOBAL_EXCLUSIONS + [target])
    full_features = [c for c in df.columns if c not in forbidden]

    model_features = {
        "SNP-only": CANDIDATE_SNPS,
        "Full clinical-genomic": full_features,
        "Strict leakage-reduced clinical-genomic": strict_features,
        "Strict leakage-reduced without rs429358/APOE proxy": [
            c for c in strict_features if c not in {"rs429358", "APOE"}
        ],
        "Non-APOE SNP-only": NON_APOE_SNPS,
    }

    for model_name, features in model_features.items():
        missing = [c for c in features if c not in df.columns]
        if missing:
            raise KeyError(f"{model_name}: missing columns {missing}")

    skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
    splits = list(skf.split(np.zeros(len(y)), y))
    folds = np.full(len(y), -1, dtype=int)
    for fold, (_, va) in enumerate(splits, start=1):
        folds[va] = fold

    bootstrap_indices = make_bootstrap_indices(y, args.bootstrap, SEED + 2026)

    performance_rows = []
    pred_frame = pd.DataFrame({
        "row_index": np.arange(len(y)),
        "Outcome": y,
        "Fold": folds,
    })
    all_perm_rows = []

    feature_type_audit = {}

    for model_name in MODEL_ORDER:
        features = model_features[model_name]
        X, numeric, categorical = prepare_X(df, features, forced_cats)
        feature_type_audit[model_name] = {
            "n_original_features": len(features),
            "features": features,
            "numeric": numeric,
            "categorical": categorical,
        }

        oof = np.full(len(y), np.nan, dtype=float)
        for fold, (tr, va) in enumerate(splits, start=1):
            pipe = build_pipeline(numeric, categorical)
            pipe.fit(X.iloc[tr], y[tr])
            oof[va] = pipe.predict_proba(X.iloc[va])[:, 1]

            if model_name == "Strict leakage-reduced clinical-genomic":
                all_perm_rows.extend(
                    permute_outer_validation(
                        pipe,
                        X.iloc[va].copy(),
                        y[va],
                        fold=fold,
                        n_repeats=args.permutation_repeats,
                    )
                )

        if np.isnan(oof).any():
            raise RuntimeError(f"Missing OOF predictions for {model_name}")

        metrics = summarize_predictions(y, oof, bootstrap_indices)
        metrics["Model"] = model_name
        metrics["N_features"] = len(features)
        performance_rows.append(metrics)
        pred_frame[model_name] = oof

    perf = pd.DataFrame(performance_rows)
    perf["Model"] = pd.Categorical(perf["Model"], categories=MODEL_ORDER, ordered=True)
    perf = perf.sort_values("Model").reset_index(drop=True)
    perf["Model"] = perf["Model"].astype(str)
    perf["AUROC_95CI"] = perf.apply(
        lambda r: f'{r.AUROC:.3f} ({r.AUROC_CI_low:.3f}–{r.AUROC_CI_high:.3f})', axis=1
    )
    perf["AUPRC_95CI"] = perf.apply(
        lambda r: f'{r.AUPRC:.3f} ({r.AUPRC_CI_low:.3f}–{r.AUPRC_CI_high:.3f})', axis=1
    )
    perf["Confusion_matrix_TN_FP_FN_TP"] = perf.apply(
        lambda r: f'{int(r.TN)}, {int(r.FP)}, {int(r.FN)}, {int(r.TP)}', axis=1
    )

    perf.to_csv(outdir / "model_performance.tsv", sep="\t", index=False)
    perf.to_excel(outdir / "model_performance.xlsx", index=False)
    pred_frame.to_csv(outdir / "oof_predictions.tsv", sep="\t", index=False)

    perm_raw = pd.DataFrame(all_perm_rows)
    perm_raw.to_csv(outdir / "strict_permutation_importance_raw.tsv", sep="\t", index=False)
    perm_summary = (
        perm_raw.groupby("Feature", as_index=False)["AUROC_decrease"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={
            "mean": "Permutation_importance_mean_AUROC_decrease",
            "std": "Permutation_importance_SD",
            "count": "N_permutations",
        })
        .sort_values("Permutation_importance_mean_AUROC_decrease", ascending=False)
        .reset_index(drop=True)
    )
    perm_summary.insert(0, "Rank", np.arange(1, len(perm_summary) + 1))
    perm_summary.to_csv(
        outdir / "strict_permutation_importance.tsv", sep="\t", index=False
    )
    perm_summary.to_excel(
        outdir / "strict_permutation_importance.xlsx", index=False
    )

    audit = {
        "seed": SEED,
        "n_splits": N_SPLITS,
        "bootstrap_resamples": args.bootstrap,
        "permutation_repeats_per_feature_per_fold": args.permutation_repeats,
        "threshold": THRESHOLD,
        "global_exclusions": GLOBAL_EXCLUSIONS,
        "target": target,
        "class_counts": {
            "0": int((y == 0).sum()),
            "1": int((y == 1).sum()),
        },
        "model_features": feature_type_audit,
        "xgboost": {
            "objective": "binary:logistic",
            "eval_metric": "logloss",
            "n_estimators": 500,
            "max_depth": 3,
            "learning_rate": 0.03,
            "subsample": 0.8,
            "colsample_bytree": 0.8,
            "min_child_weight": 5,
            "reg_lambda": 1.0,
            "n_jobs": -1,
            "random_state": SEED,
            "gamma": "XGBoost default (not explicitly set)",
            "reg_alpha": "XGBoost default (not explicitly set)",
            "tree_method": "XGBoost default (not explicitly set)",
        },
    }
    (outdir / "analysis_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print(perf[
        ["Model", "N", "N_features", "AUROC_95CI", "AUPRC_95CI",
         "Accuracy", "Balanced_accuracy", "Sensitivity", "Specificity", "Brier",
         "Confusion_matrix_TN_FP_FN_TP"]
    ].to_string(index=False))
    print(f"\nOutputs written under: {outdir}")


if __name__ == "__main__":
    main()
