# Manuscript reproducibility code

This directory contains the minimal analysis and figure/table-generation code aligned to the revised manuscript (V25 analysis specification). It intentionally excludes exploratory one-off scripts, raw participant identifiers, raw sequencing files, intermediate audit files, and obsolete holdout-only analyses.

## Included code

| Script | Purpose | Main manuscript output |
|---|---|---|
| `01_wes_candidate_association.py` | Additive candidate-variant logistic regression in the WES aMCI cohort; covariate adjustment and five-test BH/Bonferroni correction; optional PC1–PC5 sensitivity | Table S10 / PC-adjusted sensitivity (Table S14) |
| `02_wgs_pooled_oof_models.py` | Five prespecified XGBoost configurations using one common stratified 5-fold pooled OOF partition; fold-local preprocessing; 2,000 paired class-stratified bootstrap CIs; strict-model outer-fold permutation importance | Table 3, Table 4, Table S11 |
| `03_wgs_candidate_logistic_sensitivity.py` | Candidate-level WGS logistic regression separated from ML attribution | Table S15 |
| `04_candidate_ld.py` | Pairwise complete dosage-based Pearson r and r² for the five candidate variants; LD table and heatmap | Table S6 / LD figure |
| `05_make_figures_tables.py` | Manuscript-ready model-performance and permutation-importance tables plus vector/raster figures from analysis outputs | Main/revised performance and attribution figures |

## Final V25 model specification

Candidate SNP dosage features:

- `rs7946`
- `rs25489`
- `rs28469095`
- `rs429358`
- `rs440446`

Global exclusions before model fitting:

- `ID`
- `FAMILY_ID`
- `PNTTM`
- `Ages` (redundant five-year age-band variable derived from `Test_Age`)

Strict leakage-reduced model (14 original features):

- `Test_Age`
- `YearsOfEducation`
- `APOE`
- `Age`
- `Job`
- `Education`
- `Sex`
- `Education_Level`
- `Handgrip`
- `rs7946`
- `rs25489`
- `rs28469095`
- `rs429358`
- `rs440446`

The APOE-sensitivity strict model removes `rs429358` and `APOE`, leaving 12 features. The Non-APOE SNP-only model contains `rs7946`, `rs25489`, and `rs28469095`.

## XGBoost settings

The pooled OOF script uses the fixed manuscript configuration:

- objective: `binary:logistic`
- eval_metric: `logloss`
- n_estimators: `500`
- max_depth: `3`
- learning_rate: `0.03`
- subsample: `0.8`
- colsample_bytree: `0.8`
- min_child_weight: `5`
- reg_lambda: `1.0`
- n_jobs: `-1`
- random_state: `42`
- 5 outer stratified folds, identical across all model configurations
- threshold: `0.50`
- 2,000 paired class-stratified bootstrap resamples for pooled OOF AUROC/AUPRC CIs
- strict-model permutation importance only in outer validation folds, 50 repeats per feature per fold

Sentinel values `66666`, `77777`, and `99999` are converted to missing values. Missing numeric values are not imputed; XGBoost handles them natively. Categorical variables are one-hot encoded using categories learned only from each training fold.

## Input expectations

### WES

A de-identified tabular file (`.tsv`, `.csv`, or `.xlsx`) containing:

- binary `aMCI_status` (0/1)
- `age`
- `sex`
- `APOE_e4`
- five additive SNP dosage columns (0/1/2)

Optional internal ancestry PCs can be supplied with `--pcs PC1 PC2 PC3 PC4 PC5`.

### WGS

A de-identified tabular file containing:

- binary outcome `Diagnose` or `Diagnosis` (0 control, 1 AD)
- five candidate SNP dosage columns
- clinical/epidemiologic variables used in the manuscript
- `APOE` genotype/proxy variable if used in the strict/full models

The script treats object/category/bool columns as categorical. If the source export stores categorical codes as integers, provide a one-column feature-name file with `--categorical-features`.

## Example commands

```bash
python 01_wes_candidate_association.py \
  --data WES_candidate_matrix.tsv \
  --outdir results/wes_primary

python 01_wes_candidate_association.py \
  --data WES_candidate_matrix_with_PCs.tsv \
  --pcs PC1 PC2 PC3 PC4 PC5 \
  --outdir results/wes_pc_adjusted

python 02_wgs_pooled_oof_models.py \
  --data WGS_epi_plus_5snps.tsv \
  --target Diagnose \
  --outdir results/wgs_oof

python 03_wgs_candidate_logistic_sensitivity.py \
  --data WGS_epi_plus_5snps.tsv \
  --target Diagnose \
  --age Age \
  --sex Sex \
  --outdir results/wgs_candidate_logistic

python 04_candidate_ld.py \
  --data WGS_epi_plus_5snps.tsv \
  --outdir results/ld

python 05_make_figures_tables.py \
  --model-performance results/wgs_oof/model_performance.tsv \
  --permutation results/wgs_oof/strict_permutation_importance.tsv \
  --outdir results/manuscript_assets
```

## Expected V25 sanity checks

When run on the final locked WGS analytic matrix (n=995), the model feature counts should be 5, 142, 14, 12, and 3 for the five-SNP, full clinical–genomic, strict leakage-reduced, strict without rs429358/APOE proxy, and Non-APOE SNP-only models, respectively. The exact numerical performance should be checked against the frozen manuscript analysis outputs before release.

## Interpretation boundary

These scripts reproduce exploratory association and internal predictive analyses. The pooled OOF estimates are internal cross-validated performance, not external or prospective validation. Permutation importance is descriptive model reliance and must not be interpreted as independent genetic association, biological mechanism, causality, or clinical utility.

## Data availability

Participant-level data are not bundled with this repository. Users must obtain the underlying datasets through the applicable cohort/data-access procedures and construct the de-identified analysis matrices described above.

## Versioning for publication

Before manuscript resubmission/publication:

1. create a GitHub release/tag (for example `v1.0.0`);
2. archive that exact release in Zenodo;
3. add the Zenodo DOI and software licence;
4. cite the frozen version in the manuscript Code Availability statement.
