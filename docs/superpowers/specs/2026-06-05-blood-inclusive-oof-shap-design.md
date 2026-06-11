# Blood-inclusive, 5-fold OOF SHAP for TabPFN classification

Date: 2026-06-05 Status: approved (design)

## Goal

Repeat the TabPFN ALT-status SHAP analysis with three changes:

1.  Include `tf_blood_rate` as a feature (switch to the dataset that carries it).
2.  Compute SHAP **out-of-fold** over a single stratified **5-fold CV** (each sample explained by a model that did not train on it), mirroring the TreeSHAP OOF pattern in `analysis/primary/02_ALT_status_full.Rmd`, for stability.
3.  Keep **both** explainers (imputation + recontextualization), writing to new output directories `output_imputation_full2/` and `output_recontext_full2/`.

The aggregated result must reconstruct into a `shapviz` object in R, with row-aligned metadata (`cancer_type`, `cancer_group`, `alt_status`, `donor_id`) preserved so plots can be colored/faceted by label later.

## 1. Data & features (R, notebook 10)

-   Source: `data/processed/PCAWG_primary.xlsx` (replaces `PCAWG_variables_PCA_woK.xlsx` for the classification flow). This file carries `tf_blood_rate` and the metadata columns the analysis needs.
-   Features (14): `tf_primary_rate`, `tf_blood_rate`, `telomere_insertion_rate`, `telomere_content_log2`, `TERT_FPKM`, and the 9 `*_singleton_dist` motifs (`ATAGGG, CTAGGG, GTAGGG, TAAGGG, TCAGGG, TGAGGG, TTCGGG, TTGGGG, TTTGGG`).
-   Label: `alt_status` encoded 1 = ALT-high, 0 = ALT-low.
-   Metadata carried (not features): `donor_id`, `cancer_type`, `cancer_group`.
-   Row filtering: drop rows with any missing value among the 14 features (notably missing `TERT_FPKM`), consistent with the PCA-with-TERT flow. Let `N` = remaining rows.
-   Side effect: the existing R-side CV-performance section (which reads `model_data`) re-computes its metrics on this new, blood-inclusive data. The cached `outputs/tabpfn_cv_results.rds` must be regenerated (delete/ignore cache to re-run).

## 2. Fold generation & export (R, notebook 10)

-   One stratified 5-fold split on `alt_status` using `caret::createFolds(..., k = 5)` with a fixed seed (`set.seed(42)`), exactly as notebook 02 does.
-   Build an integer fold-id vector `folds` of length `N` (values 1..5), in the original row order of the filtered data.
-   Export to `outputs/shap_input_cv/`:
    -   `X_full.csv` — N x 14 feature matrix (column names = feature names).
    -   `y_full.csv` — single column `y` (0/1), length N.
    -   `folds.csv` — single column `fold` (1..5), length N, row-aligned to `X_full`.
    -   `meta.csv` — columns `donor_id, cancer_type, cancer_group, alt_status`, length N, row-aligned to `X_full`.

## 3. Python OOF mode (`scripts/SHAP_tabpfn_local.py`)

Add an **optional** `--folds PATH` argument. Behavior:

-   When `--folds` is NOT given: current single-split behavior is unchanged (reads `X_train.csv`, `X_test.csv`, `y_train.csv`). Backward compatible.
-   When `--folds` IS given (OOF mode):
    -   Read `X_full.csv` + `y_full.csv` from `--input-dir`, and the fold vector from `--folds`. Optionally read `meta.csv` if present (to copy into output).
    -   Loop folds `i = 1..K`:
        -   `train = rows where fold != i`, `test = rows where fold == i`.
        -   Fit `TabPFNClassifier()` on `(X[train], y[train])`.
        -   `proba_col = clf.classes_.index(class_index)`.
        -   Build the chosen explainer using `X[train]` as background:
            -   imputation: `shapiq.TabularExplainer(model=proba_predict, data=X[train],   index="SV", max_order=1, imputer="marginal")`, where `proba_predict(Z) = clf.predict_proba(Z)[:, proba_col]`.
            -   recontextualization: `get_tabpfn_explainer(model=clf, data=X[train],   labels=y[train], index="SV", max_order=1, class_index=class_index)`.
        -   `explain_X(X[test], budget=256, n_jobs=4, random_state=seed)`.
        -   Collect, for each test row: SHAP row, `baseline_value`, the OOF prediction `proba_predict(X[test])`, the fold id, and the row's original index.
    -   Reassemble all collected rows back into **original sample order** (by original index).

### Per-fold verification (OOF mode)

-   imputation: for each fold, assert `|baseline_fold - mean(proba_predict(X[train]))| <=   BASELINE_TOL` (abort if any fold fails — wrong class guard). After assembly, report global additivity: `baseline_row + sum(shap_row)` vs OOF prediction (mean/max error).
-   recontextualization: after assembly, report rank correlation between `baseline_row + sum(shap_row)` and the OOF prediction (must be strongly positive; warn if negative — sign/class flip).

### Outputs (per explainer dir: `output_imputation_full2/`, `output_recontext_full2/`)

Filenames kept identical to the single-split outputs so the R helper works unchanged: - `shap_values.csv` — N x 14, original order. - `X_test_explained.csv` — N x 14 feature values, original order (name kept for R reuse). - `shap_baseline.txt` — scalar = mean of per-row baselines (for `shapviz`). - `feature_names.txt`, `shap_space.txt`.

Additional OOF files (additive, ignored by the current R helper): - `baselines.csv` — per-row baseline (length N), for exact per-row reconstruction. - `oof_pred.csv` — per-row OOF predicted P(class=1). - `fold_id.csv` — per-row fold id (1..5). - `meta.csv` — copied through from input, row-aligned.

## 4. R reconstruction with labels (`scripts/reconstruct_shapviz.R`)

Extend the existing helper (do not create a new file): - `load_shap(dir)` stays as-is (reads the 4 core files, builds `shapviz`). - Add optional `meta` loading: if `meta.csv` exists in `dir`, read it and return it alongside, or expose a small `load_meta(dir)` helper. Provide: - `sv_imp2 <- load_shap(here("outputs","shap_output","output_imputation_full2"))` - `sv_rec2 <- load_shap(here("outputs","shap_output","output_recontext_full2"))` - `meta2  <- load_meta(here("outputs","shap_output","output_imputation_full2"))` - The helper only builds objects + meta; plotting stays manual (per user preference). - Keep the existing `sv_imp` / `sv_rec` (full, non-CV) lines intact.

## Runtime & launcher

Budget 256, both explainers, OOF over all N samples. **Run on the HPC (chiron) GPU queue with `--n-jobs 16`** (not carmela).

-   Launcher: `scripts/run_shap_full2.sh` — a PBS job that activates the micromamba env (`source /ngs/software/micromamba/env.sh; micromamba activate   /ngs/software/conda/envs/tabpfn2`), runs both explainers (imputation then recontext) on `outputs/shap_input_cv/`, writing to `output_imputation_full2/` and `output_recontext_full2/`.
-   Logging: do the timestamped internal `tee` (`exec > >(tee -a logs/...) 2>&1`), NOT `#PBS -o ...${PBS_JOBID}...` (that variable does not expand in directive lines and the logs are lost). Each step prints its own exit code + a summary.
-   Submit: `qsubmit.pl -g 1 -n 16 -s /ngs/iflores/andrea/run_shap_full2.sh`.
-   The conda env activation must come BEFORE python (an inline qsubmit without activation fails with `ModuleNotFoundError: No module named 'tabpfn'`).

## Out of scope

-   No change to regression flows or to `run_tabpfn_shap.py`.
-   No repeated (5x5) CV; a single stratified 5-fold OOF pass, matching notebook 02.
-   No automated plotting (user plots manually for flexibility).
