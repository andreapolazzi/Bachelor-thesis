# OOF SHAP for tTF & bTF TabPFN regression (notebook 10)

Date: 2026-06-08  Status: approved (design)

## Goal

Bring the TabPFN **regression** interpretability in `analysis/primary/10_TabPFN.qmd` up to
the same standard as the `full2` classification SHAP flow. Currently the two regression
"Interpretability" sections (tTF = `tf_primary_rate`, bTF = `tf_blood_rate`) are
placeholders: they export the full feature matrix as `X_train.csv` and the first 50 rows
as `X_test.csv`, and point at `SHAP_tabpfn_local.py`, which only supports
`TabPFNClassifier`. Replace this with a proper **out-of-fold (OOF)** SHAP computation over a
single stratified 5-fold CV, using both explainers, in the natural target (log1p) space.

This generalizes the already-approved `2026-06-05-blood-inclusive-oof-shap` design from
classification to regression. Same pipeline shape: single stratified 5-fold OOF → both
explainers → row-aligned `shapviz` object in R with preserved metadata.

## Decisions (from brainstorming)

- **Explainers**: keep **both** — imputation (primary, self-verifiable) and
  recontextualization (cross-check / TabPFN-author-recommended).
- **Fold scheme**: single stratified 5-fold OOF, stratified on `alt_status`, `set.seed(42)`.
  One SHAP value per sample. (Not the 5×5 repeated CV used for the metrics sections.)
- **Compute**: local `TabPFNRegressor` on a GPU box (carmela, the dedicated no-scheduler
  TabPFN/SHAP box; nohup). The `tabpfn_client` API used by the CV-metrics sections is too
  slow / rate-limited for SHAP's hundreds of evals per sample. The launcher also works
  unchanged under PBS on chiron.

## 1. Python script — `--task {classification,regression}` flag

`scripts/SHAP_tabpfn_local.py`. New `--task` argument, default `classification`
(fully backward-compatible; the existing classification runs are unaffected).

When `--task regression`:

- Fit `TabPFNRegressor()` instead of `TabPFNClassifier()`.
- Drop all classifier-specific logic for this path: no `predict_proba`, no `proba_col`,
  no `class_index`. The per-row prediction function is simply `reg.predict`.
- **imputation** path: wrap `reg.predict` in the existing custom
  `shapiq.TabularExplainer(model=predict, data=X_bg, index="SV", max_order=1,
  imputer="marginal")`. Output space label = `"target"` (log1p of the rate).
  Self-verification, identical in spirit to today's classification check but against
  predictions:
  - per fold: `|baseline_fold − mean(reg.predict(X[train]))| ≤ BASELINE_TOL` (abort if any
    fold fails);
  - after assembly: global additivity `baseline_row + Σ shap_row ≈ reg.predict(row)`
    (report mean/max error; warn above `ADDITIVITY_WARN`).
- **recontextualization** path:
  `get_tabpfn_explainer(model=reg, data=X[train], labels=y[train], index="SV",
  max_order=1)` (no `class_index` — ignored for regressors). Output space label =
  `"target_recontext"`. Verify by additivity + rank-corr(reconstruction, OOF prediction),
  mirroring today's recontext check (must be strongly positive).

The OOF loop, original-order reassembly, and the full set of output files are unchanged:
`shap_values.csv`, `X_test_explained.csv`, `shap_baseline.txt`, `baselines.csv`,
`oof_pred.csv`, `fold_id.csv`, `meta.csv`, `feature_names.txt`, `shap_space.txt`.

Implementation note: branch on `args.task` inside `build_explainer` and the OOF loop;
factor the "predict function + expected baseline" per task so the verification code stays
shared. Single-split (non-folds) mode should honor `--task` too for symmetry, but OOF is
the path used here.

## 2. Notebook 10 — replace the two export chunks with CV exporters

Replace `ttf-export-shap` and `btf-export-shap` (currently `eval=FALSE`,
single-split dumps) with CV exporters mirroring the classification `export-shap-cv` chunk.

For each target:

- Rebuild from `reg_data_raw` using the **same filtering** the target already applies:
  - tTF: `drop_na()` then `filter(tf_primary_rate < 10)`;
  - bTF: `drop_na()` (no outlier filter).
- Apply the same `log1p` transforms as the existing `ttf-data` / `btf-data` chunks.
- **Keep metadata aside**, row-aligned: `donor_id, cancer_type, cancer_group, alt_status`
  (these are dropped from `reg_data_ttf/btf`, so capture them from the filtered raw frame
  before the `select(-...)`).
- Features = the existing `x_ttf_r` / `x_btf_r` (everything except the target and meta;
  tTF features include `tf_blood_rate`, bTF features include `tf_primary_rate`).
- One stratified 5-fold split: `set.seed(42)`, `caret::createFolds(factor(alt_status),
  k = 5)` → integer `folds` vector of length N in original row order.
- Write to `outputs/shap_input_ttf_cv/` and `outputs/shap_input_btf_cv/`:
  - `X_full.csv` — N × p features (column names = feature names);
  - `y_full.csv` — single column (log1p target), length N;
  - `folds.csv` — single column `fold` (1..5), length N, row-aligned;
  - `meta.csv` — `donor_id, cancer_type, cancer_group, alt_status`, length N, row-aligned.

## 3. Notebook 10 — replace the two load chunks

Replace `ttf-shap-load` / `btf-shap-load` to:

- Load OOF outputs from `outputs/shap_output_ttf_imp/`, `outputs/shap_output_ttf_rec/`,
  and the bTF equivalents, via the shared R helper.
- Build `shapviz` objects and plot **manually** (per user preference).
- Apply the rank01 recolor fix on beeswarm to counter right-skewed feature distributions,
  consistent with the classification section.
- Keep the "not yet available" guard + printed run command, updated to the regression
  launcher path.

## 4. R helper `scripts/reconstruct_shapviz.R`

Already provides `load_shap(dir)` and `load_meta(dir)`. Just add the regression object
lines (no new functions):

- `sv_ttf_imp <- load_shap(here("outputs","shap_output_ttf_imp"))`
- `sv_ttf_rec <- load_shap(here("outputs","shap_output_ttf_rec"))`
- `meta_ttf   <- load_meta(here("outputs","shap_output_ttf_imp"))`
- bTF trio analogously.

Output values are in target (log1p) space; no special handling needed beyond labeling.

## 5. Launcher `scripts/run_shap_reg.sh`

Modeled on `run_shap_full2.sh`. Runs four passes:

1. tTF imputation  → `outputs/shap_output_ttf_imp/`
2. tTF recontext   → `outputs/shap_output_ttf_rec/`
3. bTF imputation  → `outputs/shap_output_btf_imp/`
4. bTF recontext   → `outputs/shap_output_btf_rec/`

Each pass: `--task regression --folds <input>/folds.csv --input-dir <input> --output-dir
<output> --explainer <imputation|recontextualization> --n-jobs 4` with
`PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` (imputation OOM'd at 16 jobs in full2;
TabPFN is GPU-bound so CPU jobs give no speedup). Timestamped internal `tee` logging,
per-step exit code + summary. Works under nohup on carmela and under PBS on chiron.

## Out of scope

- No change to the CV-metrics sections — they keep the `tabpfn_client` 5×5 API flow.
  Only interpretability moves to local-GPU OOF.
- No repeated 5×5 SHAP; single 5-fold OOF, one SHAP value per sample.
- No automated plotting (user plots manually).
- No change to `run_tabpfn_shap.py`.
