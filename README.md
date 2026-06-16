# Telomere biology, ALT status and telomere‑fusion rate in cancer

Bachelor thesis project. A reproducible analysis pipeline that uses **telomere‑derived
genomic features** to (1) classify **ALT status** (Alternative Lengthening of Telomeres)
and (2) model **telomere‑fusion (TF) rate** in primary and metastatic tumours, with an
emphasis on *honest* model comparison under heavy class imbalance and *interpretable*
predictions via SHAP.

The work spans two cohorts, classical statistical/ML models in R, and a foundation
model for tabular data (**TabPFN**) interpreted with out‑of‑fold SHAP.

---

## Research questions

1. **Can telomere features predict ALT status?** ALT‑high tumours are rare
   (~6 % of complete cases, 49 / 817 in the primary cohort), so the central
   challenge is class imbalance and avoiding metrics that flatter a majority‑class
   classifier.
2. **What drives the telomere‑fusion rate?** Modelled both for tumour TF rate
   (`tf_primary_rate`, "tTF") and blood TF rate (`tf_blood_rate`, "bTF") as a
   zero‑inflated, skewed regression target.
3. **Which features actually carry signal, and which are artifacts?** Several
   notebooks specifically probe confounding — `TERT_FPKM` missingness, cancer‑type
   (mesenchymal vs non‑mesenchymal) structure, and gene‑level covariates.

## Cohorts and features

| Cohort | Source | Shape | Used for |
|--------|--------|-------|----------|
| **Primary** | PCAWG | wide, one row per sample (`data/processed/PCAWG_primary.xlsx`) | ALT classification + TF‑rate regression |
| **Metastatic** | Hartwig | long‑format, dedup by `patient_id` before per‑sample analysis | EDA, TVR comparison primary vs metastatic |

The core predictor set (`features14`) is:

- Telomere‑fusion rates: `tf_primary_rate`, `tf_blood_rate`
- `telomere_insertion_rate`, `telomere_content_log2`
- `TERT_FPKM` (TERT expression)
- Nine **telomere variant repeat (TVR) singleton‑distance** motifs:
  `ATAGGG`, `CTAGGG`, `GTAGGG`, `TAAGGG`, `TCAGGG`, `TGAGGG`, `TTCGGG`, `TTGGGG`, `TTTGGG`

## Methods

- **Shared cross‑validation.** Every modelling notebook uses one definition of the
  folds — `make_folds()` (`scripts/make_folds.R`), a **5‑fold × 5‑repeat stratified**
  scheme with `set.seed(21)`. SHAP uses repeat 1. Because the partitions are
  byte‑identical across models, the comparisons are *paired*.
- **Classification models.** Logistic regression, random forest (with `treeshap`),
  LASSO (`glmnet`), and TabPFN — all scored with the same metric code
  (`scripts/clf_metrics.R`): AUC, **PR‑AUC / average precision** (the honest summary
  under imbalance), sensitivity, specificity, precision, F1, balanced accuracy.
- **Regression models.** Linear, LASSO, random forest, and a **two‑part model**
  (presence logistic + intensity Gamma GLM with Duan smearing) on a `log1p`‑transformed
  target. Metrics via `scripts/reg_metrics.R` (RMSE, MAE, Spearman).
- **Fair comparisons.** Notebooks `20_classification_comparison` and
  `21_regression_comparison` put every model on identical footing (same cohort, same
  features, same folds, same metrics) — the single source of truth for the headline
  numbers.
- **Interpretability.** TabPFN predictions are explained with **out‑of‑fold SHAP**
  (each sample is explained by a model that was *not* trained on it). Two explainer
  paradigms are supported (see `scripts/SHAP_tabpfn_local.py`):
  - *imputation* — `shapiq.TabularExplainer` in **probability space**; the script
    self‑verifies additivity and baseline and aborts if the class/scale is wrong.
  - *recontextualization* — `tabpfn_extensions` remove‑and‑recontextualize in **raw
    (logit) space**; validate by additivity + rank correlation, not baseline magnitude.

## Repository layout

```
.
├── data/
│   ├── raw/                  # source spreadsheets (PCAWG, Hartwig, drivers, coverage…)
│   └── processed/            # cleaned, joined tables consumed by the notebooks
├── analysis/
│   ├── primary/              # PCAWG primary-tumour notebooks (.Rmd / .qmd)
│   ├── metastatic/           # Hartwig metastatic notebooks
│   └── .venv-tabpfn/         # Python venv for TabPFN + SHAP (reticulate)
├── scripts/                  # shared R helpers, Python TabPFN/SHAP, HPC run scripts
├── outputs/                  # figures, comparison tables, SHAP inputs/outputs, OOF CSVs
├── docs/superpowers/         # design specs and plans for the OOF‑SHAP work
└── graphify-out/             # knowledge-graph snapshot of the codebase
```

> **Note on what is tracked.** `data/`, `outputs/`, rendered `*.html`, and the
> `*_files` / `*_cache` render artifacts are git‑ignored (see `.gitignore`). A fresh
> clone therefore contains **source notebooks and scripts only** — the input
> spreadsheets are not redistributed and must be supplied separately.

## Notebooks

### Primary (`analysis/primary/`)

| Notebook | Topic |
|----------|-------|
| `01_ALT_status_TVR_only`, `02_ALT_status_full` | ALT classification from TVR motifs / full feature set (RF + treeshap) |
| `06_feature_ablation` | feature ablation for ALT prediction |
| `08_cell_lines` | cell-line PCA / biplots |
| `09_LASSO` | LASSO classification |
| `10_TabPFN` | **TabPFN classification + regression with out‑of‑fold SHAP** (core notebook) |
| `13_*`, `15_*` | tTF / bTF regression (primary, and with the other compartment's TF added) |
| `16_*_noTERT`, `17_*_noTERT` | regressions excluding `TERT_FPKM` |
| `14_TERT_missingness` | TERT_FPKM missingness as a confounder/artifact |
| `18_ALT_status_genes`, `22_tTF_genes` | gene-level covariates (`lme4` mixed models) |
| `19_ALT_status_blood_TF` | blood-TF‑inclusive ALT classification |
| `20_classification_comparison`, `21_regression_comparison` | **fair head‑to‑head model comparisons** |
| `mesenchymal_groups`, `features_distribution`, `TF-blood_vs_primary` | confounding, distributions, blood vs tumour TF |

### Metastatic (`analysis/metastatic/`)

`01_EDA` → `02_duplicate_columns` → `03_further_EDA` → `4_metastatic_dataset` build and
explore the Hartwig dataset; `6_TVRs_comparison` compares TVR singleton distances
between primary and metastatic tumours.

## Key scripts (`scripts/`)

- **Data prep:** `process_raw_data.R`, `clean_columns.R` (`clean_pcawg_data`,
  `clean_drivers_data`) — `data/raw/*` → `data/processed/*`.
- **Modelling helpers:** `make_folds.R` (shared folds), `clf_metrics.R`,
  `reg_metrics.R`, `useful_functions.R`, `useful_plots.R`,
  `check_tumor_group_discordance.R`.
- **CV / SHAP I/O:** `export_cv_inputs.R` (writes `X_full.csv` / `y_full.csv` /
  `folds_5x5.csv` for the GPU runs), `reconstruct_shapviz.R` (rebuild `shapviz`
  objects in R from the Python SHAP outputs).
- **TabPFN + SHAP (Python):**
  - `SHAP_tabpfn_local.py` — local TabPFN fit + SHAP (classification & regression,
    both explainer paradigms); writes `shap_values.csv`, `shap_baseline.txt`,
    `shap_space.txt`, etc.
  - `tabpfn_cv_oof.py` — out‑of‑fold TabPFN predictions for the fair comparisons.
  - `run_tabpfn_shap.py` — legacy hosted‑client (single split) SHAP workflow.
  - `diag_recontext.py`, `test_oof_smoke.py` — diagnostics / smoke tests.
- **HPC / GPU runners (`run_*.sh`):** wrap the Python scripts for the `chiron` HPC GPU
  queue or the schedulerless `carmela` GPU box (local `tabpfn`, no API token;
  `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True`, `--n-jobs 4`).

## Reproducing the analysis

### 1. R environment

Open `R.Rproj` (RStudio) and render the notebooks (`rmarkdown::render` / Quarto).
Notebooks use `here::here()` for paths, so run from the project root. Main R packages:

`tidyverse`, `here`, `readxl`, `writexl`, `caret`, `randomForest`, `glmnet`, `pROC`,
`treeshap`, `shapviz`, `lme4`, `broom`, `patchwork`, `ggpubr`, `rstatix`, `naniar`,
`corrplot`, `umap`, `reticulate`.

Render the data‑prep step first:

```bash
Rscript scripts/process_raw_data.R   # data/raw/* -> data/processed/*
```

### 2. Python / TabPFN environment

`10_TabPFN.qmd` calls Python through `reticulate` against `analysis/.venv-tabpfn`
(Python 3.9). Key pins:

```
tabpfn==6.4.1   tabpfn-extensions (PriorLabs)   shapiq==1.2.1   shap==0.49.1
torch==2.8.0    scikit-learn==1.6.1   pandas==2.3.3   numpy==2.0.2
```

The hosted‑client path (`run_tabpfn_shap.py`) needs a `TABPFN_TOKEN` env var; the
local/HPC path does not.

### 3. TabPFN CV + SHAP on GPU

```bash
# 1) export model inputs locally
Rscript scripts/export_cv_inputs.R

# 2) sync outputs/shap_input_* and cmp_reg_* dirs to the GPU box, then:
bash scripts/run_tabpfn_cv.sh        # OOF predictions for notebooks 20 & 21
bash scripts/run_shap_reg.sh         # regression OOF SHAP   (and run_shap_*.sh variants)

# 3) sync the *_oof.csv / shap_output_* back into outputs/, then re-render
#    notebooks 10, 20, 21 for the final tables and figures.
```

## Outputs

`outputs/` collects the rendered artifacts: PCA biplots, the
`classification_comparison` / `regression_comparison` tables and figures, TabPFN OOF
prediction CSVs (`tabpfn_cv_oof.csv`, `tabpfn_reg_{ttf,btf}_oof.csv`), and the SHAP
input/output bundles (`shap_input_*`, `shap_output_*`).

## Caveats worth knowing

- **`TERT_FPKM` is a deceptive predictor.** Its apparent importance is partly a
  shared‑variance + missingness artifact; the `*_noTERT` and `14_TERT_missingness`
  notebooks exist to quantify this.
- **Cancer type confounds the signal.** Mesenchymal vs non‑mesenchymal origin tracks
  with ALT/feature signal — see `mesenchymal_groups`.
- **Use PR‑AUC, not accuracy/AUC alone**, given ~6 % positives.
- The metastatic dataset is **long‑format**: deduplicate by `patient_id`
  (23316 → 3934 rows) before any per‑sample analysis.

---

*Author: Andrea Polazzi. Code: R (Quarto/R Markdown) + Python (TabPFN, shapiq).*
