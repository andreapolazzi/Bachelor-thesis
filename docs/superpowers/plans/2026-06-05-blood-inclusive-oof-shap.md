# Blood-inclusive 5-fold OOF SHAP Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Repeat the TabPFN ALT-status SHAP analysis with `tf_blood_rate` included, computed out-of-fold over a stratified 5-fold CV, for both imputation and recontextualization explainers, with label-preserving outputs reconstructable in R.

**Architecture:** R (notebook 10) switches to `PCAWG_primary.xlsx`, builds the 14-feature matrix + stratified folds + metadata and exports them to `outputs/shap_input_cv/`. The Python script `SHAP_tabpfn_local.py` gains an optional `--folds` OOF mode that trains TabPFN per fold, explains the held-out fold, and concatenates SHAP in original order. A PBS launcher runs both explainers on the HPC GPU queue. The R helper reconstructs `shapviz` objects plus a row-aligned metadata frame.

**Tech Stack:** R (dplyr, caret, readr, shapviz), Python (tabpfn, shapiq, tabpfn_extensions, numpy, pandas), PBS/micromamba on chiron HPC.

**Spec:** `docs/superpowers/specs/2026-06-05-blood-inclusive-oof-shap-design.md`

---

## File structure

- Modify: `scripts/SHAP_tabpfn_local.py` — add `--folds` OOF mode (new functions + main branch).
- Create: `scripts/test_oof_smoke.py` — fast synthetic smoke test for OOF mode (runs in the `tabpfn2` env on the server).
- Modify: `analysis/primary/10_TabPFN.qmd` — data-source swap, 14-feature set, fold + metadata export chunk.
- Create: `scripts/run_shap_full2.sh` — PBS launcher, both explainers, internal `tee` logging.
- Modify: `scripts/reconstruct_shapviz.R` — add `load_meta()` and `full2` object lines.

> **Note on test environment:** `tabpfn`/`shapiq` are NOT installed on the local mac. The smoke test (Task 1) must be run on a machine with the `tabpfn2` env (carmela or an HPC interactive/GPU node). All `Run:` commands in Task 1 are server commands run after `micromamba activate /ngs/software/conda/envs/tabpfn2`.

---

## Task 1: Python OOF mode in `SHAP_tabpfn_local.py`

**Files:**
- Create: `scripts/test_oof_smoke.py`
- Modify: `scripts/SHAP_tabpfn_local.py`

- [ ] **Step 1: Write the failing smoke test**

Create `scripts/test_oof_smoke.py`:

```python
#!/usr/bin/env python3
"""Fast synthetic smoke test for the OOF (--folds) mode of SHAP_tabpfn_local.py.
Run in the tabpfn2 env (server/GPU box):
    python scripts/test_oof_smoke.py
Exits 0 on success, non-zero on failure.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "SHAP_tabpfn_local.py"


def make_inputs(d: Path, n=40, n_feat=4, k=5, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, n_feat))
    # signal: positive class when feature 0 is high
    logits = 2.5 * X[:, 0] - 1.0
    y = (rng.uniform(size=n) < 1 / (1 + np.exp(-logits))).astype(int)
    # guarantee both classes present
    y[0], y[1] = 0, 1
    feats = [f"f{j}" for j in range(n_feat)]
    pd.DataFrame(X, columns=feats).to_csv(d / "X_full.csv", index=False)
    pd.DataFrame({"y": y}).to_csv(d / "y_full.csv", index=False)
    folds = (np.arange(n) % k) + 1
    pd.DataFrame({"fold": folds}).to_csv(d / "folds.csv", index=False)
    pd.DataFrame({"donor_id": [f"D{i}" for i in range(n)],
                  "cancer_type": rng.choice(["A", "B"], n),
                  "cancer_group": rng.choice(["M", "N"], n),
                  "alt_status": np.where(y == 1, "ALT-high", "ALT-low")}
                 ).to_csv(d / "meta.csv", index=False)
    return feats, n


def check_outputs(out: Path, feats, n):
    for f in ["shap_values.csv", "X_test_explained.csv", "shap_baseline.txt",
              "feature_names.txt", "shap_space.txt", "baselines.csv",
              "oof_pred.csv", "fold_id.csv", "meta.csv"]:
        assert (out / f).exists(), f"missing output {f}"
    shap = pd.read_csv(out / "shap_values.csv")
    assert shap.shape == (n, len(feats)), f"shap shape {shap.shape}"
    assert not shap.isna().any().any(), "NaN in shap values (a fold did not fill)"
    base = pd.read_csv(out / "baselines.csv")["baseline"].to_numpy()
    pred = pd.read_csv(out / "oof_pred.csv")["oof_pred"].to_numpy()
    recon = base + shap.to_numpy().sum(axis=1)
    return recon, pred


def run(mode, inp, out, extra):
    cmd = [sys.executable, str(SCRIPT), "--explainer", mode,
           "--input-dir", str(inp), "--output-dir", str(out),
           "--folds", str(inp / "folds.csv"), "--budget", "16",
           "--n-jobs", "1", "--class-index", "1", "--seed", "0"] + extra
    print("RUN:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        feats, n = make_inputs(d)

        out_imp = d / "out_imp"
        run("imputation", d, out_imp, [])
        recon, pred = check_outputs(out_imp, feats, n)
        add_err = np.abs(recon - pred)
        print(f"imputation additivity mean/max: {add_err.mean():.3g}/{add_err.max():.3g}")
        assert add_err.max() < 0.05, "imputation additivity too large"
        assert pd.read_csv(out_imp / "shap_space.txt", header=None).iloc[0, 0] == "probability"

        out_rec = d / "out_rec"
        run("recontextualization", d, out_rec, [])
        recon_r, pred_r = check_outputs(out_rec, feats, n)
        ar = np.argsort(np.argsort(pred_r)); rr = np.argsort(np.argsort(recon_r))
        rho = np.corrcoef(ar, rr)[0, 1]
        print(f"recontext rank corr(recon, pred): {rho:.3f}")
        assert rho > 0.3, "recontext reconstruction does not track prediction"

    print("SMOKE TEST PASSED")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the smoke test to verify it fails**

Run (server, tabpfn2 env): `python scripts/test_oof_smoke.py`
Expected: FAIL — `SHAP_tabpfn_local.py` does not yet accept `--folds` (argparse error: `unrecognized arguments: --folds`).

- [ ] **Step 3: Add the `--folds` argument**

In `scripts/SHAP_tabpfn_local.py`, in `parse_args()`, add after the `--explainer` argument block:

```python
    p.add_argument("--folds", default=None,
                   help="Path to folds.csv (single column of 1..K fold ids, row-aligned "
                        "to X_full.csv). When given, runs OOF mode and reads X_full.csv + "
                        "y_full.csv (and meta.csv if present) from --input-dir instead of "
                        "the X_train/X_test/y_train single-split files.")
```

- [ ] **Step 4: Add the shared explainer/extraction helpers**

In `scripts/SHAP_tabpfn_local.py`, add these module-level functions just above `def main():`:

```python
def build_explainer(explainer_kind, clf, X_bg, y_bg, proba_col, class_index, budget):
    """Return (explainer, space) for the chosen paradigm, with X_bg as background."""
    if explainer_kind == "imputation":
        def proba_predict(Z):
            return clf.predict_proba(Z)[:, proba_col]
        explainer = shapiq.TabularExplainer(
            model=proba_predict, data=X_bg, index="SV", max_order=1, imputer="marginal",
        )
        return explainer, "probability"
    explainer = get_tabpfn_explainer(
        model=clf, data=X_bg, labels=y_bg, index="SV", max_order=1,
        class_index=class_index,
    )
    return explainer, "logit"


def extract_shap(iv_list, n_features):
    """Return (shap_matrix [n x n_features], baselines [n]) from a shapiq iv list."""
    shap_matrix = np.array(
        [[float(iv[(j,)]) for j in range(n_features)] for iv in iv_list]
    )
    baselines = np.array([float(iv.baseline_value) for iv in iv_list])
    return shap_matrix, baselines


def run_oof(args, input_dir, output_dir):
    """Out-of-fold SHAP: train per fold, explain held-out fold, assemble in order."""
    X_df = pd.read_csv(input_dir / "X_full.csv")
    y = pd.read_csv(input_dir / "y_full.csv").iloc[:, 0].to_numpy().astype(int)
    folds = pd.read_csv(args.folds).iloc[:, 0].to_numpy().astype(int)
    feature_names = list(X_df.columns)
    n_features = len(feature_names)
    X = X_df.to_numpy().astype(float)
    N = X.shape[0]
    if not (len(y) == N and len(folds) == N):
        print(f"ERROR: row mismatch X={N} y={len(y)} folds={len(folds)}", file=sys.stderr)
        sys.exit(1)

    class_index = args.class_index
    meta_path = input_dir / "meta.csv"
    meta_df = pd.read_csv(meta_path) if meta_path.exists() else None

    shap_all = np.full((N, n_features), np.nan)
    base_all = np.full(N, np.nan)
    pred_all = np.full(N, np.nan)
    space = None

    print(f"OOF mode: N={N}, features={n_features}, folds={sorted(np.unique(folds))}")
    print(f"Explainer: {args.explainer}  budget={args.budget}  class_index={class_index}")

    for i in sorted(np.unique(folds)):
        tr = folds != i
        te = folds == i
        print(f"\n--- fold {i}: train={tr.sum()} explain={te.sum()} ---", flush=True)
        clf = TabPFNClassifier()
        clf.fit(X[tr], y[tr])
        classes = list(clf.classes_)
        if class_index not in classes:
            print(f"ERROR fold {i}: class {class_index} not in {classes}", file=sys.stderr)
            sys.exit(1)
        proba_col = classes.index(class_index)

        def proba_predict(Z, _clf=clf, _c=proba_col):
            return _clf.predict_proba(Z)[:, _c]

        if args.explainer == "imputation":
            expected = float(proba_predict(X[tr]).mean())

        explainer, space = build_explainer(
            args.explainer, clf, X[tr], y[tr], proba_col, class_index, args.budget
        )
        iv_list = explainer.explain_X(
            X[te], budget=args.budget, n_jobs=args.n_jobs,
            random_state=args.seed, verbose=True,
        )
        shap_fold, base_fold = extract_shap(iv_list, n_features)
        pred_fold = proba_predict(X[te])

        if args.explainer == "imputation":
            bdiff = abs(float(base_fold.mean()) - expected)
            print(f"  fold {i} baseline {base_fold.mean():.4f} vs expected {expected:.4f} "
                  f"(|diff|={bdiff:.4f})")
            if bdiff > BASELINE_TOL:
                print(f"ERROR fold {i}: baseline mismatch -> wrong class. Aborting.",
                      file=sys.stderr)
                sys.exit(2)

        idx = np.where(te)[0]
        shap_all[idx] = shap_fold
        base_all[idx] = base_fold
        pred_all[idx] = pred_fold

    if np.isnan(shap_all).any():
        print("ERROR: some rows were never explained (NaN remains).", file=sys.stderr)
        sys.exit(1)

    recon = base_all + shap_all.sum(axis=1)
    print(f"\nVerification ({space} space):")
    if args.explainer == "imputation":
        add_err = np.abs(recon - pred_all)
        print(f"  additivity error (mean/max): {add_err.mean():.4g} / {add_err.max():.4g}")
    else:
        ar = np.argsort(np.argsort(pred_all))
        rr = np.argsort(np.argsort(recon))
        rho = float(np.corrcoef(ar, rr)[0, 1])
        print(f"  rank corr(reconstruction, OOF P(class={class_index})): {rho:.4f}")
        if rho < 0:
            print("  WARNING: negative rank correlation -- possible wrong class/sign.")

    baseline_scalar = float(np.nanmean(base_all))
    pd.DataFrame(shap_all, columns=feature_names).to_csv(
        output_dir / "shap_values.csv", index=False)
    X_df.to_csv(output_dir / "X_test_explained.csv", index=False)
    (output_dir / "shap_baseline.txt").write_text(str(baseline_scalar))
    pd.DataFrame({"baseline": base_all}).to_csv(output_dir / "baselines.csv", index=False)
    pd.DataFrame({"oof_pred": pred_all}).to_csv(output_dir / "oof_pred.csv", index=False)
    pd.DataFrame({"fold": folds}).to_csv(output_dir / "fold_id.csv", index=False)
    (output_dir / "feature_names.txt").write_text("\n".join(feature_names))
    (output_dir / "shap_space.txt").write_text(space)
    if meta_df is not None:
        meta_df.to_csv(output_dir / "meta.csv", index=False)

    print(f"\nSaved OOF outputs to {output_dir}/  ({space} space, baseline={baseline_scalar:.4f})")
```

- [ ] **Step 5: Branch `main()` into OOF mode**

In `scripts/SHAP_tabpfn_local.py`, at the very top of `main()` (right after `args = parse_args()` and the `output_dir.mkdir(...)` lines), insert the OOF branch. Find:

```python
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
```

and add immediately after it:

```python
    if args.folds is not None:
        run_oof(args, input_dir, output_dir)
        return
```

- [ ] **Step 6: Run the smoke test to verify it passes**

Run (server, tabpfn2 env): `python scripts/test_oof_smoke.py`
Expected: PASS — prints `imputation additivity mean/max: ...` (max < 0.05), `recontext rank corr ...` (> 0.3), and `SMOKE TEST PASSED`.

- [ ] **Step 7: Commit**

```bash
git add scripts/SHAP_tabpfn_local.py scripts/test_oof_smoke.py
git commit -m "feat: add OOF (--folds) mode to SHAP_tabpfn_local.py with smoke test"
```

---

## Task 2: Notebook 10 — data swap, 14 features, fold/metadata export

**Files:**
- Modify: `analysis/primary/10_TabPFN.qmd:36-50` (data prep chunk)
- Modify: `analysis/primary/10_TabPFN.qmd` (add export chunk after the `split` chunk, ~line 70)

- [ ] **Step 1: Swap data source and build the 14-feature `model_data`**

Replace `analysis/primary/10_TabPFN.qmd` lines 36-50 (from `data <- read_xlsx(...)` through `summary(model_data)`), which currently read:

```r
data <- read_xlsx(here('data', 'raw', 'PCAWG_variables_PCA_woK.xlsx'))

clean_data <- clean_pcawg_data(data)

model_data <- clean_data %>%
  select(-c(1, 2)) %>%
  filter(complete.cases(.)) %>%
  mutate(
    alt_status = factor(
      alt_status,
      levels = c("ALT-low", "ALT-high")
    )
  )

summary(model_data)
```

with:

```r
data <- read_xlsx(here('data', 'processed', 'PCAWG_primary.xlsx'))

# 14 telomere features, now INCLUDING tf_blood_rate
features14 <- c(
  "tf_primary_rate", "tf_blood_rate", "telomere_insertion_rate",
  "telomere_content_log2", "TERT_FPKM",
  "ATAGGG_singleton_dist", "CTAGGG_singleton_dist", "GTAGGG_singleton_dist",
  "TAAGGG_singleton_dist", "TCAGGG_singleton_dist", "TGAGGG_singleton_dist",
  "TTCGGG_singleton_dist", "TTGGGG_singleton_dist", "TTTGGG_singleton_dist"
)

# keep features + label + metadata; drop rows missing any feature (notably TERT_FPKM)
model_full <- data %>%
  select(all_of(features14), alt_status, donor_id, cancer_type, cancer_group) %>%
  filter(complete.cases(across(all_of(features14)))) %>%
  mutate(alt_status = factor(alt_status, levels = c("ALT-low", "ALT-high")))

# model_data keeps the same shape the rest of the notebook expects: features + alt_status
model_data <- model_full %>% select(all_of(features14), alt_status)

summary(model_data)
```

- [ ] **Step 2: Add the CV export chunk**

In `analysis/primary/10_TabPFN.qmd`, immediately after the `split` chunk (after line 70, before the "Here we have split..." prose), add a new chunk:

````markdown
# Export for out-of-fold SHAP (5-fold CV)

```{r export-shap-cv}
library(caret)

# one stratified 5-fold split on alt_status (mirrors notebook 02)
set.seed(42)
fold_list <- createFolds(model_data$alt_status, k = 5)
folds_vec <- integer(nrow(model_data))
for (i in seq_along(fold_list)) folds_vec[fold_list[[i]]] <- i

shap_cv_dir <- here('outputs', 'shap_input_cv')
dir.create(shap_cv_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(model_data %>% select(all_of(features14)),
                 file.path(shap_cv_dir, 'X_full.csv'))
readr::write_csv(data.frame(y = as.integer(model_data$alt_status == "ALT-high")),
                 file.path(shap_cv_dir, 'y_full.csv'))
readr::write_csv(data.frame(fold = folds_vec),
                 file.path(shap_cv_dir, 'folds.csv'))
readr::write_csv(model_full %>%
                   transmute(donor_id, cancer_type, cancer_group, alt_status),
                 file.path(shap_cv_dir, 'meta.csv'))

cat(sprintf("Exported %d samples x %d features to %s\n",
            nrow(model_data), length(features14), shap_cv_dir))
table(folds_vec, model_data$alt_status)
```
````

- [ ] **Step 3: Run the two chunks and verify exports**

In RStudio, run the data-prep chunk and the `export-shap-cv` chunk. Then verify from a terminal:

Run:
```bash
cd /Users/andreapolazzi/Desktop/CBM/project
wc -l outputs/shap_input_cv/X_full.csv outputs/shap_input_cv/y_full.csv \
      outputs/shap_input_cv/folds.csv outputs/shap_input_cv/meta.csv
head -1 outputs/shap_input_cv/X_full.csv
```
Expected: all four files have the same row count (N+1 with header); `X_full.csv` header lists the 14 feature names including `tf_blood_rate`; the `table(folds_vec, alt_status)` print shows all 5 folds populated with both classes.

- [ ] **Step 4: Commit**

```bash
git add analysis/primary/10_TabPFN.qmd
git commit -m "feat: notebook 10 uses PCAWG_primary (incl tf_blood_rate), exports 5-fold CV SHAP inputs"
```

---

## Task 3: PBS launcher `run_shap_full2.sh`

**Files:**
- Create: `scripts/run_shap_full2.sh`

- [ ] **Step 1: Write the launcher**

Create `scripts/run_shap_full2.sh`:

```bash
#!/bin/bash
#
# PRODUCTION OOF SHAP run on the chiron HPC GPU queue (16 cores).
# Both explainers, 5-fold out-of-fold, all samples, budget 256.
#   imputation     -> output_imputation_full2/   (probability space)
#   recontextualiz -> output_recontext_full2/    (logit-like space)
#
# Logging is internal (tee to a timestamped file) -- do NOT rely on #PBS -o with
# ${PBS_JOBID} (it does not expand in directive lines and the log is lost).
#
# Submit (GPU, 16 cores):
#   qsubmit.pl -g 1 -n 16 -s /ngs/iflores/andrea/run_shap_full2.sh

set -u

BASE=/ngs/iflores/andrea
mkdir -p "$BASE/logs" "$BASE/output_imputation_full2" "$BASE/output_recontext_full2"

TS=$(date +%Y%m%d_%H%M%S)
LOG="$BASE/logs/shap_full2_${TS}.log"
exec > >(tee -a "$LOG") 2>&1

echo "######################################################################"
echo "# SHAP OOF full2 (both explainers)  |  started $(date)"
echo "# log file: $LOG"
echo "######################################################################"

# env activation MUST precede python (else ModuleNotFoundError: tabpfn)
source /ngs/software/micromamba/env.sh
micromamba activate /ngs/software/conda/envs/tabpfn2
cd "$BASE"

echo "host        : $(hostname)"
echo "which python: $(which python)"
nvidia-smi || echo "(no nvidia-smi / no GPU visible)"
echo "-----------------------------------"

run_step () {
  local name="$1"; shift
  echo
  echo "=================================================================="
  echo "START $name  |  $(date)"
  echo "cmd: python -u SHAP_tabpfn_local.py $*"
  echo "=================================================================="
  python -u SHAP_tabpfn_local.py "$@"
  local rc=$?
  echo "------------------------------------------------------------------"
  echo "END   $name  |  exit code = $rc  |  $(date)"
  echo "=================================================================="
  return $rc
}

run_step "imputation-oof" \
  --explainer   imputation \
  --input-dir   "$BASE/shap_input_cv" \
  --output-dir  "$BASE/output_imputation_full2" \
  --folds       "$BASE/shap_input_cv/folds.csv" \
  --budget      256 \
  --n-jobs      16 \
  --class-index 1 \
  --seed        42
IMP_RC=$?

run_step "recontext-oof" \
  --explainer   recontextualization \
  --input-dir   "$BASE/shap_input_cv" \
  --output-dir  "$BASE/output_recontext_full2" \
  --folds       "$BASE/shap_input_cv/folds.csv" \
  --budget      256 \
  --n-jobs      16 \
  --class-index 1 \
  --seed        42
RECON_RC=$?

echo
echo "######################################################################"
echo "# SUMMARY"
echo "#   imputation-oof exit code : $IMP_RC   (0 = OK)"
echo "#   recontext-oof  exit code : $RECON_RC (0 = OK)"
echo "#   finished $(date)  |  log: $LOG"
echo "######################################################################"
```

- [ ] **Step 2: Lint the script**

Run: `bash -n scripts/run_shap_full2.sh`
Expected: no output (syntax OK).

- [ ] **Step 3: Commit**

```bash
git add scripts/run_shap_full2.sh
git commit -m "feat: add HPC PBS launcher run_shap_full2.sh for OOF SHAP (both explainers)"
```

---

## Task 4: Extend `reconstruct_shapviz.R` for `full2` + metadata

**Files:**
- Modify: `scripts/reconstruct_shapviz.R`

- [ ] **Step 1: Add a `load_meta()` helper and the full2 object lines**

In `scripts/reconstruct_shapviz.R`, after the existing `load_shap <- function(dir) { ... }` definition, add:

```r
# Load the row-aligned metadata (donor_id, cancer_type, cancer_group, alt_status)
# written alongside the OOF outputs. Returns NULL if meta.csv is absent.
load_meta <- function(dir) {
  path <- file.path(dir, "meta.csv")
  if (!file.exists(path)) {
    message("No meta.csv in ", dir)
    return(NULL)
  }
  read.csv(path, check.names = FALSE)
}
```

Then, after the existing `sv_imp` / `sv_rec` lines at the bottom, add:

```r
# --- 5-fold out-of-fold, blood-inclusive (full2) ------------------------
# Each row is explained out-of-fold; baseline used here is the mean of the
# per-row (per-fold) baselines. Per-row baselines are in baselines.csv if you
# need exact per-sample reconstruction.
sv_imp2  <- load_shap(here("outputs", "shap_output", "output_imputation_full2"))  # probability
sv_rec2  <- load_shap(here("outputs", "shap_output", "output_recontext_full2"))   # logit-like

# Row-aligned labels for coloring/faceting plots (same order as sv_imp2$X rows):
meta2 <- load_meta(here("outputs", "shap_output", "output_imputation_full2"))
# e.g. sv_dependence(sv_imp2, v = "tf_blood_rate", color_var = meta2$alt_status)
```

- [ ] **Step 2: Verify the file parses (after outputs are downloaded)**

This step is run once the `full2` output dirs are downloaded into `outputs/shap_output/`. Until then it will message about missing dirs, which is expected.

Run (in R, from project root): `source("scripts/reconstruct_shapviz.R")`
Expected (once outputs present): `sv_imp2` and `sv_rec2` are `shapviz` objects with N rows × 14 features; `meta2` is a data.frame with N rows and columns `donor_id, cancer_type, cancer_group, alt_status`; `nrow(meta2) == nrow(sv_imp2$X)`.

- [ ] **Step 3: Commit**

```bash
git add scripts/reconstruct_shapviz.R
git commit -m "feat: reconstruct_shapviz.R loads full2 OOF objects + row-aligned metadata"
```

---

## Task 5: Deploy & run on HPC (manual)

This task has no code; it is the run procedure. Do it after Tasks 1-4 are committed and the smoke test (Task 1) passed.

- [ ] **Step 1: Copy the updated script, launcher, and CV inputs to the HPC**

```bash
scp scripts/SHAP_tabpfn_local.py scripts/run_shap_full2.sh iflores@chiron1:/ngs/iflores/andrea/
scp -r outputs/shap_input_cv iflores@chiron1:/ngs/iflores/andrea/
```

- [ ] **Step 2: Submit the GPU job (16 cores)**

On chiron1:
```bash
qsubmit.pl -g 1 -n 16 -s /ngs/iflores/andrea/run_shap_full2.sh
```

- [ ] **Step 3: Monitor**

```bash
qstat -u iflores
tail -f /ngs/iflores/andrea/logs/shap_full2_*.log
```
Expected per fold: `--- fold i: train=... explain=... ---`; imputation prints per-fold `baseline ... vs expected ...` with small `|diff|`; final `additivity error` tiny for imputation and `rank corr ... ~ +0.9x` for recontext; `SUMMARY` with both exit codes 0.

- [ ] **Step 4: Download the outputs for R**

On the local mac:
```bash
mkdir -p outputs/shap_output
scp -r iflores@chiron1:/ngs/iflores/andrea/output_imputation_full2 outputs/shap_output/
scp -r iflores@chiron1:/ngs/iflores/andrea/output_recontext_full2  outputs/shap_output/
```

- [ ] **Step 5: Reconstruct in R**

In R from the project root: `source("scripts/reconstruct_shapviz.R")`, then confirm `nrow(meta2) == nrow(sv_imp2$X)` and that `colnames(sv_imp2$X)` contains `tf_blood_rate`.

---

## Self-review notes

- Spec section 1 (data/features) → Task 2 Step 1. Section 2 (folds/export) → Task 2 Step 2. Section 3 (Python OOF) → Task 1. Section 4 (R helper) → Task 4. Runtime/launcher → Task 3 + Task 5. All spec sections covered.
- Filenames `shap_values.csv`, `X_test_explained.csv`, `shap_baseline.txt`, `shap_space.txt` are produced by `run_oof` (Task 1 Step 4) and consumed unchanged by `load_shap` (existing) — type/name consistency holds.
- `--folds` is defined (Task 1 Step 3) and used by both the smoke test (Task 1 Step 1) and the launcher (Task 3) consistently.
