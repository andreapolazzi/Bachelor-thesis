#!/usr/bin/env python
"""
Out-of-fold TabPFN predictions for the model comparisons (classification + regression).

Mirrors the shared 5-fold x 5-repeat stratified CV (scripts/make_folds.R; notebooks
10, 20, 21): for each of the 25 (rep, fold) test partitions a fresh TabPFN model is
fit on the other folds and predicts the held-out rows. The result feeds the fair
comparisons, so TabPFN is scored with the same folds and the same metric code as the
in-R models.

    --task classification -> TabPFNClassifier, yhat = P(positive class = ALT-high)
    --task regression     -> TabPFNRegressor,  yhat = predicted target (log1p space)

Designed for the carmela standalone GPU box (local tabpfn, no API token), same env as
the SHAP scripts (see scripts/SHAP_tabpfn_local.py / run_shap_reg.sh):
    source /ngs/software/micromamba/env.sh
    micromamba activate /ngs/software/conda/envs/tabpfn2
    cd "$BASE"
    python tabpfn_cv_oof.py --task classification \
        --input-dir shap_input_cmp_cv   --output outputs/tabpfn_cv_oof.csv
    python tabpfn_cv_oof.py --task regression \
        --input-dir cmp_reg_ttf_cv      --output outputs/tabpfn_reg_ttf_oof.csv
    python tabpfn_cv_oof.py --task regression \
        --input-dir cmp_reg_btf_cv      --output outputs/tabpfn_reg_btf_oof.csv

Inputs (written by the .qmd export chunks, in --input-dir):
    X_full.csv      n x p design matrix (column order preserved)
    y_full.csv      column `y` (0/1 for classification, continuous for regression)
    folds_5x5.csv   columns rep1..rep5, each the fold id (1..k) of every row

Output (--output): long CSV with columns rep, fold, row, yhat, y
"""

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd


def score(oof, task):
    """Per-fold + summary metrics, matching clf_metrics.R / reg_metrics.R so the
    TabPFN rows are directly comparable to the in-R models."""
    import numpy as np
    rows = []
    for (rep, fold), g in oof.groupby(["rep", "fold"]):
        y = g["y"].to_numpy(); yh = g["yhat"].to_numpy()
        if task == "classification":
            from sklearn.metrics import roc_auc_score, average_precision_score
            yb = y.astype(int); pred = (yh >= 0.5).astype(int)
            tp = int(((pred == 1) & (yb == 1)).sum()); fp = int(((pred == 1) & (yb == 0)).sum())
            tn = int(((pred == 0) & (yb == 0)).sum()); fn = int(((pred == 0) & (yb == 1)).sum())
            sens = tp / (tp + fn) if (tp + fn) else float("nan")
            spec = tn / (tn + fp) if (tn + fp) else float("nan")
            prec = tp / (tp + fp) if (tp + fp) else float("nan")
            f1 = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else float("nan")
            rows.append({"rep": rep, "fold": fold,
                         "auc": roc_auc_score(yb, yh),
                         "pr_auc": average_precision_score(yb, yh),
                         "sensitivity": sens, "specificity": spec,
                         "precision": prec, "f1": f1,
                         "balanced_accuracy": 0.5 * (sens + spec)})
        else:
            from scipy.stats import spearmanr
            rows.append({"rep": rep, "fold": fold,
                         "rmse": float(np.sqrt(np.mean((y - yh) ** 2))),
                         "mae": float(np.mean(np.abs(y - yh))),
                         "spearman": float(spearmanr(y, yh)[0])})
    per = pd.DataFrame(rows)
    cols = [c for c in per.columns if c not in ("rep", "fold")]
    return per, per[cols].agg(["mean", "std"]).T


def make_model(task, n_jobs, seed):
    if task == "classification":
        from tabpfn import TabPFNClassifier
        ctor = TabPFNClassifier
    else:
        from tabpfn import TabPFNRegressor
        ctor = TabPFNRegressor
    try:
        return ctor(n_jobs=n_jobs, random_state=seed)
    except TypeError:          # older signatures
        return ctor()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--task", choices=["classification", "regression"], required=True)
    ap.add_argument("--input-dir", required=True,
                    help="dir with X_full.csv, y_full.csv, folds_5x5.csv")
    ap.add_argument("--output", required=True, help="output CSV path")
    ap.add_argument("--n-jobs", type=int, default=4)
    ap.add_argument("--seed", type=int, default=21)
    args = ap.parse_args()

    in_dir = Path(args.input_dir)
    X = pd.read_csv(in_dir / "X_full.csv").to_numpy(dtype=float)
    y = pd.read_csv(in_dir / "y_full.csv")["y"].to_numpy()
    y = y.astype(int) if args.task == "classification" else y.astype(float)
    folds = pd.read_csv(in_dir / "folds_5x5.csv")
    rep_cols = list(folds.columns)

    print(f"task={args.task}  X={X.shape}  y={y.shape}  reps={len(rep_cols)}", flush=True)

    records = []
    for r_idx, rcol in enumerate(rep_cols, start=1):
        fold_ids = folds[rcol].to_numpy(dtype=int)
        for f in sorted(np.unique(fold_ids)):
            te = np.where(fold_ids == f)[0]
            tr = np.where(fold_ids != f)[0]

            t0 = time.time()
            model = make_model(args.task, args.n_jobs, args.seed)
            model.fit(X[tr], y[tr])
            if args.task == "classification":
                yhat = model.predict_proba(X[te])[:, 1]
            else:
                yhat = model.predict(X[te])

            for row, yh in zip(te, yhat):
                records.append({"rep": r_idx, "fold": int(f), "row": int(row),
                                "yhat": float(yh), "y": float(y[row])})
            print(f"  rep {r_idx} fold {f}: n_test={len(te)} "
                  f"({time.time() - t0:.1f}s)", flush=True)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    oof = pd.DataFrame(records)
    oof.to_csv(out, index=False)
    print(f"wrote {len(oof)} rows -> {out}", flush=True)

    # per-fold + summary metrics (same definitions as the in-R models)
    per_fold, summary = score(oof, args.task)
    summ_path = out.with_name(out.stem + "_summary.csv")
    summary.to_csv(summ_path)
    n_folds = per_fold.shape[0]
    print(f"\n=== {args.task} TabPFN CV — mean ± SD over {n_folds} folds ===", flush=True)
    print(summary.to_string(float_format=lambda v: f"{v:.4f}"), flush=True)
    print(f"wrote summary -> {summ_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
