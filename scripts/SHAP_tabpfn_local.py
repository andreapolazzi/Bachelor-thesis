#!/usr/bin/env python3
"""
HPC script: fit TabPFN locally and compute SHAP values via shapiq.

Input files (--input-dir):
  X_train.csv   — training feature matrix (rows = samples, cols = features)
  X_test.csv    — test feature matrix to explain
  y_train.csv   — training labels (integer 0/1 in first column)

Output files (--output-dir):
  shap_values.csv       — SHAP matrix (n_test x n_features)
  X_test_explained.csv  — feature values for explained rows
  shap_baseline.txt     — scalar base value (mean expected prediction)
  feature_names.txt     — ordered feature names, one per line

Reconstruct in R:
  sv <- shapviz(
    object   = as.matrix(read.csv("shap_values.csv")),
    X        = read.csv("X_test_explained.csv"),
    baseline = as.numeric(readLines("shap_baseline.txt"))
  )
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from tabpfn import TabPFNClassifier
from tabpfn_extensions.interpretability.shapiq import (
    get_tabpfn_explainer,
    get_tabpfn_imputation_explainer,
)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir",  default="outputs/shap_input",
                   help="Directory containing X_train.csv, X_test.csv, y_train.csv")
    p.add_argument("--output-dir", default="outputs/shap_output",
                   help="Directory where output files are written")
    p.add_argument("--budget", type=int, default=256,
                   help="Shapley budget: model evaluations per sample. "
                        "Higher = more accurate but slower. Default: 256")
    p.add_argument("--n-explain", type=int, default=None,
                   help="Explain only the first N rows of X_test. Default: all rows")
    p.add_argument("--n-jobs", type=int, default=-1,
                   help="Parallel jobs for explain_X (-1 = all cores). Default: -1")
    p.add_argument("--explainer", choices=["imputation", "recontextualization"],
                   default="recontextualization",
                   help="Explainer type. 'recontextualization' is theoretically correct for "
                        "TabPFN but slower; 'imputation' is faster. Default: recontextualization")
    p.add_argument("--class-index", type=int, default=1,
                   help="Class index for SHAP values (1 = ALT-high). Default: 1")
    return p.parse_args()


def main():
    args = parse_args()
    input_dir  = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Load data ---
    paths = {
        "X_train": input_dir / "X_train.csv",
        "X_test":  input_dir / "X_test.csv",
        "y_train": input_dir / "y_train.csv",
    }
    for name, path in paths.items():
        if not path.exists():
            print(f"ERROR: missing file: {path}", file=sys.stderr)
            sys.exit(1)

    X_train_df = pd.read_csv(paths["X_train"])
    X_test_df  = pd.read_csv(paths["X_test"])
    y_train    = pd.read_csv(paths["y_train"]).iloc[:, 0].to_numpy().astype(int)

    feature_names = list(X_train_df.columns)
    n_features    = len(feature_names)
    X_train = X_train_df.to_numpy().astype(float)
    X_test  = X_test_df.to_numpy().astype(float)

    if args.n_explain is not None:
        X_test    = X_test[:args.n_explain]
        X_test_df = X_test_df.iloc[:args.n_explain].reset_index(drop=True)

    print(f"X_train : {X_train.shape}")
    print(f"X_test  : {X_test.shape}  (rows to explain)")
    print(f"Features: {n_features}")

    # --- Fit TabPFN (local, no API token required) ---
    print("\nFitting TabPFN...", flush=True)
    t0 = time.time()
    clf = TabPFNClassifier()
    clf.fit(X_train, y_train)
    print(f"Done in {time.time() - t0:.1f}s")

    # --- Build explainer ---
    # index="SV" + max_order=1 gives standard first-order Shapley values (equivalent to SHAP)
    # recontextualization: theoretically correct for TabPFN (refits on each coalition subset)
    # imputation: faster but less principled for TabPFN
    print(f"\nBuilding explainer (type={args.explainer}, budget={args.budget}, "
          f"class={args.class_index})...", flush=True)
    if args.explainer == "recontextualization":
        explainer = get_tabpfn_explainer(
            model=clf,
            data=X_train,
            labels=y_train,
            index="SV",
            max_order=1,
            class_index=args.class_index,
        )
    else:
        explainer = get_tabpfn_imputation_explainer(
            model=clf,
            data=X_train,
            index="SV",
            max_order=1,
            class_index=args.class_index,
        )

    # --- Compute Shapley values for all test rows ---
    print(f"Computing Shapley values (n_jobs={args.n_jobs}, verbose=True)...", flush=True)
    t0 = time.time()
    iv_list = explainer.explain_X(
        X_test,
        budget=args.budget,
        n_jobs=args.n_jobs,
        verbose=True,
    )
    print(f"Done in {time.time() - t0:.1f}s")

    # --- Assemble SHAP matrix from InteractionValues list ---
    # Each iv in iv_list is an InteractionValues object; iv[(feat_idx,)] gives
    # the first-order Shapley value for that feature.
    shap_matrix = np.array(
        [[float(iv[(j,)]) for j in range(n_features)] for iv in iv_list]
    )
    baseline = float(np.mean([iv.baseline_value for iv in iv_list]))

    # --- Save outputs ---
    pd.DataFrame(shap_matrix, columns=feature_names).to_csv(
        output_dir / "shap_values.csv", index=False
    )
    X_test_df.to_csv(output_dir / "X_test_explained.csv", index=False)
    (output_dir / "shap_baseline.txt").write_text(str(baseline))
    (output_dir / "feature_names.txt").write_text("\n".join(feature_names))

    print(f"\nSaved to: {output_dir}/")
    print(f"  shap_values.csv       {shap_matrix.shape}")
    print(f"  X_test_explained.csv  {X_test.shape}")
    print(f"  shap_baseline.txt     {baseline:.4f}")
    print(f"  feature_names.txt     {n_features} features")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise
