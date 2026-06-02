#!/usr/bin/env python3
"""
HPC script: fit TabPFN locally and compute SHAP values via shapiq.

Two explainer paradigms are available via --explainer:

  imputation (default)
    Uses shapiq.TabularExplainer on an explicit predict_proba[:, class] function.
    Output is in PROBABILITY SPACE and the target class is unambiguous. The script
    self-verifies and ABORTS before writing if the class/scale is wrong:
      * baseline must equal mean(predict_proba[:, class]) over the background, and
      * baseline + sum(SHAP_row) must reconstruct each row's predicted probability.

  recontextualization
    Uses tabpfn_extensions get_tabpfn_explainer (remove-and-recontextualize). This is
    the paradigm the TabPFN authors recommend and avoids off-manifold imputation, but
    output is in the model's RAW (logit-like) space and class handling is delegated to
    the extension. The script prints diagnostics to identify the explained class but
    does NOT hard-abort, since the output space is not probability.

Input files (--input-dir):
  X_train.csv   - training feature matrix (rows = samples, cols = features)
  X_test.csv    - test feature matrix to explain
  y_train.csv   - training labels (integer 0/1 in first column)

Output files (--output-dir):
  shap_values.csv       - SHAP matrix (n_test x n_features)
  X_test_explained.csv  - feature values for explained rows
  shap_baseline.txt     - scalar base value
  feature_names.txt     - ordered feature names, one per line
  shap_space.txt        - "probability" or "logit" (which space the values are in)

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
import shapiq
from tabpfn import TabPFNClassifier
from tabpfn_extensions.interpretability.shapiq import get_tabpfn_explainer

# Tolerances for the imputation-mode self-verification step.
BASELINE_TOL = 0.02     # |baseline - mean(predict_proba)| must be below this
ADDITIVITY_WARN = 0.05  # warn (do not fail) if mean reconstruction error exceeds this


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--input-dir", default="outputs/shap_input",
                   help="Directory containing X_train.csv, X_test.csv, y_train.csv")
    p.add_argument("--output-dir", default="outputs/shap_output",
                   help="Directory where output files are written")
    p.add_argument("--explainer", choices=["imputation", "recontextualization"],
                   default="imputation",
                   help="Explainer paradigm. Default: imputation (probability space, "
                        "self-verified).")
    p.add_argument("--budget", type=int, default=256,
                   help="Shapley budget: model evaluations per sample. "
                        "Higher = more accurate but slower. Default: 256")
    p.add_argument("--n-explain", type=int, default=None,
                   help="Explain only the first N rows of X_test. Default: all rows")
    p.add_argument("--n-jobs", type=int, default=-1,
                   help="Parallel jobs for explain_X (-1 = all cores). Default: -1")
    p.add_argument("--class-index", type=int, default=1,
                   help="Class to explain (1 = ALT-high minority class). Default: 1")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for the Shapley approximator. Default: 42")
    return p.parse_args()


def load_data(input_dir):
    paths = {
        "X_train": input_dir / "X_train.csv",
        "X_test": input_dir / "X_test.csv",
        "y_train": input_dir / "y_train.csv",
    }
    for name, path in paths.items():
        if not path.exists():
            print(f"ERROR: missing file: {path}", file=sys.stderr)
            sys.exit(1)
    X_train_df = pd.read_csv(paths["X_train"])
    X_test_df = pd.read_csv(paths["X_test"])
    y_train = pd.read_csv(paths["y_train"]).iloc[:, 0].to_numpy().astype(int)
    if list(X_test_df.columns) != list(X_train_df.columns):
        print("ERROR: X_test columns do not match X_train columns.", file=sys.stderr)
        sys.exit(1)
    return X_train_df, X_test_df, y_train


def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    X_train_df, X_test_df, y_train = load_data(input_dir)
    feature_names = list(X_train_df.columns)
    n_features = len(feature_names)
    X_train = X_train_df.to_numpy().astype(float)
    X_test = X_test_df.to_numpy().astype(float)

    class_index = args.class_index
    n_classes = int(y_train.max()) + 1
    if not (0 <= class_index < n_classes):
        print(f"ERROR: --class-index {class_index} out of range for {n_classes} classes.",
              file=sys.stderr)
        sys.exit(1)

    if args.n_explain is not None:
        X_test = X_test[:args.n_explain]
        X_test_df = X_test_df.iloc[:args.n_explain].reset_index(drop=True)

    print(f"Explainer: {args.explainer}")
    print(f"X_train : {X_train.shape}")
    print(f"X_test  : {X_test.shape}  (rows to explain)")
    print(f"Features: {n_features}")
    print(f"Target class index: {class_index} "
          f"(prevalence in train: {(y_train == class_index).mean():.4f})")

    # --- Fit TabPFN (local, no API token required) ---
    print("\nFitting TabPFN...", flush=True)
    t0 = time.time()
    clf = TabPFNClassifier()
    clf.fit(X_train, y_train)
    print(f"Done in {time.time() - t0:.1f}s")

    classes = list(clf.classes_)
    if class_index not in classes:
        print(f"ERROR: class label {class_index} not in model classes {classes}.",
              file=sys.stderr)
        sys.exit(1)
    proba_col = classes.index(class_index)
    print(f"Model classes: {classes} -> target predict_proba column {proba_col}")

    def proba_predict(Z):
        return clf.predict_proba(Z)[:, proba_col]

    # --- Build explainer ---
    if args.explainer == "imputation":
        space = "probability"
        expected_baseline = float(proba_predict(X_train).mean())
        print(f"Expected baseline (mean predict_proba over train): {expected_baseline:.4f}")
        print(f"\nBuilding shapiq.TabularExplainer (budget={args.budget})...", flush=True)
        explainer = shapiq.TabularExplainer(
            model=proba_predict,
            data=X_train,
            index="SV",        # standard Shapley values (first order)
            max_order=1,
            imputer="marginal",
        )
    else:  # recontextualization
        space = "logit"
        expected_baseline = None
        print(f"\nBuilding TabPFN recontextualization explainer (budget={args.budget})...",
              flush=True)
        explainer = get_tabpfn_explainer(
            model=clf,
            data=X_train,
            labels=y_train,
            index="SV",
            max_order=1,
            class_index=class_index,
        )

    # --- Compute Shapley values for all test rows ---
    print(f"Computing Shapley values (n_jobs={args.n_jobs})...", flush=True)
    t0 = time.time()
    iv_list = explainer.explain_X(
        X_test,
        budget=args.budget,
        n_jobs=args.n_jobs,
        random_state=args.seed,
        verbose=True,
    )
    print(f"Done in {time.time() - t0:.1f}s")

    shap_matrix = np.array(
        [[float(iv[(j,)]) for j in range(n_features)] for iv in iv_list]
    )
    baselines = np.array([float(iv.baseline_value) for iv in iv_list])
    baseline = float(baselines.mean())

    # --- Verification / diagnostics ---
    print(f"\nVerification ({space} space):")
    print(f"  baseline (from SHAP)          : {baseline:.4f}")
    if args.explainer == "imputation":
        baseline_err = abs(baseline - expected_baseline)
        print(f"  expected (mean predict_proba) : {expected_baseline:.4f}")
        print(f"  |difference|                  : {baseline_err:.4f} (tol {BASELINE_TOL})")
        if baseline_err > BASELINE_TOL:
            print("ERROR: baseline does not match the target class probability. SHAP values "
                  "would be for the WRONG class. Aborting without writing output.",
                  file=sys.stderr)
            sys.exit(2)
        actual = proba_predict(X_test)
        recon = baselines + shap_matrix.sum(axis=1)
        add_err = np.abs(recon - actual)
        print(f"  additivity error (mean/max)   : {add_err.mean():.4g} / {add_err.max():.4g}")
        if add_err.mean() > ADDITIVITY_WARN:
            print(f"  WARNING: mean additivity error {add_err.mean():.4g} exceeds "
                  f"{ADDITIVITY_WARN}. Consider increasing --budget.")
        print("  Class check PASSED.")
    else:
        # Diagnostics to identify which class/space the output represents.
        p_target = float(proba_predict(X_train).mean())
        p_other = 1.0 - p_target
        eps = 1e-9
        logit_target = float(np.log((p_target + eps) / (p_other + eps)))
        print(f"  mean predict_proba target class : {p_target:.4f}")
        print(f"  logit(mean prob target class)   : {logit_target:.4f}")
        print(f"  logit(mean prob OTHER class)    : {-logit_target:.4f}")
        print("  (Compare baseline above to these to confirm the explained class/sign.)")

    # --- Save outputs ---
    pd.DataFrame(shap_matrix, columns=feature_names).to_csv(
        output_dir / "shap_values.csv", index=False
    )
    X_test_df.to_csv(output_dir / "X_test_explained.csv", index=False)
    (output_dir / "shap_baseline.txt").write_text(str(baseline))
    (output_dir / "feature_names.txt").write_text("\n".join(feature_names))
    (output_dir / "shap_space.txt").write_text(space)

    print(f"\nSaved to: {output_dir}/  (values in {space} space)")
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
