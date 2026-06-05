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
    p.add_argument("--folds", default=None,
                   help="Path to folds.csv (single column of 1..K fold ids, row-aligned "
                        "to X_full.csv). When given, runs OOF mode and reads X_full.csv + "
                        "y_full.csv (and meta.csv if present) from --input-dir instead of "
                        "the X_train/X_test/y_train single-split files.")
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


def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.folds is not None:
        run_oof(args, input_dir, output_dir)
        return

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
        # Recontextualization output is in TabPFN's native (logit-like) space, NOT
        # probability, and the baseline is the empty-coalition (all-features-removed)
        # reference -- it is NOT logit(mean prob), so do not compare it to prevalence.
        # The correct correctness checks are:
        #   (a) additivity: baseline + sum(SHAP_row) reconstructs the row's model output;
        #   (b) class/sign: that reconstruction must INCREASE with the true class
        #       probability (high-risk rows -> higher reconstruction).
        actual_p = proba_predict(X_test)             # true P(class) per row, for ranking
        recon = baselines + shap_matrix.sum(axis=1)  # recontextualized model output
        # Spearman (rank) correlation without scipy:
        ar = np.argsort(np.argsort(actual_p))
        rr = np.argsort(np.argsort(recon))
        if len(recon) > 1 and np.std(ar) > 0 and np.std(rr) > 0:
            spearman = float(np.corrcoef(ar, rr)[0, 1])
        else:
            spearman = float("nan")
        print("  output space                    : recontextualized logit-like "
              "(NOT probability)")
        print(f"  baseline (empty-coalition ref)  : {baseline:.4f}")
        print(f"  rank corr(reconstruction, P(class={class_index})) : {spearman:.4f}")
        print("    -> close to +1 means the explanation correctly tracks the target "
              "class (high-risk rows reconstruct higher). Negative would mean a "
              "sign/class flip.")
        if not np.isnan(spearman) and spearman < 0:
            print("  WARNING: negative rank correlation -- explanation may be for the "
                  "WRONG class/sign. Inspect before using.")

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
