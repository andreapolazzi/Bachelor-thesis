#!/usr/bin/env python3
"""
Diagnostic: understand exactly what get_tabpfn_explainer (recontextualization)
computes, so we can make it explain the CORRECT class.

It answers four questions and prints them clearly:

  Q1. SOURCE  - the actual source of get_tabpfn_explainer in THIS install,
                so we see how `class_index` is used and what is returned.
  Q2. MODEL   - what does explainer.model(X) output? shape, values, and how it
                compares to clf.predict_proba(X). Tells us the space (proba vs
                logit) and which class / sign.
  Q3. ADDIT.  - per row: baseline + sum(shap) vs explainer.model(row) vs the
                true class-1 probability. This is the real correctness test.
  Q4. OBJECT  - the explainer type and its model-related attributes.

Run on carmela:
  cd /ngs/iflores/andrea
  python -u diag_recontext.py --input-dir /ngs/iflores/andrea/input | tee logs/diag_recontext.log
"""
import argparse
import inspect
from pathlib import Path

import numpy as np
import pandas as pd
import shapiq
from tabpfn import TabPFNClassifier
from tabpfn_extensions.interpretability import shapiq as tex_shapiq


def sig(x):
    return 1.0 / (1.0 + np.exp(-np.asarray(x, dtype=float)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", default="/ngs/iflores/andrea/input")
    ap.add_argument("--class-index", type=int, default=1)
    ap.add_argument("--n-rows", type=int, default=5)
    ap.add_argument("--budget", type=int, default=64)
    args = ap.parse_args()

    d = Path(args.input_dir)
    X_train = pd.read_csv(d / "X_train.csv").to_numpy().astype(float)
    X_test = pd.read_csv(d / "X_test.csv").to_numpy().astype(float)[: args.n_rows]
    y_train = pd.read_csv(d / "y_train.csv").iloc[:, 0].to_numpy().astype(int)
    n_feat = X_train.shape[1]
    ci = args.class_index

    print("=" * 70)
    print("Q1. SOURCE of get_tabpfn_explainer (this installed version)")
    print("=" * 70)
    print("file:", inspect.getsourcefile(tex_shapiq.get_tabpfn_explainer))
    try:
        print(inspect.getsource(tex_shapiq.get_tabpfn_explainer))
    except Exception as e:
        print("could not get source:", e)

    print("\nFitting TabPFN...")
    clf = TabPFNClassifier()
    clf.fit(X_train, y_train)
    classes = list(clf.classes_)
    proba_col = classes.index(ci)
    print(f"classes={classes}  target class {ci} -> predict_proba column {proba_col}")

    print("\n" + "=" * 70)
    print("Q4. OBJECT: build explainer, inspect type & attributes")
    print("=" * 70)
    explainer = tex_shapiq.get_tabpfn_explainer(
        model=clf,
        data=X_train,
        labels=y_train,
        index="SV",
        max_order=1,
        class_index=ci,
    )
    print("type:", type(explainer))
    print("attrs:", [a for a in dir(explainer) if not a.startswith("__")])
    model_fn = getattr(explainer, "model", None)
    print("explainer.model:", model_fn, "type:", type(model_fn))
    # shapiq stores the prediction function; how it maps to classes matters.
    for attr in ("predict_function", "_predict_function", "class_index", "_class_index"):
        if hasattr(explainer, attr):
            print(f"  explainer.{attr} =", getattr(explainer, attr))

    print("\n" + "=" * 70)
    print("Q2. MODEL OUTPUT: explainer.model(X) vs clf.predict_proba(X)")
    print("=" * 70)
    proba = clf.predict_proba(X_test)  # (n, n_classes)
    out = None
    if callable(model_fn):
        try:
            out = np.asarray(model_fn(X_test))
            print("explainer.model(X) shape:", out.shape)
            print("explainer.model(X) values:\n", out)
        except Exception as e:
            print("calling explainer.model(X) failed:", e)
    print("\nclf.predict_proba(X):\n", proba)
    print("\nFor reference, per row:")
    print(f"{'row':>3} {'proba[:,0]':>10} {'proba[:,1]':>10} "
          f"{'logit(c1)':>10} {'model_out':>12}")
    for i in range(len(X_test)):
        p1 = proba[i, proba_col]
        logit_c1 = np.log((p1 + 1e-9) / (1 - p1 + 1e-9))
        mo = "" if out is None else np.array2string(out[i], precision=4)
        print(f"{i:>3} {proba[i,0]:>10.4f} {proba[i,1]:>10.4f} "
              f"{logit_c1:>10.4f} {mo:>12}")

    print("\n" + "=" * 70)
    print(f"Q3. ADDITIVITY: explain {len(X_test)} rows (budget={args.budget})")
    print("=" * 70)
    ivs = explainer.explain_X(
        X_test, budget=args.budget, n_jobs=1, random_state=42, verbose=False
    )
    print(f"{'row':>3} {'baseline':>10} {'sum_shap':>10} {'base+sum':>10} "
          f"{'model_out':>12} {'proba[:,1]':>10} {'sig(b+s)':>10}")
    for i, iv in enumerate(ivs):
        s = float(sum(float(iv[(j,)]) for j in range(n_feat)))
        b = float(iv.baseline_value)
        recon = b + s
        mo = "" if out is None else float(np.ravel(out[i])[-1])
        mo_str = "" if out is None else f"{mo:>12.4f}"
        print(f"{i:>3} {b:>10.4f} {s:>10.4f} {recon:>10.4f} "
              f"{mo_str} {proba[i,proba_col]:>10.4f} {sig(recon):>10.4f}")

    print("\nINTERPRETATION GUIDE:")
    print("  * If 'base+sum' matches 'model_out' per row -> additivity holds (good).")
    print("  * If 'model_out' matches 'proba[:,1]' (or its logit) -> correct class.")
    print("  * If 'sig(b+s)' tracks 'proba[:,1]' (low for low-risk rows) -> correct class,")
    print("    and the baseline is just the recontextualized reference (NOT logit(mean p)).")
    print("  * If 'sig(b+s)' is ~1 - proba[:,1] -> it's explaining the OTHER class (sign flip).")


if __name__ == "__main__":
    main()
