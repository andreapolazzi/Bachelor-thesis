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


def make_inputs(d: Path, n=60, n_feat=4, k=5, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, n_feat))
    # signal: positive class when feature 0 is high (strong, so even a tiny-budget
    # recontextualization run tracks it on this toy data)
    logits = 4.0 * X[:, 0] - 1.0
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
        # Smoke test only checks the pipeline runs and there is no sign/class flip
        # (a flip would give a strongly NEGATIVE correlation). Statistical fidelity
        # is validated on the real data, where this correlation is ~0.95.
        assert rho > 0.1, "recontext reconstruction does not track prediction (possible sign flip)"

    print("SMOKE TEST PASSED")


if __name__ == "__main__":
    main()
