import os
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tabpfn_client
from tabpfn_client import TabPFNClassifier
from tabpfn_extensions import interpretability


def main() -> None:
    token = os.getenv("TABPFN_TOKEN")
    if not token:
        raise RuntimeError("TABPFN_TOKEN non trovato nelle variabili d'ambiente.")

    tabpfn_client.set_access_token(token)

    base_dir = Path(__file__).resolve().parent.parent
    input_dir = base_dir / "outputs" / "shap_input"

    output_dir = base_dir / "analysis" / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    x_train_path = input_dir / "X_train.csv"
    x_test_path = input_dir / "X_test_small.csv"
    y_train_path = input_dir / "y_train.csv"

    if not x_test_path.exists():
        x_test_path = input_dir / "X_test.csv"

    if not x_train_path.exists() or not x_test_path.exists() or not y_train_path.exists():
        raise FileNotFoundError("File input SHAP mancanti in outputs/shap_input/")

    X_train = pd.read_csv(x_train_path)
    X_test = pd.read_csv(x_test_path)
    y_train = pd.read_csv(y_train_path).iloc[:, 0].to_numpy().astype(int)

    feature_names = list(X_train.columns)

    X_train_np = X_train.to_numpy().astype(float)
    X_test_np = X_test.to_numpy().astype(float)

    print(f"X_train shape: {X_train_np.shape}")
    print(f"X_test shape:  {X_test_np.shape}")

    clf = TabPFNClassifier()
    print("Fitting TabPFN...")
    clf.fit(X_train_np, y_train)

    print("Starting SHAP...")
    t0 = time.time()

    shap_values = interpretability.shap.get_shap_values(
        estimator=clf,
        test_x=X_test_np,
        attribute_names=feature_names,
        algorithm="permutation",
    )

    print(f"SHAP completed in {time.time() - t0:.2f} seconds")

    interpretability.shap.plot_shap(shap_values)

    fig_nums = plt.get_fignums()
    for i, num in enumerate(fig_nums):
        fig = plt.figure(num)
        fig.savefig(output_dir / f"tabpfn_shap{i+1}.png", bbox_inches="tight", dpi=150)

    plt.close("all")

    print(f"Saved SHAP plots to: {output_dir}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise
