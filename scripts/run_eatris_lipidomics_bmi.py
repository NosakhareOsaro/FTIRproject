"""
run_eatris_lipidomics_bmi.py

Method comparison (PLS, Ridge, LASSO, elastic net) predicting BMI from
EATRIS-Plus lipidomics feature matrices (125 healthy human adults),
leave-one-out cross-validation, same evaluation discipline as the DGRP
scripts (run_regularised_regression.py, run_pls_analysis.py): StandardScaler
and hyperparameter selection both fitted inside each training fold only.

Written as a new standalone script rather than adapting the DGRP-specific
ones: this dataset is already one row per subject (no DGRP-line-averaging
step to replicate), phenotype coverage is complete (no missing-value
handling or --sex/--spectral-sex split needed, unlike run_dgrpool_phenotype.py),
and the input is a plain samples x features TSV rather than the .dat
spectral format read by ftir_loader.py. Reusing the FTIR-specific scripts
would mean carrying along DGRP-line and FTIR-file machinery that does
nothing useful here.

Usage:
  .venv/bin/python scripts/run_eatris_lipidomics_bmi.py \\
      phenotype-data/EATRIS_Lipidomics_positive.tsv --condition "Positive mode"

  .venv/bin/python scripts/run_eatris_lipidomics_bmi.py \\
      phenotype-data/EATRIS_Lipidomics_negative.tsv --condition "Negative mode"

  .venv/bin/python scripts/run_eatris_lipidomics_bmi.py \\
      phenotype-data/EATRIS_Lipidomics_combined.tsv --condition "Combined"

Treated as three separate test conditions, the same way
run_dgrpool_phenotype.py is invoked once per Unckless diet condition:
each run is independent and appends its own rows to the summary CSV.

Each run appends one row per method to
results/EATRIS/lipidomics_bmi_summary.csv (created with a header if
absent), so results accumulate across conditions the same way
dgrpool_phenotype_summary.csv does for the DGRP work.
"""

import argparse
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.cross_decomposition import PLSRegression
from sklearn.linear_model import ElasticNetCV, LassoCV, RidgeCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).parent.parent
BMI_PATH = REPO / "phenotype-data" / "EATRIS_BMI.tsv"
OUT_DIR = REPO / "results" / "EATRIS"
SUMMARY_PATH = OUT_DIR / "lipidomics_bmi_summary.csv"
SUMMARY_COLS = [
    "condition", "method", "param", "n_samples", "n_features",
    "cv_r2", "rmse", "spearman_rho", "date_run",
]

RIDGE_GRID = np.logspace(-3, 6, 100)   # RidgeCV: alpha candidates, GCV-selected
N_ALPHAS = 30                            # LassoCV / ElasticNetCV: alpha grid size
L1_RATIOS = [0.5, 0.7, 0.9, 0.95, 1.0]  # ElasticNetCV: L1/L2 mixing sweep
PLS_N_COMPONENTS = [1, 2, 3, 5, 10, 15, 20]

loo = LeaveOneOut()


def metrics(y_true, y_pred):
    r2 = r2_score(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    rho = spearmanr(y_true, y_pred).statistic
    return r2, rmse, rho


def loo_cv(X, y, model_fn):
    """LOO-CV with per-fold StandardScaling. model_fn returns a fresh estimator."""
    y_pred = np.zeros(len(y))
    for train_idx, test_idx in loo.split(X):
        sc = StandardScaler()
        X_tr = sc.fit_transform(X[train_idx])
        X_te = sc.transform(X[test_idx])
        mdl = model_fn()
        mdl.fit(X_tr, y[train_idx])
        y_pred[test_idx] = np.asarray(mdl.predict(X_te)).ravel()
    return y_pred


def main():
    parser = argparse.ArgumentParser(
        description="LOO-CV method comparison (PLS, Ridge, LASSO, elastic net) "
        "predicting BMI from an EATRIS-Plus lipidomics feature matrix."
    )
    parser.add_argument(
        "features_tsv", type=Path,
        help="Path to an EATRIS_Lipidomics_*.tsv file (sample_id + feature columns)",
    )
    parser.add_argument(
        "--condition", required=True,
        help="Label for this run, e.g. 'Positive mode', 'Negative mode', 'Combined'",
    )
    args = parser.parse_args()

    X_df = pd.read_csv(args.features_tsv, sep="\t")
    bmi_df = pd.read_csv(BMI_PATH, sep="\t")

    merged = X_df.merge(bmi_df, on="sample_id", how="inner")
    n_dropped = len(X_df) - len(merged)
    if n_dropped:
        print(f"Note: {n_dropped} sample(s) in {args.features_tsv.name} had no "
              f"matching BMI value and were dropped.")

    feature_cols = [c for c in X_df.columns if c != "sample_id"]
    X = merged[feature_cols].values
    y = merged["BMI"].values.astype(float)
    n_samples, n_features = X.shape

    print("=" * 68)
    print(f"Condition : {args.condition}")
    print(f"File      : {args.features_tsv}")
    print(f"Samples   : {n_samples}   Features: {n_features}")
    print(f"BMI range : {y.min():.1f} - {y.max():.1f}  "
          f"(mean {y.mean():.2f}, SD {y.std():.2f})")
    print("=" * 68)
    print()

    results = []  # (method, param, r2, rmse, rho)

    # ── PLS ──────────────────────────────────────────────────────────────────
    comp_grid = [k for k in PLS_N_COMPONENTS if k < n_samples - 1]
    print(f"Running PLS LOO-CV (n_components sweep {comp_grid}) ...")
    pls_scores = {}
    for k in comp_grid:
        y_pred = loo_cv(X, y, lambda k=k: PLSRegression(n_components=k))
        pls_scores[k] = metrics(y, y_pred)
        r2, rmse, rho = pls_scores[k]
        print(f"  n_components={k:2d}  R²={r2:+.3f}  RMSE={rmse:.4f}  ρ={rho:+.3f}")
    best_k = max(pls_scores, key=lambda k: pls_scores[k][0])
    r2, rmse, rho = pls_scores[best_k]
    print(f"  -> optimal n_components={best_k}")
    results.append(("PLS", f"n_components={best_k}", r2, rmse, rho))

    # ── Ridge ────────────────────────────────────────────────────────────────
    # cv=None uses the GCV formula: exact LOO within the training fold, no refits.
    print("Running Ridge (RidgeCV, GCV) LOO-CV ...")
    y_pred = loo_cv(X, y, lambda: RidgeCV(alphas=RIDGE_GRID, cv=None))
    r2, rmse, rho = metrics(y, y_pred)
    print(f"  R²={r2:+.3f}  RMSE={rmse:.4f}  ρ={rho:+.3f}")
    results.append(("Ridge", "GCV", r2, rmse, rho))

    # ── LASSO ────────────────────────────────────────────────────────────────
    print("Running LASSO (LassoCV, cv=3) LOO-CV ...")
    y_pred = loo_cv(X, y, lambda: LassoCV(cv=3, max_iter=5000, tol=0.01, alphas=N_ALPHAS))
    r2, rmse, rho = metrics(y, y_pred)
    print(f"  R²={r2:+.3f}  RMSE={rmse:.4f}  ρ={rho:+.3f}")
    results.append(("LASSO", "cv=3", r2, rmse, rho))

    # ── Elastic net ──────────────────────────────────────────────────────────
    print(f"Running elastic net (ElasticNetCV, cv=3, l1_ratio={L1_RATIOS}) LOO-CV ...")
    y_pred = loo_cv(
        X, y,
        lambda: ElasticNetCV(cv=3, l1_ratio=L1_RATIOS, alphas=N_ALPHAS,
                              max_iter=5000, tol=0.01),
    )
    r2, rmse, rho = metrics(y, y_pred)
    print(f"  R²={r2:+.3f}  RMSE={rmse:.4f}  ρ={rho:+.3f}")
    results.append(("Elastic net", f"l1_ratio sweep {L1_RATIOS}", r2, rmse, rho))

    # ── Summary table ────────────────────────────────────────────────────────
    print()
    print("=" * 68)
    print(f"METHOD COMPARISON - LOO-CV, {args.condition} "
          f"({n_samples} samples x {n_features} features)")
    print("=" * 68)
    print(f"  {'Method':<14} {'param':<26} {'CV R²':>7}  {'RMSE':>7}  {'Spearman ρ':>10}")
    print("  " + "-" * 64)
    best_method = max(results, key=lambda r: r[2])[0]
    for method, param, r2, rmse, rho in results:
        marker = "  <-" if method == best_method else ""
        print(f"  {method:<14} {param:<26} {r2:>7.3f}  {rmse:>7.4f}  {rho:>+10.3f}{marker}")
    print("=" * 68)

    # ── Append to summary CSV ────────────────────────────────────────────────
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = pd.DataFrame([
        {
            "condition": args.condition,
            "method": method,
            "param": param,
            "n_samples": n_samples,
            "n_features": n_features,
            "cv_r2": r2,
            "rmse": rmse,
            "spearman_rho": rho,
            "date_run": date.today().isoformat(),
        }
        for method, param, r2, rmse, rho in results
    ])
    write_header = not SUMMARY_PATH.exists()
    rows.to_csv(SUMMARY_PATH, mode="a", header=write_header, index=False, columns=SUMMARY_COLS)
    print(f"\nAppended {len(rows)} rows to {SUMMARY_PATH.relative_to(REPO)}")


if __name__ == "__main__":
    main()
