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

Each condition is run through three analyses, to answer whether the
Sex/Age demographic confound (previously checked only as a standalone
trivial baseline) changes the lipidomics-BMI result when incorporated
directly into the same modelling pipeline:

  1. lipidomics_only         - the original analysis: BMI ~ lipidomics features.
  2. covariates_appended     - Sex (0/1) and Age appended as two extra columns
                                to the feature matrix before running the same
                                four methods, same LOO-CV/per-fold scaling.
  3. residualized            - a simple BMI ~ Sex + Age regression is fit
                                *inside each LOO training fold* (never on the
                                full data), and the four methods predict the
                                residual BMI (actual BMI minus that fold's
                                Sex+Age prediction) from the lipidomics
                                features alone. Because the covariate model is
                                refit per fold, the "true" residual for the
                                held-out point also depends on that fold's fit
                                - loo_cv_residualized returns both the
                                predicted and the actual (per-fold) residual.

Usage:
  .venv/bin/python scripts/run_eatris_lipidomics_bmi.py \\
      phenotype-data/EATRIS_Lipidomics_positive.tsv --condition "Positive mode"

  .venv/bin/python scripts/run_eatris_lipidomics_bmi.py \\
      phenotype-data/EATRIS_Lipidomics_negative.tsv --condition "Negative mode"

  .venv/bin/python scripts/run_eatris_lipidomics_bmi.py \\
      phenotype-data/EATRIS_Lipidomics_combined.tsv --condition "Combined"

Treated as three separate test conditions, the same way
run_dgrpool_phenotype.py is invoked once per Unckless diet condition:
each run is independent and appends its own rows to the summary CSV
(12 rows per run: 3 analyses x 4 methods).

Sex/Age come from phenotype-data/EATRIS_covariates.tsv (Sex encoded
FEMALE=0, MALE=1; only these two categories are present in the cohort).
"""

import argparse
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.cross_decomposition import PLSRegression
from sklearn.linear_model import ElasticNetCV, LassoCV, LinearRegression, RidgeCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).parent.parent
BMI_PATH = REPO / "phenotype-data" / "EATRIS_BMI.tsv"
COVARIATES_PATH = REPO / "phenotype-data" / "EATRIS_covariates.tsv"
OUT_DIR = REPO / "results" / "EATRIS"
SUMMARY_PATH = OUT_DIR / "lipidomics_bmi_summary.csv"
SUMMARY_COLS = [
    "condition", "analysis", "method", "param", "n_samples", "n_features",
    "cv_r2", "rmse", "spearman_rho", "date_run",
]
SEX_MAP = {"FEMALE": 0, "MALE": 1}

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


def loo_cv_residualized(X_lipid, X_cov, y, model_fn):
    """LOO-CV predicting covariate-adjusted BMI residuals from lipidomics features.

    Within each fold: fit plain LinearRegression BMI ~ Sex + Age on the
    training fold only, subtract its prediction from BMI on both the training
    and held-out rows to get residual BMI, then fit model_fn (with its own
    per-fold StandardScaler, same discipline as loo_cv) on the lipidomics
    features to predict the training residuals, and predict the held-out
    residual. The covariate model is refit fresh every fold, so the "actual"
    residual for the held-out point is fold-specific too - returned alongside
    the prediction rather than reusing a single full-data residual vector.
    """
    resid_pred = np.zeros(len(y))
    resid_actual = np.zeros(len(y))
    for train_idx, test_idx in loo.split(X_lipid):
        cov_mdl = LinearRegression()
        cov_mdl.fit(X_cov[train_idx], y[train_idx])
        resid_train = y[train_idx] - cov_mdl.predict(X_cov[train_idx])
        resid_actual[test_idx] = y[test_idx] - cov_mdl.predict(X_cov[test_idx])

        sc = StandardScaler()
        X_tr = sc.fit_transform(X_lipid[train_idx])
        X_te = sc.transform(X_lipid[test_idx])
        mdl = model_fn()
        mdl.fit(X_tr, resid_train)
        resid_pred[test_idx] = np.asarray(mdl.predict(X_te)).ravel()
    return resid_pred, resid_actual


def run_suite(loo_runner, n_samples, n_features, label):
    """Run PLS/Ridge/LASSO/elastic net LOO-CV via loo_runner(model_fn) -> (y_pred, y_true).

    Kept generic over loo_runner so the same four-method sweep serves the
    plain (X, y), covariates-appended, and residualized analyses without
    duplicating the method-by-method code three times.
    """
    results = []  # (method, param, r2, rmse, rho)

    comp_grid = [k for k in PLS_N_COMPONENTS if k < n_samples - 1]
    print(f"  [{label}] Running PLS LOO-CV (n_components sweep {comp_grid}) ...")
    pls_scores = {}
    for k in comp_grid:
        y_pred, y_true = loo_runner(lambda k=k: PLSRegression(n_components=k))
        pls_scores[k] = metrics(y_true, y_pred)
        r2, rmse, rho = pls_scores[k]
        print(f"    n_components={k:2d}  R²={r2:+.3f}  RMSE={rmse:.4f}  ρ={rho:+.3f}")
    best_k = max(pls_scores, key=lambda k: pls_scores[k][0])
    r2, rmse, rho = pls_scores[best_k]
    print(f"    -> optimal n_components={best_k}")
    results.append(("PLS", f"n_components={best_k}", r2, rmse, rho))

    print(f"  [{label}] Running Ridge (RidgeCV, GCV) LOO-CV ...")
    y_pred, y_true = loo_runner(lambda: RidgeCV(alphas=RIDGE_GRID, cv=None))
    r2, rmse, rho = metrics(y_true, y_pred)
    print(f"    R²={r2:+.3f}  RMSE={rmse:.4f}  ρ={rho:+.3f}")
    results.append(("Ridge", "GCV", r2, rmse, rho))

    print(f"  [{label}] Running LASSO (LassoCV, cv=3) LOO-CV ...")
    y_pred, y_true = loo_runner(lambda: LassoCV(cv=3, max_iter=5000, tol=0.01, alphas=N_ALPHAS))
    r2, rmse, rho = metrics(y_true, y_pred)
    print(f"    R²={r2:+.3f}  RMSE={rmse:.4f}  ρ={rho:+.3f}")
    results.append(("LASSO", "cv=3", r2, rmse, rho))

    print(f"  [{label}] Running elastic net (ElasticNetCV, cv=3, l1_ratio={L1_RATIOS}) LOO-CV ...")
    y_pred, y_true = loo_runner(
        lambda: ElasticNetCV(cv=3, l1_ratio=L1_RATIOS, alphas=N_ALPHAS,
                              max_iter=5000, tol=0.01),
    )
    r2, rmse, rho = metrics(y_true, y_pred)
    print(f"    R²={r2:+.3f}  RMSE={rmse:.4f}  ρ={rho:+.3f}")
    results.append(("Elastic net", f"l1_ratio sweep {L1_RATIOS}", r2, rmse, rho))

    return results


def print_summary_table(results, condition, label, n_samples, n_features):
    print()
    print("=" * 68)
    print(f"METHOD COMPARISON - LOO-CV, {condition} [{label}] "
          f"({n_samples} samples x {n_features} features)")
    print("=" * 68)
    print(f"  {'Method':<14} {'param':<26} {'CV R²':>7}  {'RMSE':>7}  {'Spearman ρ':>10}")
    print("  " + "-" * 64)
    best_method = max(results, key=lambda r: r[2])[0]
    for method, param, r2, rmse, rho in results:
        marker = "  <-" if method == best_method else ""
        print(f"  {method:<14} {param:<26} {r2:>7.3f}  {rmse:>7.4f}  {rho:>+10.3f}{marker}")
    print("=" * 68)


def build_rows(results, condition, analysis, n_samples, n_features):
    return pd.DataFrame([
        {
            "condition": condition,
            "analysis": analysis,
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


def ensure_summary_schema():
    """Backfill an 'analysis' column into a pre-existing summary CSV.

    Rows appended before this script supported multiple analyses are all
    the original lipidomics-only run; tagged as such so old and new rows
    stay comparable in one file rather than silently misaligning columns.
    """
    if not SUMMARY_PATH.exists():
        return
    existing = pd.read_csv(SUMMARY_PATH)
    if "analysis" in existing.columns:
        return
    existing.insert(existing.columns.get_loc("method"), "analysis", "lipidomics_only")
    existing.to_csv(SUMMARY_PATH, index=False, columns=SUMMARY_COLS)
    print(f"Migrated {SUMMARY_PATH.relative_to(REPO)} to add an 'analysis' column "
          f"(existing rows backfilled as 'lipidomics_only').")


def main():
    parser = argparse.ArgumentParser(
        description="LOO-CV method comparison (PLS, Ridge, LASSO, elastic net) "
        "predicting BMI from an EATRIS-Plus lipidomics feature matrix, run as "
        "lipidomics-only, with Sex/Age covariates appended, and against "
        "Sex/Age-adjusted BMI residuals."
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
    cov_df = pd.read_csv(COVARIATES_PATH, sep="\t")[["sample_id", "Sex", "Age"]]

    merged = X_df.merge(bmi_df, on="sample_id", how="inner")
    n_dropped = len(X_df) - len(merged)
    if n_dropped:
        print(f"Note: {n_dropped} sample(s) in {args.features_tsv.name} had no "
              f"matching BMI value and were dropped.")

    merged = merged.merge(cov_df, on="sample_id", how="inner")
    n_dropped_cov = (len(X_df) - n_dropped) - len(merged)
    if n_dropped_cov:
        print(f"Note: {n_dropped_cov} sample(s) had no matching Sex/Age covariate "
              f"row and were dropped.")

    unknown_sex = sorted(set(merged["Sex"].unique()) - set(SEX_MAP))
    if unknown_sex:
        raise ValueError(f"Unrecognized Sex value(s) in {COVARIATES_PATH.name}: "
                          f"{unknown_sex} (expected only {sorted(SEX_MAP)})")

    feature_cols = [c for c in X_df.columns if c != "sample_id"]
    X = merged[feature_cols].values
    y = merged["BMI"].values.astype(float)
    sex_enc = merged["Sex"].map(SEX_MAP).values.astype(float)
    age = merged["Age"].values.astype(float)
    X_cov = np.column_stack([sex_enc, age])
    n_samples, n_features = X.shape

    print("=" * 68)
    print(f"Condition : {args.condition}")
    print(f"File      : {args.features_tsv}")
    print(f"Samples   : {n_samples}   Features: {n_features}")
    print(f"BMI range : {y.min():.1f} - {y.max():.1f}  "
          f"(mean {y.mean():.2f}, SD {y.std():.2f})")
    print(f"Sex       : {int(sex_enc.sum())} male, {int(len(sex_enc) - sex_enc.sum())} female "
          f"(encoded {SEX_MAP})")
    print(f"Age range : {age.min():.0f} - {age.max():.0f}  (mean {age.mean():.1f})")
    print("=" * 68)

    all_rows = []

    # ── Analysis 1/3: lipidomics only (original analysis) ──────────────────────
    print("\n--- Analysis 1/3: lipidomics only ---")
    results_lip = run_suite(
        lambda model_fn: (loo_cv(X, y, model_fn), y),
        n_samples, n_features, "lipidomics only",
    )
    print_summary_table(results_lip, args.condition, "lipidomics only", n_samples, n_features)
    all_rows.append(build_rows(results_lip, args.condition, "lipidomics_only", n_samples, n_features))

    # ── Analysis 2/3: Sex + Age appended to the feature matrix ─────────────────
    print("\n--- Analysis 2/3: lipidomics + Sex/Age covariates appended ---")
    X_aug = np.column_stack([X, X_cov])
    n_features_aug = X_aug.shape[1]
    results_cov = run_suite(
        lambda model_fn: (loo_cv(X_aug, y, model_fn), y),
        n_samples, n_features_aug, "covariates appended",
    )
    print_summary_table(results_cov, args.condition, "covariates appended", n_samples, n_features_aug)
    all_rows.append(build_rows(results_cov, args.condition, "covariates_appended", n_samples, n_features_aug))

    # ── Analysis 3/3: lipidomics predicting Sex/Age-adjusted BMI residuals ─────
    print("\n--- Analysis 3/3: lipidomics predicting Sex/Age-adjusted BMI residuals ---")
    results_resid = run_suite(
        lambda model_fn: loo_cv_residualized(X, X_cov, y, model_fn),
        n_samples, n_features, "covariate-adjusted residuals",
    )
    print_summary_table(results_resid, args.condition, "covariate-adjusted residuals", n_samples, n_features)
    all_rows.append(build_rows(results_resid, args.condition, "residualized", n_samples, n_features))

    # ── Cross-analysis comparison ───────────────────────────────────────────────
    print()
    print("=" * 68)
    print(f"CROSS-ANALYSIS COMPARISON (CV R²) - {args.condition}")
    print("=" * 68)
    print(f"  {'Method':<14} {'Lipidomics only':>16} {'+ Sex/Age':>10} {'Residualized':>13}")
    for (m, _, r2_lip, *_), (_, _, r2_cov, *_), (_, _, r2_res, *_) in zip(
        results_lip, results_cov, results_resid
    ):
        print(f"  {m:<14} {r2_lip:>16.3f} {r2_cov:>10.3f} {r2_res:>13.3f}")
    print("=" * 68)

    # ── Append to summary CSV ────────────────────────────────────────────────
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ensure_summary_schema()
    rows = pd.concat(all_rows, ignore_index=True)
    write_header = not SUMMARY_PATH.exists()
    rows.to_csv(SUMMARY_PATH, mode="a", header=write_header, index=False, columns=SUMMARY_COLS)
    print(f"\nAppended {len(rows)} rows to {SUMMARY_PATH.relative_to(REPO)}")


if __name__ == "__main__":
    main()
