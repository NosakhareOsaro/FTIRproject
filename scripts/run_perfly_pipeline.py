"""
run_perfly_pipeline.py

Per-fly evaluation of four regression methods with line-stratified
GroupKFold cross-validation.

Training unit: individual fly spectra (~1,772 females, 108 DGRP lines).
Target: each fly's line-level starvation-resistance EMMean — repeated for
every fly in the same line because the FTIR and starvation assays are both
destructive (fly-level pairing is impossible).

CV design (critical):
  Outer: GroupKFold(n_splits=10) with DGRP line as the group.
  No fly from a given line ever appears in both train and test in the same
  fold.  Violating this would allow the line-level target to leak across
  folds and inflate all metrics — it is the single most important
  correctness requirement in this pipeline.

Evaluation:
  After predicting all test-fold flies, average predictions within each
  test-fold DGRP line → one predicted value per line.  Compute R², RMSE,
  Spearman ρ, Pearson r against that line's true EMMean.  This is
  line-level evaluation on per-fly-trained models.

  Per-fly R² is also reported for the comparison table but is expected to
  be much lower than line R²: averaging cancels within-line spectral noise
  that the model cannot learn.

Methods:
  PLS        — PLSRegression; n_components ∈ {1,2,3,5,10} selected by
               inner GroupKFold(n_splits=5) on the outer training fold.
               GroupKFold used for the inner sweep to prevent the same
               line appearing on both sides of an inner split.
  Ridge      — RidgeCV, GCV α selection (no explicit inner loop needed —
               GCV solves the full α path analytically on the training fold).
  LASSO      — LassoCV, cv=3, n_alphas=30.  Random 3-fold for α selection
               only; test lines are already excluded by the outer fold.
  ElasticNet — ElasticNetCV, cv=3, l1_ratio=[0.5,0.9,1.0], n_alphas=30.

StandardScaler is fitted on the outer-fold training flies only inside
each outer fold.

Outputs written to results/DGRP/:
  perfly_metrics.csv
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.cross_decomposition import PLSRegression
from sklearn.linear_model import ElasticNetCV, LassoCV, RidgeCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, cross_val_score
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))
from ftir_loader import load_ftir

OUT_DIR = REPO / "results" / "DGRP"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_OUTER     = 10
N_INNER_PLS = 5
N_COMP_LIST = [1, 2, 3, 5, 10]
L1_RATIOS   = [0.5, 0.9, 1.0]
N_ALPHAS    = 30
RIDGE_GRID  = np.logspace(-3, 6, 100)

# Line-mean LOO-CV results carried forward from run_pls_analysis.py and
# run_regularised_regression.py (same 108-line dataset, different CV scheme).
LINE_MEAN_R2 = {
    "PLS":        0.623,
    "Ridge":      0.635,
    "LASSO":      0.669,
    "ElasticNet": 0.673,
}

# ── Load and prepare data ─────────────────────────────────────────────────────

meta, spectra = load_ftir(REPO / "FTIR-data" / "DGRPFTIR.dat")
female_mask = meta["Sex"] == "F"
meta    = meta[female_mask].reset_index(drop=True)
spectra = spectra[female_mask].reset_index(drop=True)

emmeans = pd.read_csv(REPO / "Emmeans.csv")
em_map  = emmeans.set_index("DGRP")["emmean"].to_dict()

# Attach EMMean to each fly; drop flies from lines absent in Emmeans.csv
fly_df          = meta[["Genot."]].copy()
fly_df["emmean"] = fly_df["Genot."].map(em_map)
keep    = fly_df["emmean"].notna()
fly_df  = fly_df[keep].reset_index(drop=True)
spectra = spectra[keep].reset_index(drop=True)

X        = spectra.values.astype(float)   # (n_flies, 1723)
y        = fly_df["emmean"].values         # repeated EMMean per fly
line_ids = fly_df["Genot."].values         # DGRP line string per fly

# Integer-encode groups for GroupKFold
line_labels, groups = np.unique(line_ids, return_inverse=True)
n_lines = len(line_labels)
n_flies = len(y)

print(f"Flies in analysis : {n_flies}  ({X.shape[1]} wavenumbers)")
print(f"DGRP lines        : {n_lines}")
print()

# ── CV helpers ────────────────────────────────────────────────────────────────

outer_cv    = GroupKFold(n_splits=N_OUTER)
inner_cv_pls = GroupKFold(n_splits=N_INNER_PLS)


def run_outer_cv(method_name, model_fn):
    """
    Outer GroupKFold loop.

    model_fn(X_tr_scaled, y_tr, groups_tr) → fitted estimator with .predict().
    Returns DataFrame with columns [line_id, fly_pred, emmean].
    """
    records = []
    for fold, (train_idx, test_idx) in enumerate(outer_cv.split(X, y, groups)):
        sc   = StandardScaler()
        X_tr = sc.fit_transform(X[train_idx])
        X_te = sc.transform(X[test_idx])

        mdl   = model_fn(X_tr, y[train_idx], groups[train_idx])
        preds = np.asarray(mdl.predict(X_te)).ravel()

        for i, pred in zip(test_idx, preds):
            records.append({
                "line_id":  line_ids[i],
                "fly_pred": pred,
                "emmean":   y[i],
            })

        n_test_lines = len(np.unique(groups[test_idx]))
        print(f"  fold {fold + 1:2d}/{N_OUTER}  "
              f"test lines={n_test_lines:3d}  "
              f"train flies={len(train_idx):5d}  "
              f"test flies={len(test_idx):4d}")

    return pd.DataFrame(records)


def line_aggregate(pred_df):
    """Average per-fly predictions to line level and compute all metrics."""
    line_df = (
        pred_df
        .groupby("line_id")
        .agg(line_pred=("fly_pred", "mean"), emmean=("emmean", "first"))
        .reset_index()
    )
    line_r2   = r2_score(line_df["emmean"], line_df["line_pred"])
    line_rmse = np.sqrt(mean_squared_error(line_df["emmean"], line_df["line_pred"]))
    line_rho  = spearmanr(line_df["emmean"], line_df["line_pred"]).statistic
    line_r    = pearsonr(line_df["emmean"], line_df["line_pred"])[0]
    perfly_r2 = r2_score(pred_df["emmean"], pred_df["fly_pred"])
    return {
        "line_r2":    line_r2,
        "line_rmse":  line_rmse,
        "line_rho":   line_rho,
        "line_r":     line_r,
        "perfly_r2":  perfly_r2,
        "n_lines":    len(line_df),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# PLS — n_components selected by inner GroupKFold(5)
# ═══════════════════════════════════════════════════════════════════════════════

print("Running PLS  (outer GroupKFold=10, inner GroupKFold=5 for n_components) …")


def pls_model_fn(X_tr, y_tr, groups_tr):
    best_k, best_score = N_COMP_LIST[0], -np.inf
    for k in N_COMP_LIST:
        scores = cross_val_score(
            PLSRegression(n_components=k),
            X_tr, y_tr,
            cv=inner_cv_pls,
            groups=groups_tr,
            scoring="r2",
        )
        if scores.mean() > best_score:
            best_score = scores.mean()
            best_k = k
    mdl = PLSRegression(n_components=best_k)
    mdl.fit(X_tr, y_tr)
    return mdl


pls_preds   = run_outer_cv("PLS", pls_model_fn)
pls_m       = line_aggregate(pls_preds)
print(f"  → line R²={pls_m['line_r2']:+.3f}  RMSE={pls_m['line_rmse']:.4f}  "
      f"ρ={pls_m['line_rho']:+.3f}  r={pls_m['line_r']:+.3f}  "
      f"per-fly R²={pls_m['perfly_r2']:+.3f}\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Ridge — RidgeCV with GCV (no explicit inner loop)
# ═══════════════════════════════════════════════════════════════════════════════

print("Running Ridge  (RidgeCV, GCV α selection) …")


def ridge_model_fn(X_tr, y_tr, groups_tr):
    mdl = RidgeCV(alphas=RIDGE_GRID, cv=None)
    mdl.fit(X_tr, y_tr)
    return mdl


ridge_preds = run_outer_cv("Ridge", ridge_model_fn)
ridge_m     = line_aggregate(ridge_preds)
print(f"  → line R²={ridge_m['line_r2']:+.3f}  RMSE={ridge_m['line_rmse']:.4f}  "
      f"ρ={ridge_m['line_rho']:+.3f}  r={ridge_m['line_r']:+.3f}  "
      f"per-fly R²={ridge_m['perfly_r2']:+.3f}\n")

# ═══════════════════════════════════════════════════════════════════════════════
# LASSO
# ═══════════════════════════════════════════════════════════════════════════════

print("Running LASSO  (LassoCV, cv=3) …")


def lasso_model_fn(X_tr, y_tr, groups_tr):
    mdl = LassoCV(cv=3, alphas=N_ALPHAS, max_iter=5000, tol=0.01)
    mdl.fit(X_tr, y_tr)
    return mdl


lasso_preds = run_outer_cv("LASSO", lasso_model_fn)
lasso_m     = line_aggregate(lasso_preds)
print(f"  → line R²={lasso_m['line_r2']:+.3f}  RMSE={lasso_m['line_rmse']:.4f}  "
      f"ρ={lasso_m['line_rho']:+.3f}  r={lasso_m['line_r']:+.3f}  "
      f"per-fly R²={lasso_m['perfly_r2']:+.3f}\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Elastic net
# ═══════════════════════════════════════════════════════════════════════════════

print(f"Running elastic net  (ElasticNetCV, cv=3, l1_ratio={L1_RATIOS}) …")


def enet_model_fn(X_tr, y_tr, groups_tr):
    mdl = ElasticNetCV(
        cv=3, l1_ratio=L1_RATIOS,
        alphas=N_ALPHAS, max_iter=5000, tol=0.01,
    )
    mdl.fit(X_tr, y_tr)
    return mdl


enet_preds = run_outer_cv("ElasticNet", enet_model_fn)
enet_m     = line_aggregate(enet_preds)
print(f"  → line R²={enet_m['line_r2']:+.3f}  RMSE={enet_m['line_rmse']:.4f}  "
      f"ρ={enet_m['line_rho']:+.3f}  r={enet_m['line_r']:+.3f}  "
      f"per-fly R²={enet_m['perfly_r2']:+.3f}\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Summary tables
# ═══════════════════════════════════════════════════════════════════════════════

all_results = [
    ("PLS",        pls_m),
    ("Ridge",      ridge_m),
    ("LASSO",      lasso_m),
    ("ElasticNet", enet_m),
]

print("=" * 66)
print("COMPARISON TABLE  —  line-mean LOO-CV  vs  per-fly GroupKFold(10)")
print("=" * 66)
print(f"  {'Method':<14}  {'Line-mean R²':>12}  {'Per-fly R²':>10}  {'Change':>8}")
print("  " + "-" * 50)
for name, m in all_results:
    lm  = LINE_MEAN_R2[name]
    pf  = m["perfly_r2"]
    ch  = pf - lm
    print(f"  {name:<14}  {lm:>12.3f}  {pf:>10.3f}  {ch:>+8.3f}")

print()
print("=" * 66)
print("PER-FLY PIPELINE  —  full metrics (line-level evaluation after averaging)")
print("=" * 66)
print(f"  {'Method':<14}  {'Line R²':>8}  {'RMSE':>7}  "
      f"{'Spearman ρ':>10}  {'Pearson r':>9}  {'n lines':>7}")
print("  " + "-" * 62)
for name, m in all_results:
    print(f"  {name:<14}  {m['line_r2']:>8.3f}  {m['line_rmse']:>7.4f}  "
          f"{m['line_rho']:>+10.3f}  {m['line_r']:>+9.3f}  {m['n_lines']:>7d}")
print("=" * 66)
print()
print("Note: 'per-fly R²' and 'line R²' differ because per-fly R² is computed")
print("      against individual-fly targets before line-averaging.  'Line R²'")
print("      is the dissertation metric: line-mean predictions vs EMMeans.")
print("      Line-mean LOO-CV rows use 108 line-mean spectra as training units.")

# ── Save CSV ──────────────────────────────────────────────────────────────────

rows = []
for name, m in all_results:
    rows.append({
        "method":         name,
        "linemean_r2":    LINE_MEAN_R2[name],
        "perfly_r2":      m["perfly_r2"],
        "line_r2":        m["line_r2"],
        "line_rmse":      m["line_rmse"],
        "line_spearman":  m["line_rho"],
        "line_pearson":   m["line_r"],
        "n_lines":        m["n_lines"],
    })

out_csv = OUT_DIR / "perfly_metrics.csv"
pd.DataFrame(rows).to_csv(out_csv, index=False, float_format="%.4f")
print(f"\nSaved : {out_csv.relative_to(REPO)}")
