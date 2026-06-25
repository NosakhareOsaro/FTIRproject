"""
run_regularised_regression.py

Ridge, LASSO, and elastic net on DGRP line-mean spectra with
leave-one-line-out cross-validation.  Regularisation strength is
selected automatically inside each training fold (correct nested CV —
the held-out test line never influences hyperparameter selection).

All CV uses LeaveOneOut over 108 DGRP lines.  StandardScaler is fitted
on the 107 training lines inside each fold.

Outputs written to results/DGRP/:
  regularised_coefficients_vs_wavenumber.pdf   (best model)
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import ElasticNetCV, LassoCV, RidgeCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))
from ftir_loader import load_ftir

OUT_DIR = REPO / "results" / "DGRP"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Results carried forward from run_compression_analysis.py and run_pls_analysis.py
# (LOO-CV on the same 108 DGRP line-mean spectra)
PRIOR_RESULTS = [
    ("PCA + Ridge  (α=1, 95% var)", "4 PCs",   0.553, 0.4927, +0.743),
    ("PLS  (n_components = 10)",    "10 comp",  0.623, 0.4524, +0.801),
]

L1_RATIOS  = [0.5, 0.9, 1.0]         # ElasticNet: L1/L2 mixing sweep
RIDGE_GRID = np.logspace(-3, 6, 100) # RidgeCV: α candidates
N_ALPHAS   = 30                       # LassoCV / ElasticNetCV: α grid size

plt.rcParams.update({
    "figure.dpi": 150,
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# ── Load and prepare data ─────────────────────────────────────────────────────

meta, spectra = load_ftir(REPO / "FTIR-data" / "DGRPFTIR.dat")
female_mask = meta["Sex"] == "F"
meta    = meta[female_mask].reset_index(drop=True)
spectra = spectra[female_mask].reset_index(drop=True)

wavenumbers = np.array(spectra.columns.tolist())   # [3900, 3898, …, 456]

spec_df = spectra.copy()
spec_df["Genot."] = meta["Genot."].values
X_line_df = spec_df.groupby("Genot.").mean()       # 108 × 1723

emmeans  = pd.read_csv(REPO / "Emmeans.csv")
line_data = (
    X_line_df.reset_index()
    .merge(emmeans[["DGRP", "emmean"]], left_on="Genot.", right_on="DGRP", how="inner")
    .sort_values("Genot.")
    .reset_index(drop=True)
)

X_line  = line_data[wavenumbers.tolist()].values   # 108 × 1723
y_line  = line_data["emmean"].values               # 108
n_lines = len(y_line)

print(f"Lines in analysis : {n_lines}  ({X_line.shape[1]} wavenumbers)")

# ── Helpers ───────────────────────────────────────────────────────────────────

loo = LeaveOneOut()


def metrics(y_true, y_pred):
    return (
        r2_score(y_true, y_pred),
        np.sqrt(mean_squared_error(y_true, y_pred)),
        spearmanr(y_true, y_pred).statistic,
    )


def loo_cv_reg(X, y, model_fn, extra_attrs=()):
    """
    LOO-CV with per-fold StandardScaling and hyperparameter tracking.

    model_fn must return a fresh fitted-CV estimator (RidgeCV, LassoCV,
    or ElasticNetCV) — hyperparameter selection happens inside each fold
    on the 107 training lines only.

    Returns y_pred (n,), alphas (n,), extras dict of tracked attributes.
    """
    n = len(y)
    y_pred = np.zeros(n)
    alphas = []
    extras = {a: [] for a in extra_attrs}
    for train_idx, test_idx in loo.split(X):
        sc = StandardScaler()
        X_tr_s = sc.fit_transform(X[train_idx])
        X_te_s = sc.transform(X[test_idx])
        mdl = model_fn()
        mdl.fit(X_tr_s, y[train_idx])
        y_pred[test_idx] = np.asarray(mdl.predict(X_te_s)).ravel()
        alphas.append(mdl.alpha_)
        for a in extra_attrs:
            extras[a].append(getattr(mdl, a))
    return y_pred, np.array(alphas), {k: np.array(v) for k, v in extras.items()}


# ═══════════════════════════════════════════════════════════════════════════════
# Part A — LOO-CV for Ridge, LASSO, elastic net
# ═══════════════════════════════════════════════════════════════════════════════

# ── Ridge ─────────────────────────────────────────────────────────────────────
# cv=None uses the GCV formula — exact LOO within the training fold, no refits.

print("\nRunning Ridge (RidgeCV, GCV) LOO-CV …")
y_pred_ridge, alphas_ridge, _ = loo_cv_reg(
    X_line, y_line,
    lambda: RidgeCV(alphas=RIDGE_GRID, cv=None),
)
r2_ridge, rmse_ridge, rho_ridge = metrics(y_line, y_pred_ridge)
alpha_ridge_med = np.median(alphas_ridge)
print(f"  R²={r2_ridge:+.3f}  RMSE={rmse_ridge:.4f}  ρ={rho_ridge:+.3f}  "
      f"median α={alpha_ridge_med:.3g}")

# ── LASSO ─────────────────────────────────────────────────────────────────────

print("Running LASSO (LassoCV, cv=3) LOO-CV …")
y_pred_lasso, alphas_lasso, _ = loo_cv_reg(
    X_line, y_line,
    lambda: LassoCV(cv=3, max_iter=5000, tol=0.01, alphas=N_ALPHAS),
)
r2_lasso, rmse_lasso, rho_lasso = metrics(y_line, y_pred_lasso)
alpha_lasso_med = np.median(alphas_lasso)
print(f"  R²={r2_lasso:+.3f}  RMSE={rmse_lasso:.4f}  ρ={rho_lasso:+.3f}  "
      f"median α={alpha_lasso_med:.4g}")

# ── Elastic net ───────────────────────────────────────────────────────────────
# Slowest step: sweeps α × l1_ratio with 5-fold inner CV across 108 outer folds.

print(f"Running elastic net (ElasticNetCV, cv=3, l1_ratio={L1_RATIOS}) LOO-CV …")
y_pred_enet, alphas_enet, extras_enet = loo_cv_reg(
    X_line, y_line,
    lambda: ElasticNetCV(
        cv=3, l1_ratio=L1_RATIOS,
        max_iter=5000, tol=0.01, alphas=N_ALPHAS,
    ),
    extra_attrs=("l1_ratio_",),
)
r2_enet, rmse_enet, rho_enet = metrics(y_line, y_pred_enet)
alpha_enet_med   = np.median(alphas_enet)
l1_ratio_enet_med = np.median(extras_enet["l1_ratio_"])
print(f"  R²={r2_enet:+.3f}  RMSE={rmse_enet:.4f}  ρ={rho_enet:+.3f}  "
      f"median α={alpha_enet_med:.4g}  median ℓ₁={l1_ratio_enet_med:.2f}")

# ═══════════════════════════════════════════════════════════════════════════════
# Part B — Coefficient plot for the best regularised model
# ═══════════════════════════════════════════════════════════════════════════════

new_results = {
    "Ridge":       (r2_ridge, rmse_ridge, rho_ridge),
    "LASSO":       (r2_lasso, rmse_lasso, rho_lasso),
    "Elastic net": (r2_enet,  rmse_enet,  rho_enet),
}
best_name = max(new_results, key=lambda k: new_results[k][0])
print(f"\nBest regularised model: {best_name} (CV R² = {new_results[best_name][0]:.3f})")

# Fit best model on all 108 lines; CV variant selects α from the full dataset.
scaler_full  = StandardScaler()
X_scaled_all = scaler_full.fit_transform(X_line)

if best_name == "Ridge":
    final_model = RidgeCV(alphas=RIDGE_GRID, cv=None)
elif best_name == "LASSO":
    final_model = LassoCV(cv=3, max_iter=5000, tol=0.01, alphas=N_ALPHAS)
else:
    final_model = ElasticNetCV(
        cv=3, l1_ratio=L1_RATIOS,
        max_iter=5000, tol=0.01, alphas=N_ALPHAS,
    )

final_model.fit(X_scaled_all, y_line)
coef      = final_model.coef_            # shape (1723,)
n_nonzero = int(np.sum(coef != 0))
alpha_final = final_model.alpha_

subtitle_parts = [f"α = {alpha_final:.4g}", f"n_nonzero = {n_nonzero} / {len(wavenumbers)}"]
if hasattr(final_model, "l1_ratio_"):
    subtitle_parts.insert(1, f"ℓ₁ = {final_model.l1_ratio_:.2f}")

fig, ax = plt.subplots(figsize=(9, 3.5))
ax.plot(wavenumbers, coef, lw=0.9, color="#2166ac")
ax.axhline(0, color="#888888", lw=0.7, ls="--")
ax.set_xlabel("Wavenumber (cm⁻¹)")
ax.set_ylabel("Coefficient (standardised units)")
ax.set_title(
    f"{best_name} coefficient vector — spectral drivers of starvation-resistance prediction\n"
    + "  ".join(subtitle_parts)
)
fig.tight_layout()
out_coef = OUT_DIR / "regularised_coefficients_vs_wavenumber.pdf"
fig.savefig(out_coef)
plt.close(fig)
print(f"Saved : {out_coef.relative_to(REPO)}")

# ═══════════════════════════════════════════════════════════════════════════════
# Full comparison table
# ═══════════════════════════════════════════════════════════════════════════════

ridge_param = f"α={alpha_ridge_med:.3g}"
lasso_param = f"α={alpha_lasso_med:.4g}"
enet_param  = f"α={alpha_enet_med:.4g}, ℓ₁={l1_ratio_enet_med:.2f}"

all_rows = PRIOR_RESULTS + [
    ("Ridge",       ridge_param, r2_ridge, rmse_ridge, rho_ridge),
    ("LASSO",       lasso_param, r2_lasso, rmse_lasso, rho_lasso),
    ("Elastic net", enet_param,  r2_enet,  rmse_enet,  rho_enet),
]

print()
print("=" * 72)
print("METHOD COMPARISON — LOO-CV on 108 DGRP line-mean spectra")
print("=" * 72)
print(f"  {'Method':<30}  {'param':<20}  {'CV R²':>6}  {'RMSE':>7}  {'Spearman ρ':>10}")
print("  " + "-" * 68)
for method, param, r2, rmse, rho in all_rows:
    marker = "  ←" if method == best_name else ""
    print(f"  {method:<30}  {param:<20}  {r2:>6.3f}  {rmse:>7.4f}  {rho:>+10.3f}{marker}")
print("=" * 72)
print()
print("Note: PCA+Ridge and PLS rows carried forward from prior scripts.")
print("      Ridge/LASSO/ElasticNet α selected by inner CV (cv=3) inside")
print("      each outer LOO fold — test line never influences α selection.")
print("      LASSO/ElasticNet use tol=0.01 (sufficient for α selection).")
print()
print("Plots saved:")
print(f"  {out_coef.relative_to(REPO)}")
