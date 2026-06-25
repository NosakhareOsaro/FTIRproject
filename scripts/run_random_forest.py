"""
run_random_forest.py

Random forest regression on DGRP line-mean spectra with leave-one-line-out
cross-validation.  Hyperparameters (n_estimators, max_features) are selected
by inner 3-fold GridSearchCV on the 107 training lines inside each outer LOO
fold — the held-out test line never influences hyperparameter selection.

n_jobs=1 throughout: parallel RF on M1 causes memory pressure.

Outputs written to results/DGRP/:
  rf_feature_importance_vs_wavenumber.pdf
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))
from ftir_loader import load_ftir

OUT_DIR = REPO / "results" / "DGRP"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PARAM_GRID = {
    "n_estimators": [100, 300],
    "max_features": ["sqrt", "log2", 0.1, 0.3],
}
RANDOM_STATE = 0

# Prior line-mean LOO-CV results (hardcoded from run_compression_analysis.py,
# run_pls_analysis.py, and run_regularised_regression.py outputs)
PRIOR_RESULTS = [
    ("PCA + Ridge  (4 PCs)",  0.553, 0.4927, +0.743),
    ("PLS  (10 comp)",        0.623, 0.4524, +0.801),
    ("Ridge",                 0.635, 0.4451, +0.809),
    ("LASSO",                 0.669, 0.4266, +0.813),
    ("Elastic net",           0.673, 0.4244, +0.816),
]

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
print(f"Param grid        : {PARAM_GRID}")
print(f"Inner CV folds    : 3  (GridSearchCV, scoring=r2)")
print(f"Outer CV          : LeaveOneOut ({n_lines} folds)")
print()

# ── LOO-CV ────────────────────────────────────────────────────────────────────

loo     = LeaveOneOut()
y_pred  = np.zeros(n_lines)
best_params_log = []

print("Running Random Forest LOO-CV …")
for i, (train_idx, test_idx) in enumerate(loo.split(X_line)):
    sc      = StandardScaler()
    X_tr    = sc.fit_transform(X_line[train_idx])
    X_te    = sc.transform(X_line[test_idx])

    gs = GridSearchCV(
        RandomForestRegressor(random_state=RANDOM_STATE, n_jobs=1),
        PARAM_GRID,
        cv=3,
        scoring="r2",
        n_jobs=1,
    )
    gs.fit(X_tr, y_line[train_idx])
    y_pred[test_idx] = gs.predict(X_te)
    best_params_log.append(gs.best_params_)

    if (i + 1) % 10 == 0 or i == 0:
        print(f"  fold {i + 1:3d}/{n_lines}  best={gs.best_params_}  "
              f"inner R²={gs.best_score_:+.3f}")

r2_rf   = r2_score(y_line, y_pred)
rmse_rf = np.sqrt(mean_squared_error(y_line, y_pred))
rho_rf  = spearmanr(y_line, y_pred).statistic

print(f"\nRandom forest LOO-CV:  R²={r2_rf:+.3f}  RMSE={rmse_rf:.4f}  ρ={rho_rf:+.3f}")

# Most-selected hyperparameter combination
params_df   = pd.DataFrame(best_params_log)
modal_n_est = params_df["n_estimators"].mode()[0]
modal_mf    = params_df["max_features"].mode()[0]
print(f"Modal best params   :  n_estimators={modal_n_est}  max_features={modal_mf}")

# ── Feature importance plot (full-data fit) ───────────────────────────────────

scaler_full  = StandardScaler()
X_scaled_all = scaler_full.fit_transform(X_line)

rf_full = RandomForestRegressor(
    n_estimators=modal_n_est,
    max_features=modal_mf,
    random_state=RANDOM_STATE,
    n_jobs=1,
)
rf_full.fit(X_scaled_all, y_line)
importances = rf_full.feature_importances_   # shape (1723,)

fig, ax = plt.subplots(figsize=(9, 3.5))
ax.plot(wavenumbers, importances, lw=0.9, color="#d6604d")
ax.set_xlabel("Wavenumber (cm⁻¹)")
ax.set_ylabel("Mean decrease in impurity")
ax.set_title(
    f"Random forest feature importances — spectral drivers of starvation-resistance prediction\n"
    f"n_estimators={modal_n_est}  max_features={modal_mf}  "
    f"n={n_lines} DGRP line-mean spectra  (full-data fit)"
)
fig.tight_layout()
out_imp = OUT_DIR / "rf_feature_importance_vs_wavenumber.pdf"
fig.savefig(out_imp)
plt.close(fig)
print(f"Saved : {out_imp.relative_to(REPO)}")

# ── Full comparison table ─────────────────────────────────────────────────────

all_rows = PRIOR_RESULTS + [
    ("Random forest", r2_rf, rmse_rf, rho_rf),
]

print()
print("=" * 68)
print("METHOD COMPARISON — LOO-CV on 108 DGRP line-mean spectra")
print("=" * 68)
print(f"  {'Method':<28}  {'CV R²':>6}  {'RMSE':>7}  {'Spearman ρ':>10}")
print("  " + "-" * 57)
for method, r2, rmse, rho in all_rows:
    marker = "  ←" if method == "Random forest" else ""
    print(f"  {method:<28}  {r2:>6.3f}  {rmse:>7.4f}  {rho:>+10.3f}{marker}")
print("=" * 68)
print()
print("Note: PCA+Ridge, PLS, Ridge, LASSO, Elastic net rows hardcoded from")
print("      prior scripts (run_compression_analysis.py, run_pls_analysis.py,")
print("      run_regularised_regression.py).  RF hyperparameters selected by")
print("      inner GridSearchCV(cv=3) inside each outer LOO fold.")
print()
print("Plots saved:")
print(f"  {out_imp.relative_to(REPO)}")
