"""
run_pls_analysis.py

Supervised counterpart to run_compression_analysis.py.  PLS regression
explicitly maximises covariance between spectra and starvation EMMeans,
whereas PCA only maximises spectral variance.

Parts:
  A — PLSRegression with LOO-CV over n_components ∈ {1,2,3,5,10}.
      Full-data fit for: LV1 scores vs EMMean scatter; LV1 weight vector plot.
  B — PCA + Ridge LOO-CV (same folds, same scaler discipline) for direct
      comparison on identical evaluation footing.

All CV uses leave-one-line-out (LeaveOneOut over 108 DGRP lines).
Scaling is fitted inside each fold on the 107 training lines to prevent
test-point information leaking into the normalisation.

Outputs written to results/DGRP/:
  pls_component1_vs_emmean.pdf
  pls_loading1_vs_wavenumber.pdf
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))
from ftir_loader import load_ftir

OUT_DIR = REPO / "results" / "DGRP"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_PCS        = 4                    # PCA: 95% variance threshold (established in compression analysis)
N_COMP_LIST  = [1, 2, 3, 5, 10]    # PLS: sweep over these
RIDGE_ALPHA  = 1.0                  # Ridge alpha — not tuned; purpose is to evaluate PCA representation

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

wavenumbers = np.array(spectra.columns.tolist())  # [3900, 3898, …, 456]

# Average spectra within each DGRP line (Order B)
spec_df = spectra.copy()
spec_df["Genot."] = meta["Genot."].values
X_line_df = spec_df.groupby("Genot.").mean()   # 108 × 1723, index = Genot.

emmeans  = pd.read_csv(REPO / "Emmeans.csv")
em_vmin  = emmeans["emmean"].min()
em_vmax  = emmeans["emmean"].max()

line_data = (
    X_line_df.reset_index()
    .merge(emmeans[["DGRP", "emmean"]], left_on="Genot.", right_on="DGRP", how="inner")
    .sort_values("Genot.")
    .reset_index(drop=True)
)

X_line   = line_data[wavenumbers.tolist()].values   # 108 × 1723
y_line   = line_data["emmean"].values               # 108
n_lines  = len(y_line)

print(f"Lines in analysis : {n_lines}  (DGRP lines with both spectra and EMMeans)")
print(f"Spectral features : {X_line.shape[1]} wavenumbers")

# ═══════════════════════════════════════════════════════════════════════════════
# Part A — PLS regression with LOO-CV
# ═══════════════════════════════════════════════════════════════════════════════

loo = LeaveOneOut()

def loo_cv(X, y, model_fn):
    """Run LOO-CV with per-fold StandardScaling.  Returns y_pred array."""
    y_pred = np.zeros(len(y))
    for train_idx, test_idx in loo.split(X):
        X_tr, X_te = X[train_idx], X[test_idx]
        y_tr = y[train_idx]
        sc = StandardScaler()
        X_tr_s = sc.fit_transform(X_tr)
        X_te_s = sc.transform(X_te)
        mdl = model_fn()
        mdl.fit(X_tr_s, y_tr)
        y_pred[test_idx] = np.asarray(mdl.predict(X_te_s)).ravel()
    return y_pred


def metrics(y_true, y_pred):
    r2   = r2_score(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    rho  = spearmanr(y_true, y_pred).statistic
    return r2, rmse, rho


print()
print("Running PLS LOO-CV …")
pls_results = {}
for k in N_COMP_LIST:
    y_pred_k = loo_cv(X_line, y_line, lambda k=k: PLSRegression(n_components=k))
    pls_results[k] = dict(zip(("r2", "rmse", "rho"), metrics(y_line, y_pred_k)))
    print(f"  n_components={k:2d}  R²={pls_results[k]['r2']:+.3f}  "
          f"RMSE={pls_results[k]['rmse']:.4f}  ρ={pls_results[k]['rho']:+.3f}")

opt_k = max(N_COMP_LIST, key=lambda k: pls_results[k]["r2"])
print(f"  → optimal n_components = {opt_k}")

# ── Full-data PLS fit (for plots) ─────────────────────────────────────────────

scaler_full  = StandardScaler()
X_scaled_all = scaler_full.fit_transform(X_line)

pls_full = PLSRegression(n_components=opt_k)
pls_full.fit(X_scaled_all, y_line)

lv1_scores   = pls_full.transform(X_scaled_all)[:, 0]   # LV1 score per line
weights_lv1  = pls_full.x_weights_[:, 0]                # weight vector: which WNs drive LV1
r_insample, _ = pearsonr(lv1_scores, y_line)

# ── Plot 1: LV1 scores vs EMMean ──────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(6, 5.5))
sc = ax.scatter(
    lv1_scores, y_line,
    c=y_line,
    cmap="plasma",
    vmin=em_vmin, vmax=em_vmax,
    s=60, alpha=0.85, linewidths=0.4, edgecolors="white",
)
fig.colorbar(sc, ax=ax, label="Starvation resistance (EMMean)")
ax.set_xlabel(f"PLS LV1 score  (n_components = {opt_k})")
ax.set_ylabel("Starvation resistance (EMMean)")
ax.set_title(
    f"PLS regression — LV1 score vs starvation EMMean\n"
    f"n = {n_lines} DGRP lines  (full-data fit)"
)
ax.text(
    0.03, 0.97,
    f"r = {r_insample:+.3f}  (in-sample)\n"
    f"CV R² = {pls_results[opt_k]['r2']:+.3f}  (LOO)",
    transform=ax.transAxes, va="top", ha="left", fontsize=9,
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, edgecolor="none"),
)
fig.tight_layout()
out1 = OUT_DIR / "pls_component1_vs_emmean.pdf"
fig.savefig(out1)
plt.close(fig)
print(f"Saved : {out1.relative_to(REPO)}")

# ── Plot 2: LV1 weight vector vs wavenumber ───────────────────────────────────

fig, ax = plt.subplots(figsize=(9, 3.5))
ax.plot(wavenumbers, weights_lv1, lw=0.9, color="#2166ac")
ax.axhline(0, color="#888888", lw=0.7, ls="--")
ax.set_xlabel("Wavenumber (cm⁻¹)")
ax.set_ylabel("PLS weight  (component 1)")
ax.set_title(
    f"PLS LV1 weight vector — spectral drivers of starvation-resistance prediction\n"
    f"n_components = {opt_k},  n = {n_lines} DGRP line-mean spectra"
)
fig.tight_layout()
out2 = OUT_DIR / "pls_loading1_vs_wavenumber.pdf"
fig.savefig(out2)
plt.close(fig)
print(f"Saved : {out2.relative_to(REPO)}")

# ═══════════════════════════════════════════════════════════════════════════════
# Part B — PCA + Ridge LOO-CV (same evaluation footing)
# ═══════════════════════════════════════════════════════════════════════════════

print()
print(f"Running PCA ({N_PCS} PCs) + Ridge (α={RIDGE_ALPHA}) LOO-CV …")

def pca_ridge():
    class _Model:
        def fit(self, X, y):
            self._pca = PCA(n_components=N_PCS).fit(X)
            self._ridge = Ridge(alpha=RIDGE_ALPHA).fit(self._pca.transform(X), y)
        def predict(self, X):
            return self._ridge.predict(self._pca.transform(X))
    return _Model()

y_pred_pca = loo_cv(X_line, y_line, pca_ridge)
r2_pca, rmse_pca, rho_pca = metrics(y_line, y_pred_pca)
print(f"  R²={r2_pca:+.3f}  RMSE={rmse_pca:.4f}  ρ={rho_pca:+.3f}")

# ═══════════════════════════════════════════════════════════════════════════════
# Summary
# ═══════════════════════════════════════════════════════════════════════════════

opt = pls_results[opt_k]

print()
print("=" * 66)
print("PLS LOO-CV performance by n_components:")
print(f"  {'n_comp':>6}  {'CV R²':>7}  {'CV RMSE':>8}  {'Spearman ρ':>10}")
print("  " + "-" * 38)
for k in N_COMP_LIST:
    marker = "  ← optimal" if k == opt_k else ""
    r = pls_results[k]
    print(f"  {k:>6}  {r['r2']:>7.3f}  {r['rmse']:>8.4f}  {r['rho']:>+10.3f}{marker}")

print()
print("=" * 66)
print("METHOD COMPARISON — LOO-CV on 108 DGRP line-mean spectra")
print("=" * 66)
print(f"  {'Method':<34} {'n_comp':>6}  {'CV R²':>6}  {'RMSE':>7}  {'Spearman ρ':>10}")
print("  " + "-" * 62)
print(f"  {'PCA + Ridge  (α=1.0, 95% var)':<34} {N_PCS:>6}  {r2_pca:>6.3f}  {rmse_pca:>7.4f}  {rho_pca:>+10.3f}")
print(f"  {'PLS regression  (optimal)':<34} {opt_k:>6}  {opt['r2']:>6.3f}  {opt['rmse']:>7.4f}  {opt['rho']:>+10.3f}")
print("=" * 66)
print()
print("Note: PCA+Ridge α not tuned — purpose is to evaluate the PCA")
print("      representation, not to optimise Ridge hyperparameters.")
print()
print("Plots saved:")
print(f"  {out1.relative_to(REPO)}")
print(f"  {out2.relative_to(REPO)}")
