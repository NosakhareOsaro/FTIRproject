"""
run_fecundity_enet.py

Elastic net regression predicting lifetime fecundity from DGRP FTIR spectra.
Second phenotype after starvation resistance — tests whether spectral signal
generalises beyond the primary training-target phenotype.

Phenotype source:
  phenotype-data/S18_LifeFecundity_mean.tsv  (DGRPool study 18, female means)
  Columns: DGRP | sex | value
  Line IDs normalised: DGRP_021 → DGRP21 (strip underscore, drop leading zeros)

Spectral data: DGRPFTIR.dat, female spectra only, averaged per line.
Overlap: 96 of 108 spectral lines have a fecundity mean; 12 dropped (no phenotype).

CV design: LeaveOneOut over 96 lines. StandardScaler fitted on 95 training
lines inside each fold. ElasticNetCV(cv=3, l1_ratio=[0.5,0.9,1.0], alphas=30,
max_iter=5000, tol=0.01) — identical hyperparameter setup to the starvation
resistance analysis in run_regularised_regression.py.

RESULT INTERPRETATION NOTE:
  The model predicts approximately the training mean for every test point
  (SD of predictions ≈ 2.8 vs true SD ≈ 19.9). Elastic net selects very high
  regularisation, driving all coefficients near zero — indicating no detectable
  spectral signal for lifetime fecundity.

  As a consequence, Spearman ρ is UNRELIABLE here and should not be reported:
  with near-constant predictions, the LOO mean-shift artifact dominates (when
  a high-fecundity line is held out, the training mean drops slightly, so the
  model predicts slightly lower → artificial monotone negative trend → ρ → -1).
  This is a numerical artefact, not a biological signal. The honest metric is
  R² = -0.109 (worse than predicting the global mean).
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))
from ftir_loader import load_ftir

L1_RATIOS = [0.5, 0.9, 1.0]
N_ALPHAS  = 30

# ── Load and normalise fecundity phenotype ────────────────────────────────────

fec = pd.read_csv(REPO / "phenotype-data" / "S18_LifeFecundity_mean.tsv", sep="\t")
# All rows are female; normalise DGRP_021 → DGRP21
fec["line"] = fec["DGRP"].str.replace("_", "").apply(
    lambda s: "DGRP" + str(int(s.replace("DGRP", "")))
)
fec_map = fec.set_index("line")["value"].to_dict()

# ── Load spectral data, female only, average per line ────────────────────────

meta, spectra = load_ftir(REPO / "FTIR-data" / "DGRPFTIR.dat")
female_mask = meta["Sex"] == "F"
meta    = meta[female_mask].reset_index(drop=True)
spectra = spectra[female_mask].reset_index(drop=True)

wavenumbers = np.array(spectra.columns.tolist())

spec_df = spectra.copy()
spec_df["Genot."] = meta["Genot."].values
X_line_df = spec_df.groupby("Genot.").mean()   # 108 × 1723

# Inner join: keep only lines present in both spectral data and fecundity file
spectral_lines = set(X_line_df.index)
fec_lines      = set(fec_map)
overlap        = sorted(spectral_lines & fec_lines)
missing_spec   = sorted(spectral_lines - fec_lines)

print(f"Spectral lines (female, 108)    : {len(spectral_lines)}")
print(f"Fecundity lines (DGRPool S18)   : {len(fec_map)}")
print(f"Overlap                         : {len(overlap)}")
print(f"Spectral lines with no fecundity: {len(missing_spec)}")
print(f"  {missing_spec}")
print()

X_line  = X_line_df.loc[overlap, wavenumbers.tolist()].values   # 96 × 1723
y_line  = np.array([fec_map[ln] for ln in overlap])             # 96 fecundity means
n_lines = len(y_line)

print(f"Analysis matrix : {n_lines} lines × {X_line.shape[1]} wavenumbers")
print(f"Fecundity range : {y_line.min():.2f} – {y_line.max():.2f}  "
      f"(mean {y_line.mean():.2f}, SD {y_line.std():.2f})")
print()

# ── LOO-CV ────────────────────────────────────────────────────────────────────

print(f"Running elastic net LOO-CV  "
      f"(l1_ratio={L1_RATIOS}, alphas={N_ALPHAS}, cv=3) …")

loo    = LeaveOneOut()
y_pred = np.zeros(n_lines)

for train_idx, test_idx in loo.split(X_line):
    sc   = StandardScaler()
    X_tr = sc.fit_transform(X_line[train_idx])
    X_te = sc.transform(X_line[test_idx])
    mdl  = ElasticNetCV(
        cv=3, l1_ratio=L1_RATIOS,
        alphas=N_ALPHAS, max_iter=5000, tol=0.01,
    )
    mdl.fit(X_tr, y_line[train_idx])
    y_pred[test_idx] = mdl.predict(X_te)

r2_fec   = r2_score(y_line, y_pred)
rmse_fec = np.sqrt(mean_squared_error(y_line, y_pred))
pred_sd  = y_pred.std()

print(f"  R²={r2_fec:+.3f}  RMSE={rmse_fec:.4f}")
print(f"  Prediction SD = {pred_sd:.3f}  (true SD = {y_line.std():.3f})")
print(f"  Predictions collapse to training mean — no spectral signal detected.")
print()

# ── Cross-phenotype comparison table ─────────────────────────────────────────

print("=" * 68)
print("CROSS-PHENOTYPE COMPARISON — Elastic net LOO-CV, line-mean spectra")
print("=" * 68)
print(f"  {'Phenotype':<26}  {'n lines':>7}  {'CV R²':>8}  {'RMSE':>8}")
print("  " + "-" * 55)
print(f"  {'Starvation resistance':<26}  {108:>7}  {0.673:>8.3f}  {0.4244:>8.4f}")
print(f"  {'Lifetime fecundity':<26}  {n_lines:>7}  {r2_fec:>8.3f}  {rmse_fec:>8.4f}")
print("=" * 68)
print()
print("Starvation resistance row from run_regularised_regression.py")
print("(LOO-CV on 108 lines, same ElasticNetCV hyperparameters).")
print()
print("Fecundity result: model predicts ≈ training mean for all test points")
print(f"(pred SD = {pred_sd:.2f} vs true SD = {y_line.std():.2f}). No spectral signal")
print("detectable for lifetime fecundity. Spearman ρ is not reported: with")
print("near-constant predictions, LOO mean-shift creates an artifactual ρ ≈ -1.")
