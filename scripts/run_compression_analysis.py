"""
run_compression_analysis.py

Stage 1 of the two-stage pipeline: compress DGRP female FTIR spectra
via PCA and examine alignment with starvation-resistance EMMeans.

Parts:
  A: PCA on individual fly spectra; explained-variance curve and
      PC1 vs PC2 scatter coloured by EMMean.
  B: Compare two collapse orderings:
        A: project all fly spectra to PCs, then average per DGRP line
        B: average raw spectra per DGRP line first, then run PCA
      Pearson r between PC1 and EMMeans is reported for both.
  C: Print summary.

Outputs written to results/DGRP/:
  pca_explained_variance.pdf
  pca_coloured_by_emmean.pdf
  pca_linemeans_coloured_by_emmean.pdf
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # must precede pyplot import
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))
from ftir_loader import load_ftir

OUT_DIR = REPO / "results" / "DGRP"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_COMPONENTS = 100   # upper cap; 95% is typically reached in ~10–20 PCs for spectral data
VAR_THRESHOLD = 0.95

plt.rcParams.update({
    "figure.dpi": 150,
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# ── Load data ─────────────────────────────────────────────────────────────────

meta, spectra = load_ftir(REPO / "FTIR-data" / "DGRPFTIR.dat")
n_total = len(meta)
female_mask = meta["Sex"] == "F"
meta    = meta[female_mask].reset_index(drop=True)
spectra = spectra[female_mask].reset_index(drop=True)
n_female = len(meta)

print(f"Spectra loaded : {n_total} total → {n_female} female (Sex=F)")

emmeans = pd.read_csv(REPO / "Emmeans.csv")
em_vmin = emmeans["emmean"].min()   # shared colour scale across all EMMean scatter plots
em_vmax = emmeans["emmean"].max()

# ── StandardScale (mean=0, std=1 per wavenumber) ──────────────────────────────

scaler = StandardScaler()
X_scaled = scaler.fit_transform(spectra.values)  # n_flies × 1723

# ═══════════════════════════════════════════════════════════════════════════════
# Part A: PCA on individual fly spectra
# ═══════════════════════════════════════════════════════════════════════════════

pca = PCA(n_components=N_COMPONENTS)
pca.fit(X_scaled)

cumvar   = np.cumsum(pca.explained_variance_ratio_)
n_pcs_95 = int(np.argmax(cumvar >= VAR_THRESHOLD)) + 1  # first PC that clears 95%
pct1     = pca.explained_variance_ratio_[0] * 100
pct2     = pca.explained_variance_ratio_[1] * 100

# ── Plot 1: cumulative explained variance ─────────────────────────────────────

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(np.arange(1, N_COMPONENTS + 1), cumvar * 100, lw=1.8, color="#2166ac")
ax.axhline(95, color="#d6604d", lw=1.2, ls="--", label="95% threshold")
ax.axvline(n_pcs_95, color="#d6604d", lw=1.2, ls=":")
ax.annotate(
    f"PC {n_pcs_95}  ({cumvar[n_pcs_95 - 1] * 100:.1f}%)",
    xy=(n_pcs_95, cumvar[n_pcs_95 - 1] * 100),
    xytext=(n_pcs_95 + 4, cumvar[n_pcs_95 - 1] * 100 - 6),
    arrowprops=dict(arrowstyle="->", color="#333333", lw=0.9),
    fontsize=9,
)
ax.set_xlabel("Number of principal components")
ax.set_ylabel("Cumulative explained variance (%)")
ax.set_title(f"PCA on DGRP female FTIR spectra (n = {n_female})")
ax.set_xlim(1, N_COMPONENTS)
ax.set_ylim(0, 102)
ax.legend(fontsize=9)
fig.tight_layout()
out1 = OUT_DIR / "pca_explained_variance.pdf"
fig.savefig(out1)
plt.close(fig)
print(f"Saved : {out1.relative_to(REPO)}")

# ── Project all flies to PC space ─────────────────────────────────────────────

scores_fly = pca.transform(X_scaled)  # n_flies × N_COMPONENTS

# Join EMMeans to each fly via its DGRP line
meta_em = meta[["Genot."]].merge(
    emmeans[["DGRP", "emmean"]],
    left_on="Genot.", right_on="DGRP",
    how="left",
)
n_missing = meta_em["emmean"].isna().sum()
if n_missing:
    print(f"  Warning: {n_missing} flies have no EMMean — excluded from colour plot")

scatter_mask = meta_em["emmean"].notna().values

# ── Plot 2: PC1 vs PC2 coloured by EMMean ────────────────────────────────────

fig, ax = plt.subplots(figsize=(7, 5.5))
sc = ax.scatter(
    scores_fly[scatter_mask, 0],
    scores_fly[scatter_mask, 1],
    c=meta_em.loc[scatter_mask, "emmean"].values,
    cmap="plasma",
    vmin=em_vmin,
    vmax=em_vmax,
    s=8,
    alpha=0.75,
    linewidths=0,
)
fig.colorbar(sc, ax=ax, label="Starvation resistance (EMMean)")
ax.set_xlabel(f"PC1 ({pct1:.1f}% var)")
ax.set_ylabel(f"PC2 ({pct2:.1f}% var)")
ax.set_title(
    "DGRP female FTIR spectra — PC1 vs PC2\n"
    "Coloured by starvation resistance EMMean"
)
fig.tight_layout()
out2 = OUT_DIR / "pca_coloured_by_emmean.pdf"
fig.savefig(out2)
plt.close(fig)
print(f"Saved : {out2.relative_to(REPO)}")

# ═══════════════════════════════════════════════════════════════════════════════
# Part B: Two orderings of collapse
# ═══════════════════════════════════════════════════════════════════════════════

# -- Order A: PCA first, then average PC scores per DGRP line
fly_df = pd.DataFrame(
    scores_fly,
    columns=[f"PC{i + 1}" for i in range(N_COMPONENTS)],
)
fly_df["Genot."] = meta["Genot."].values
line_A = fly_df.groupby("Genot.").mean().reset_index()
line_A = line_A.merge(
    emmeans[["DGRP", "emmean"]], left_on="Genot.", right_on="DGRP", how="inner"
)
r_A, p_A = pearsonr(line_A["PC1"], line_A["emmean"])

# -- Order B: average raw (scaled) spectra per DGRP line, then PCA on 108 line means
spec_df = pd.DataFrame(X_scaled, columns=spectra.columns)
spec_df["Genot."] = meta["Genot."].values
X_line = spec_df.groupby("Genot.").mean()          # 108 × 1723, index = Genot.
n_lines = len(X_line)
pca_B    = PCA(n_components=min(n_lines - 1, N_COMPONENTS))
scores_B = pca_B.fit_transform(X_line.values)       # 108 × n_components
line_B   = pd.DataFrame({
    "Genot.": X_line.index.values,
    "PC1_B":  scores_B[:, 0],
    "PC2_B":  scores_B[:, 1],
})
line_B   = line_B.merge(
    emmeans[["DGRP", "emmean"]], left_on="Genot.", right_on="DGRP", how="inner"
)
r_B, p_B = pearsonr(line_B["PC1_B"], line_B["emmean"])

# ── Plot 3: Order B: line-mean PC1 vs PC2 coloured by EMMean ────────────────

pct1_B = pca_B.explained_variance_ratio_[0] * 100
pct2_B = pca_B.explained_variance_ratio_[1] * 100

fig, ax = plt.subplots(figsize=(7, 5.5))
sc = ax.scatter(
    line_B["PC1_B"],
    line_B["PC2_B"],
    c=line_B["emmean"].values,
    cmap="plasma",
    vmin=em_vmin,
    vmax=em_vmax,
    s=60,
    alpha=0.85,
    linewidths=0.4,
    edgecolors="white",
)
fig.colorbar(sc, ax=ax, label="Starvation resistance (EMMean)")
ax.set_xlabel(f"PC1 ({pct1_B:.1f}% var)")
ax.set_ylabel(f"PC2 ({pct2_B:.1f}% var)")
ax.set_title(
    f"DGRP line-mean FTIR spectra — PC1 vs PC2  (n = {len(line_B)} lines)\n"
    "Order B: average per line first, then PCA"
)
ax.text(
    0.03, 0.97,
    f"r(PC1, EMMean) = {r_B:+.3f}  (p = {p_B:.2e})",
    transform=ax.transAxes,
    va="top", ha="left",
    fontsize=9,
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7, edgecolor="none"),
)
fig.tight_layout()
out3 = OUT_DIR / "pca_linemeans_coloured_by_emmean.pdf"
fig.savefig(out3)
plt.close(fig)
print(f"Saved : {out3.relative_to(REPO)}")

# ═══════════════════════════════════════════════════════════════════════════════
# Part C: Summary
# ═══════════════════════════════════════════════════════════════════════════════

print()
print("=" * 60)
print("PCA COMPRESSION ANALYSIS — DGRPFTIR.dat (female flies)")
print("=" * 60)
print(f"  Female spectra  : {n_female}  across {n_lines} DGRP lines")
print(f"  PCs to 95% var  : {n_pcs_95}  (PC1 = {pct1:.1f}%,  PC2 = {pct2:.1f}%)")
print()
print("PC1 vs starvation EMMean (Pearson r):")
print(f"  Order A — PCA first, average per line  : r = {r_A:+.3f}  p = {p_A:.4f}  n = {len(line_A)}")
print(f"  Order B — average per line, then PCA   : r = {r_B:+.3f}  p = {p_B:.4f}  n = {len(line_B)}")
print()

# Plain-English verdict
sig   = p_A < 0.05 or p_B < 0.05
strong = max(abs(r_A), abs(r_B)) >= 0.30

if strong and sig:
    signal_verdict = (
        "YES — PC1 carries meaningful spectral structure "
        "aligned with starvation resistance."
    )
elif sig:
    signal_verdict = (
        "MARGINAL — PC1 shows weak but statistically "
        "significant alignment with starvation resistance."
    )
else:
    signal_verdict = (
        "NO — PC1 does not significantly align with "
        "starvation resistance at p < 0.05."
    )

diff_r = abs(abs(r_A) - abs(r_B))
if diff_r < 0.03:
    order_verdict = (
        "Order of operations makes no meaningful difference (|Δr| < 0.03).\n"
        "  Averaging flies before or after projection gives equivalent PC1 scores."
    )
else:
    better = "A (PCA → average)" if abs(r_A) > abs(r_B) else "B (average → PCA)"
    order_verdict = (
        f"Order {better} gives stronger PC1–EMMean alignment (|Δr| = {diff_r:.3f})."
    )

print("Verdict:")
print(f"  Spectral signal: {signal_verdict}")
print(f"  Collapse order : {order_verdict}")
print()
print("Note: PC1 direction is sign-arbitrary; |r| is the relevant comparison.")
print()
print("Plots saved:")
print(f"  {out1.relative_to(REPO)}")
print(f"  {out2.relative_to(REPO)}")
print(f"  {out3.relative_to(REPO)}")
