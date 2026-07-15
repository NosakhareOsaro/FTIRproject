"""
Inspect the Morgante et al. 2015 (DGRPool Study 24) starvation resistance file,
check overlap with our PSM EMMeans, and sanity-check the correlation between
the two independent measures of the same phenotype.
"""

import matplotlib
matplotlib.use("Agg")

import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt

root = Path(__file__).parent.parent

morgante_path = root / "phenotype-data" / "S24_StarvationRes_summary_mean.tsv"
emmeans_path  = root / "Emmeans.csv"
out_path      = root / "results" / "DGRP" / "emmeans_vs_morgante_correlation.pdf"

# Step 1: inspect the raw file
print("=" * 60)
print("STEP 1 — Raw file inspection")
print("=" * 60)

morg_raw = pd.read_csv(morgante_path, sep="\t")
print(f"Columns      : {list(morg_raw.columns)}")
print(f"Total rows   : {len(morg_raw)}")
print(f"Sex values   : {sorted(morg_raw['sex'].unique())}")
print(f"Unique lines : {morg_raw['DGRP'].nunique()}")
print(f"\nFirst 5 rows:")
print(morg_raw.head().to_string(index=False))
print(f"\nDGRP ID format examples: {morg_raw['DGRP'].iloc[:3].tolist()}")

# Step 2: overlap with our 108 EMMeans lines
print()
print("=" * 60)
print("STEP 2 — Line overlap with Emmeans.csv")
print("=" * 60)

emmeans = pd.read_csv(emmeans_path)
print(f"EMMeans lines : {len(emmeans)}")
print(f"EMMeans ID format examples: {emmeans['DGRP'].iloc[:3].tolist()}")

# Normalise: Morgante uses DGRP_021 (zero-padded), ours uses DGRP21.
# Stripping the underscore alone is not enough -- DGRP_021 -> DGRP021 still
# fails to match DGRP21. Round-tripping through int() drops the leading
# zeros too (same normalisation as scripts/run_dgrpool_phenotype.py).
def _normalise_dgrp_id(s):
    return "DGRP" + str(int(s.replace("DGRP", "").replace("_", "")))

morg_f = morg_raw[morg_raw["sex"] == "F"].copy()
morg_f["DGRP_norm"] = morg_f["DGRP"].apply(_normalise_dgrp_id)
emmeans["DGRP_norm"] = emmeans["DGRP"].apply(_normalise_dgrp_id)

our_lines     = set(emmeans["DGRP_norm"])
morgante_lines = set(morg_f["DGRP_norm"])
overlap       = our_lines & morgante_lines

print(f"\nOur lines (EMMeans)        : {len(our_lines)}")
print(f"Morgante female lines      : {len(morgante_lines)}")
print(f"Overlap                    : {len(overlap)}")

morg_overlap = morg_f[morg_f["DGRP_norm"].isin(overlap)].copy()
valid = morg_overlap["value"].notna()
print(f"Of those, with non-NA value: {valid.sum()}")

missing_from_morgante = our_lines - morgante_lines
if missing_from_morgante:
    print(f"Our lines absent in Morgante ({len(missing_from_morgante)}): "
          f"{sorted(missing_from_morgante)}")

vals = morg_overlap.loc[valid, "value"]
print(f"\nMorgante female means (overlapping lines):")
print(f"  Min    : {vals.min():.2f}")
print(f"  Median : {vals.median():.2f}")
print(f"  Max    : {vals.max():.2f}")
print(f"  SD     : {vals.std():.2f}")

# Step 3: correlation
print()
print("=" * 60)
print("STEP 3 — Correlation: PSM EMMeans vs Morgante female means")
print("=" * 60)

# Build merged frame on the overlapping lines
merged = (
    emmeans[["DGRP_norm", "emmean"]]
    .merge(
        morg_overlap.loc[valid, ["DGRP_norm", "value"]],
        on="DGRP_norm",
        how="inner",
    )
    .rename(columns={"emmean": "psm_emmean", "value": "morgante_mean"})
)
print(f"Lines in merged frame: {len(merged)}")

pr, pp = pearsonr(merged["psm_emmean"], merged["morgante_mean"])
sr, sp = spearmanr(merged["psm_emmean"], merged["morgante_mean"])
print(f"\nPearson  r  = {pr:.3f}  (p = {pp:.2e})")
print(f"Spearman rho = {sr:.3f}  (p = {sp:.2e})")

if pr > 0.4 and pp < 0.05:
    verdict = "POSITIVE and significant — the two measures agree."
elif pr > 0 and pp < 0.05:
    verdict = "Weakly positive and significant."
elif pr <= 0:
    verdict = "WARNING: negative or zero correlation — unexpected, investigate."
else:
    verdict = "Positive but not significant — small overlap or noisy data."
print(f"\nVerdict: {verdict}")

# Scatter plot
fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(merged["psm_emmean"], merged["morgante_mean"],
           color="steelblue", alpha=0.7, edgecolors="white", linewidths=0.4)
ax.set_xlabel("PSM EMMean (log-scale survival, our R model)", fontsize=12)
ax.set_ylabel("Morgante mean starvation resistance (hours)", fontsize=12)
ax.set_title(
    f"EMMeans vs Morgante et al. 2015  (n={len(merged)} lines)\n"
    f"Pearson r = {pr:.3f}, Spearman ρ = {sr:.3f}",
    fontsize=11,
)
# Annotate a few extreme points for orientation
for _, row in merged.nlargest(3, "psm_emmean").iterrows():
    ax.annotate(row["DGRP_norm"], (row["psm_emmean"], row["morgante_mean"]),
                fontsize=7, xytext=(4, 2), textcoords="offset points")
for _, row in merged.nsmallest(3, "psm_emmean").iterrows():
    ax.annotate(row["DGRP_norm"], (row["psm_emmean"], row["morgante_mean"]),
                fontsize=7, xytext=(4, 2), textcoords="offset points")

fig.tight_layout()
fig.savefig(out_path)
print(f"\nScatter plot saved to: {out_path.relative_to(root)}")
