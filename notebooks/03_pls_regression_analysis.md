# 03 — PLS Regression Analysis

**Script:** `scripts/run_pls_analysis.py`
**Status:** Complete
**Key result:** PLS regression (10 components, LOO-CV) achieves CV R² = 0.623, Spearman ρ = +0.801 on 108 DGRP line-mean spectra. The LV1 weight vector identifies the C-H stretching region (~2,900-3,000 cm⁻¹) as the main spectral driver — consistent with lipid content as the biochemical correlate of starvation resistance.

---

## What this step does and why

PCA (markdown 02) is unsupervised — it finds directions of maximum spectral variance without any knowledge of the starvation phenotype. PC1 ended up correlating with starvation resistance (r = +0.685), but this was because the dominant source of between-line spectral variation happens to be lipid content, which also predicts starvation resistance. There is no guarantee PCA will find phenotypically relevant directions in general.

Partial Least Squares (PLS) regression is the supervised counterpart. Instead of maximising spectral variance, PLS finds latent directions that maximise covariance between the spectra and the phenotype. It compresses the 1,723 wavenumbers into a small number of latent variables (LVs) that are optimally aligned with the starvation resistance target. PLS is the standard method in the metabolomics and chemometrics field for exactly this type of problem — high-dimensional spectra, continuous phenotype, small sample size.

This step answers two questions:

1. Does a supervised compression (PLS) improve on unsupervised compression (PCA + Ridge regression) for predicting starvation resistance?
2. Which spectral regions drive the prediction — and does this make biological sense?

---

## Data

108 DGRP line-mean spectra (Order B from markdown 02 — spectra averaged within each line first). Joined to `Emmeans.csv` to get one starvation EMMean per line. StandardScaler applied to the spectra inside each cross-validation fold (fitted on training lines only — see evaluation design below).

---

## Evaluation design — leave-one-line-out cross-validation (LOO-CV)

With only 108 lines, the choice of cross-validation strategy matters. K-fold CV with k=10 would give folds of roughly 10 lines each — too small for stable performance estimates. LOO-CV holds out one line at a time and trains on the remaining 107, repeating this 108 times. Each line is predicted exactly once by a model that never saw it.

This is computationally more expensive than k-fold but gives the most honest performance estimate at this sample size.

**Inside each LOO fold:**
1. StandardScaler fitted on the 107 training lines (not the held-out test line)
2. PLSRegression fitted on the scaled training lines
3. Prediction made for the single held-out test line

After all 108 folds, performance metrics are computed on the 108 held-out predictions versus the true EMMeans.

---

## Model selection — how many components?

PLS has one main hyperparameter: the number of latent components. Too few and the model underfits; too many and it overfits on the small sample. LOO-CV was run across five values of n_components:

| n_components | CV R² | Spearman ρ |
|---|---|---|
| 1 | 0.479 | +0.682 |
| 2 | 0.534 | +0.707 |
| 3 | 0.578 | +0.761 |
| 5 | 0.614 | +0.787 |
| **10** | **0.623** | **+0.801** |
| 15 | 0.550 | +0.784 |
| 20 | 0.515 | +0.750 |

Performance peaks at n=10 and then declines. This is a classic overfitting pattern — with n=108 lines, adding more than 10 latent components gives the model too many degrees of freedom to fit training noise. The decline at n=15 and n=20 is the model memorising the 107 training lines rather than learning a signal that generalises.

**n=10 is confirmed as the optimal number of components.** This was found empirically, not assumed.

---

## Results

| Metric | Value |
|---|---|
| CV R² (LOO) | 0.623 |
| CV RMSE (LOO) | 0.4524 |
| Spearman ρ (LOO) | +0.801 |
| In-sample r (full-data fit) | +0.714 |

The model explains 62.3% of between-line variation in starvation resistance from the FTIR spectra alone, on held-out lines. The Spearman rank correlation of +0.801 means that when lines are ranked by their predicted starvation resistance, the ranking closely matches their true ranking.

---

## Comparison with PCA + Ridge

To put PLS on the same evaluation footing as PCA, a PCA + Ridge regression baseline was also run with LOO-CV. PCA (4 components, 95% variance threshold from markdown 02) was fitted on the 107 training lines inside each fold, and Ridge regression was fitted on the 4 PC scores.

| Method | CV R² | Spearman ρ |
|---|---|---|
| PCA + Ridge (4 PCs) | 0.553 | +0.743 |
| PLS (10 components) | **0.623** | **+0.801** |

PLS outperforms PCA + Ridge. This is expected — PLS uses the phenotype labels during compression, so its latent variables are by construction more aligned with starvation resistance than the unsupervised PCA components. The gap (R² +0.07) is moderate but meaningful.

---

## Spectral interpretation — what drives the prediction?

The PLS weight vector for latent variable 1 (LV1) shows which wavenumbers contribute most to the first latent variable — the component most directly predictive of starvation resistance.

**`results/DGRP/pls_loading1_vs_wavenumber.pdf`** shows this weight vector plotted against wavenumber. The dominant feature is a strong positive peak in the **C-H stretching region at approximately 2,900-3,000 cm⁻¹**. This region corresponds to the absorption of alkyl C-H bonds in lipids and fatty acids.

Additional features appear at:
- ~1,650 cm⁻¹ (Amide I — protein C=O stretching, negative weight)
- ~1,550 cm⁻¹ (Amide II — protein N-H bending, positive weight)
- ~1,000-1,200 cm⁻¹ (carbohydrate and polysaccharide fingerprint region, moderate negative weights)

The loading vector is smooth and interpretable — it does not look like random noise. This is a good sign that the model is capturing real chemical structure rather than overfitting to measurement artefacts.

**Biological interpretation:** Lines with higher starvation resistance appear to have greater lipid content in their cuticular chemical fingerprint, as measured by FTIR. This makes strong biological sense. Fat reserves (primarily stored as triacylglycerols) are the primary energy source during starvation in Drosophila. Lines that can store more lipid are better buffered against starvation. The FTIR measurement at the cuticle surface picks up this variation in lipid content and uses it to predict how long flies will survive without food.

This connects directly to Adam's comment after seeing the results: the FTIR prediction method may be most useful for traits where lipid metabolism is involved — an observation that is directly supported by the loading vector showing the lipid band as the dominant spectral feature.

---

## Technical details

- `sklearn.cross_decomposition.PLSRegression` was used throughout
- `PLSRegression.predict()` gives the scalar regression output used for CV metrics
- `PLSRegression.transform()` gives the LV scores used for the scatter plot — these are different and easy to confuse in the sklearn API
- `x_weights_` (not `x_loadings_`) was used for the weight vector plot — the weight vector is what PLS actually uses to project wavenumbers into the latent space; loadings are a derived quantity

---

## Output files

- `results/DGRP/pls_component1_vs_emmean.pdf` — LV1 scores vs EMMean scatter (in-sample, full-data fit; r = +0.714 annotated; CV R² = 0.623 also shown)
- `results/DGRP/pls_loading1_vs_wavenumber.pdf` — LV1 weight vector vs wavenumber (C-H peak visible at ~2,900-3,000 cm⁻¹)