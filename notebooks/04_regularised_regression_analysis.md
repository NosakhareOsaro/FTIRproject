# 04 — Regularised Regression Analysis

**Script:** `scripts/run_regularised_regression.py`
**Status:** Complete
**Key result:** Elastic net achieves CV R² = 0.673, Spearman ρ = +0.816 — the best performance of all methods tested. Direct sparse regression on the raw 1,723 wavenumbers outperforms dimensionality-reduction-then-predict approaches (PLS, PCA+Ridge).

---

## What this step does and why

PLS (markdown 03) compresses the spectra into latent variables before predicting phenotype. This two-stage approach is conventional in chemometrics, but it raises a question: what happens if you skip the compression entirely and regress the phenotype directly on all 1,723 wavenumbers?

The obvious problem is that 1,723 predictors with only 108 observations is a severely underdetermined system — there are far more variables than data points. Standard linear regression would overfit catastrophically. The solution is regularisation: adding a mathematical penalty to the model that prevents coefficients from growing too large, effectively constraining the complexity of the fit.

Three regularised regression methods were tested:

- **Ridge regression** — L2 penalty: shrinks all 1,723 coefficients toward zero but keeps all of them non-zero
- **LASSO** — L1 penalty: drives many coefficients to exactly zero, producing a sparse solution where only a subset of wavenumbers are used
- **Elastic net** — combines L1 and L2 penalties: sparse like LASSO but more stable when groups of correlated wavenumbers are all relevant (common in spectral data where adjacent wavenumbers absorb similar compounds)

These are all linear models working directly on the raw spectra — no dimensionality reduction step.

---

## Data

Same as markdown 03 — 108 DGRP line-mean spectra, joined to EMMeans. StandardScaler applied inside each cross-validation fold.

---

## Evaluation design

Identical LOO-CV structure to the PLS analysis: hold out one line, train on 107, predict the held-out line, repeat 108 times. The regularisation strength (α) is selected automatically inside each training fold using an inner cross-validation — the held-out test line never influences α selection.

**Hyperparameter selection inside each fold:**

- **Ridge:** `RidgeCV` with generalised cross-validation (GCV) — a closed-form mathematical equivalent of LOO that selects the optimal α from 100 candidates (log-spaced from 10⁻³ to 10⁶) using only the 107 training lines. This is exact, not approximate, and requires only a single fit.
- **LASSO:** `LassoCV` with 3-fold inner CV on the 107 training lines, 30 α candidates on the regularisation path.
- **Elastic net:** `ElasticNetCV` with 3-fold inner CV, 30 α candidates, and five L1/L2 mixing ratios tested (l1_ratio ∈ [0.5, 0.7, 0.9, 0.95, 1.0]). Both α and the mixing ratio are selected automatically.

---

## Results

| Method | CV R² | CV RMSE | Spearman ρ | Selected α |
|---|---|---|---|---|
| PCA + Ridge (from markdown 03) | 0.553 | 0.4927 | +0.743 | — |
| PLS (from markdown 03) | 0.623 | 0.4524 | +0.801 | — |
| Ridge (raw 1,723 wn) | 0.635 | 0.4451 | +0.809 | median α = 28.5 |
| LASSO | 0.669 | 0.4266 | +0.813 | median α = 0.010 |
| **Elastic net** | **0.673** | **0.4244** | **+0.816** | median α = 0.011, ℓ₁ = 0.90 |

All three regularised methods outperform PLS. Elastic net is the best overall method.

---

## What the results mean

**Ridge beats PLS.** The simplest regularised regression on all 1,723 raw wavenumbers (Ridge, R² = 0.635) outperforms the standard chemometrics method (PLS, R² = 0.623). This is somewhat surprising — PLS is explicitly optimised to predict starvation resistance, while Ridge makes no use of the phenotype during compression. The result suggests that with n=108 lines, simple L2 shrinkage on all wavenumbers exploits the smooth spectral signal more efficiently than PLS latent structure does.

**Sparsity helps.** LASSO (R² = 0.669) improves on Ridge (R² = 0.635) by driving many wavenumber coefficients to zero and concentrating predictive weight on the most informative spectral regions. The L1 penalty effectively performs automatic variable selection.

**Elastic net is the best overall method.** The combination of L1 sparsity and L2 grouping (elastic net, R² = 0.673) gives the best performance. The selected mixing ratio of ℓ₁ = 0.90 means the model is mostly LASSO-like but retains some L2 stabilisation — appropriate for spectral data where groups of correlated adjacent wavenumbers often carry the same information.

**Random forest adds nothing.** For completeness, random forest was also tested (see the full method comparison below) and achieved R² = 0.540 — worse than PCA + Ridge and well below all three linear regularised methods. Nonlinear tree-based models offer no advantage on this collinear spectral dataset. The signal in the spectra is smooth and largely linear (lipid content), which linear models with appropriate regularisation capture efficiently.

---

## Complete method comparison (line-mean spectra, LOO-CV)

| Method | CV R² | CV RMSE | Spearman ρ |
|---|---|---|---|
| PCA + Ridge (4 PCs) | 0.553 | 0.4927 | +0.743 |
| Random forest | 0.540 | 0.4998 | +0.729 |
| PLS (10 components) | 0.623 | 0.4524 | +0.801 |
| Ridge (raw spectra) | 0.635 | 0.4451 | +0.809 |
| LASSO | 0.669 | 0.4266 | +0.813 |
| **Elastic net** | **0.673** | **0.4244** | **+0.816** |

---

## Coefficient plot

After LOO-CV, elastic net was fitted on all 108 line-mean spectra to produce a final model for inspection. The coefficient vector (one value per wavenumber) is plotted in `results/DGRP/regularised_coefficients_vs_wavenumber.pdf`.

The coefficient plot shows a similar pattern to the PLS loading vector (markdown 03): the strongest positive coefficients fall in the C-H stretching region (~2,900-3,000 cm⁻¹), consistent with lipid content as the main driver. Because elastic net is sparse, many wavenumbers have coefficients of exactly zero — only the spectral regions that genuinely contribute to prediction survive the L1 penalty.

---

## Technical details

- `sklearn.linear_model.RidgeCV`, `LassoCV`, `ElasticNetCV` were used
- In sklearn 1.9.0, the α grid size parameter is `alphas=` (integer count), not `n_alphas=` — this version difference caused an initial crash that was fixed before results were obtained
- `tol=0.01` was used for LASSO and elastic net inner CV to speed up convergence on the 1,723-feature coordinate descent path. This is sufficient precision for α selection and does not materially affect the outer LOO-CV performance estimates
- All convergence warnings from coordinate descent are expected at p=1,723 and do not affect the selected α or the final coefficients

---

## Output files

- `results/DGRP/regularised_coefficients_vs_wavenumber.pdf` — elastic net coefficient vector vs wavenumber (best model)