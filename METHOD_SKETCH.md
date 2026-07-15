# Method-comparison structure: prep for meeting with Vinny Davies (Mon 29 Jun 2026, 2 pm)

## Framing

The core question is whether FTIR spectra of individual *Drosophila melanogaster*
encode enough information about whole-organism phenotype to predict that phenotype in
flies not seen during training. The data: ~1,700 spectra from 108 inbred DGRP lines
(~16 flies per line), 1,723 wavenumber features, and a continuous starvation-resistance
score per line (EMMean from a parametric survival model) as the primary target. Because
assays are destructive, spectra and phenotype are only linked at the line level, not the
individual fly, which drives a specific cross-validation design (line-stratified, see
§4 of PROJECT_NOTES.md). The question is fundamentally one of high-dimensional regression
on a small, structured dataset (n=108 lines, p=1,723), and the aim is to map the
method space thoroughly rather than to advocate for one approach.

---

## Method table

### (a) Dimensionality reduction only

| Method | What it does | Continuous phenotype output? | Wavenumber interpretability | scikit-learn? |
|--------|-------------|-----------------------------|-----------------------------|---------------|
| PCA | Linear projection onto orthogonal axes of maximum spectral variance | No: requires a downstream regression step on PC scores | Yes: PC loadings reveal which wavenumbers drive each component | `decomposition.PCA` |
| tSNE | Non-linear embedding for 2–3D visualisation of spectral similarity | No: for visualisation only; cannot generalise to new samples | No: stochastic; no stable mapping from wavenumbers to embedding axes | `manifold.TSNE` |

### (b) Both dimensionality reduction and prediction (hybrid)

| Method | What it does | Continuous phenotype output? | Wavenumber interpretability | scikit-learn? |
|--------|-------------|-----------------------------|-----------------------------|---------------|
| PLS-DA | Finds latent components that maximise covariance between spectra and class labels; predicts discrete class from latent scores | Latent scores are continuous; final prediction is discrete (class) | Yes: loadings and VIP scores rank wavenumber contributions | `cross_decomposition.PLSRegression` (DA variant: apply threshold to continuous score) |
| sPLS-DA | PLS-DA with L1 sparsity on loadings, explicitly selecting a subset of wavenumbers per component | Latent scores are continuous; final prediction is discrete | Yes, and sparser than PLS-DA: selected wavenumbers only | Not in scikit-learn; `mixOmics` (R) or custom Python |

### (c) Direct prediction methods

| Method | What it does | Continuous phenotype output? | Wavenumber interpretability | scikit-learn? |
|--------|-------------|-----------------------------|-----------------------------|---------------|
| LASSO | L1-regularised linear regression; drives most coefficients to exactly zero | Yes | Yes: non-zero coefficients identify important wavenumbers | `linear_model.Lasso` |
| Ridge regression | L2-regularised linear regression; shrinks all coefficients toward zero without zeroing them | Yes | Partial: all wavenumbers retained; magnitude indicates relative contribution | `linear_model.Ridge` |
| Elastic net | Mixes L1 and L2 penalties; more stable than LASSO in correlated (spectral) data | Yes | Yes: like LASSO but with group selection behaviour in correlated regions | `linear_model.ElasticNet` |
| Random forest | Ensemble of decision trees trained on bootstrap samples; prediction is average of trees | Yes | Partial: feature importances give wavenumber ranking but not directionality | `ensemble.RandomForestRegressor` |
| SVR | Support vector machine adapted for regression; uses ε-insensitive loss with kernel trick | Yes | No (with RBF kernel) | `svm.SVR` (already used in the classification baseline as SVM) |

---

## Evaluation

Two parallel evaluation settings, both required:

**(a) Per-fly spectra with line-stratified cross-validation**: ~1,700 spectra as
training units; GroupKFold splits by DGRP line so no fly's line appears in both train
and test. Per-fly predictions within each test fold are averaged to the line level before
computing metrics. This is the primary setting and the most important to get right
(see §4 of PROJECT_NOTES.md for the correctness argument).

**(b) Line-mean spectra**: collapse the ~16 per-line spectra to a single mean before
fitting; train and evaluate on the 108 line means directly. Lower N but conceptually
simpler; serves as a sanity check against setting (a).

**Metrics for continuous targets:** R² (explained variance), RMSE, Spearman ρ (rank
correlation, less sensitive to outlier lines). All computed at the line level.

**Positive control:** starvation resistance EMMeans, with signal confirmed in the pre-print
and in our reproduction (r=0.46 vs Morgante 2015, 93/108 lines overlap). If no method
achieves meaningful R² on starvation resistance, the spectral signal is absent and
nothing downstream is credible.

---

## Open questions for Vinny

1. **Comparing across method groups:** Methods in group (a) do not produce phenotype
   predictions directly: they require a second-stage regression on PC or embedding
   scores. Their performance therefore depends on how many components are retained and
   which regression method is applied in the second stage. Is there a principled way to
   place group (a) on the same prediction-performance axis as group (c) methods that
   produce predictions end-to-end? Or should group (a) be evaluated differently, for
   example by how much phenotypic variance is captured by the retained components,
   rather than by held-out prediction R²?

2. **Bayesian/GP methods and small effective N:** The effective sample size is n=108
   lines (not ~1,700 flies), and between-line phenotype variance is the signal to
   capture. Gaussian process regression or Bayesian linear models with shrinkage priors
   might handle the high-p, small-n structure better than frequentist penalised
   regression, particularly if the spectral signal is weak. Are these worth adding to
   the comparison? If so, which formulation is most natural: GP regression directly on
   spectra, or Bayesian variable selection (e.g., spike-and-slab) on wavenumber
   coefficients?
