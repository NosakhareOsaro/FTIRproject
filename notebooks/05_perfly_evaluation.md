# 05: Per-Fly Evaluation with Line-Stratified GroupKFold CV

**Script:** `scripts/run_perfly_pipeline.py`
**Status:** Complete
**Key result:** PLS wins the per-fly setting (line-level R² = 0.534, ρ = +0.794), despite losing the line-mean setting. The ranking of methods inverts between the two evaluation settings.

---

## What this step does and why

All analyses in markdowns 02, 03, and 04 used **line-mean spectra**: the 108 DGRP lines were each represented by a single average spectrum computed from ~16 individual fly spectra. This is a clean and principled approach, but it discards within-line information and reduces the sample size from 1,772 flies to 108 lines.

A more challenging and arguably more realistic evaluation is to train models on **individual fly spectra** and ask whether the predictions, when averaged to the line level, can recover the line-level starvation EMMeans. This is harder because:

1. The model has to learn the between-line signal from noisy within-line data
2. The cross-validation split must be done carefully to avoid data leakage

This step runs four methods (PLS, Ridge, LASSO, elastic net) on individual fly spectra with line-stratified cross-validation, then averages predictions to the line level for evaluation.

---

## The data leakage problem and how it is solved

In the line-mean pipeline, one line is held out at a time (LOO-CV). There is no leakage because each line appears in either training or test, never both.

With per-fly spectra, the same line-level discipline must be maintained. If flies from the same DGRP line appear in both training and test folds, two problems arise:

1. **Label leakage:** all flies from the same line share the same EMMean label. A fly in the test fold has the same label as its siblings in the training fold, so the model can partially infer the test label by recognising the line.
2. **Spectral leakage:** flies from the same inbred line have very similar spectra. A model trained on some flies from DGRP100 will generalise trivially to other flies from DGRP100, inflating the apparent accuracy.

**The solution is GroupKFold cross-validation with DGRP line as the group.** `GroupKFold(n_splits=10)` splits the 1,772 flies into 10 folds such that all flies from a given DGRP line always land in the same fold, either all in training or all in test, never split across folds. This is the correct design given the destructive-assay constraint.

---

## Evaluation design

**Outer CV:** `GroupKFold(n_splits=10)` with DGRP line as the group. 10 folds, each with approximately 10-11 test lines and 97-98 training lines (~1,590-1,600 training flies, ~170-180 test flies).

**Inside each outer fold:**
1. `StandardScaler` fitted on the training flies only (never the test flies)
2. The model is trained on the scaled training fly spectra, with each fly's EMMean used as the training target (the same EMMean value repeated for all ~16 flies from the same line)
3. The model predicts a score for each test fly

**Line-level aggregation:**
After all 10 folds, predictions for all test flies are collected. For each DGRP line in the test set, the per-fly predictions are averaged to produce one line-level prediction. This line-level prediction is then compared to the true EMMean. All performance metrics (R², RMSE, Spearman ρ) are computed at the line level.

**For PLS only, inner GroupKFold for component selection:**
PLS requires selecting the number of components. Rather than fixing this in advance, an inner `GroupKFold(n_splits=5)` is run on the outer training flies to select the best n_components from {1, 2, 3, 5, 10}. This inner CV also uses line-stratified groups: even within the training fold, flies from the same line are kept together. This avoids a secondary leakage at the inner level.

---

## Results

### Line-level R² after per-fly averaging

| Method | Line-level R² | Line RMSE | Spearman ρ | Pearson r |
|---|---|---|---|---|
| Ridge | 0.515 | 0.5130 | +0.767 | +0.764 |
| LASSO | 0.517 | 0.5121 | +0.781 | +0.768 |
| Elastic net | 0.518 | 0.5113 | +0.783 | +0.771 |
| **PLS** | **0.534** | **0.5029** | **+0.794** | **+0.802** |

### Per-fly R² (before line-level averaging)

| Method | Per-fly R² |
|---|---|
| PLS | 0.360 |
| Ridge | 0.309 |
| LASSO | 0.327 |
| Elastic net | 0.332 |

---

## The ranking flip

In the line-mean setting (markdown 04), elastic net was the best method (R² = 0.673) and PLS was third (R² = 0.623). In the per-fly setting, PLS is the best method (line R² = 0.534) and elastic net, LASSO, and Ridge are essentially tied behind it (R² = 0.515-0.518).

**Why does the ranking flip?**

In the line-mean setting, within-line noise has already been removed by averaging. The models are trained on clean, noise-free spectra and the only variation is between-line genetic variation. Sparse linear methods (elastic net, LASSO) are well suited to this because they can identify the small number of spectral regions (lipid bands) that carry the between-line signal and ignore the rest.

In the per-fly setting, within-line noise is present in the training data. Each fly's spectrum contains both the genuine between-line genetic signal and random within-fly measurement noise. Sparse methods that rely on sharp coefficient selection become less stable: they pick up on individual-fly spectral quirks that do not generalise to the line average. PLS's latent variable structure, by finding directions of maximum covariance with the phenotype across many flies simultaneously, is more robust to this noise.

In short: PLS is better at filtering noise during training; sparse linear methods are better when noise has already been removed.

---

## Comparison between evaluation settings

| Method | Line-mean LOO-CV R² | Per-fly GroupKFold line R² | Drop |
|---|---|---|---|
| PLS | 0.623 | 0.534 | -0.089 |
| Ridge | 0.635 | 0.515 | -0.120 |
| LASSO | 0.669 | 0.517 | -0.152 |
| Elastic net | 0.673 | 0.518 | -0.155 |

All methods perform worse in the per-fly setting than the line-mean setting. This is expected: the per-fly setting is harder, both because within-line noise is present during training and because GroupKFold(10) is a less generous evaluation than LOO. The drop is not a failure of the models; it quantifies the cost of within-line spectral noise and the benefit of pre-averaging.

**The line-mean setting represents the performance ceiling**, the best achievable when noise is removed. **The per-fly setting represents the realistic bound**, what a model can achieve when trained and evaluated on the same terms as a real deployment scenario (where you would have individual fly spectra, not pre-computed line means).

---

## Per-fly R² vs line R²: an important distinction

The per-fly R² values (0.31-0.36) are much lower than the line R² values (0.51-0.53) for all methods. This is expected and not a problem.

Per-fly R² is computed before any averaging: it measures how well the model predicts an individual fly's starvation EMMean from that fly's spectrum alone. Because within-line spectral variation is large relative to between-line variation, individual predictions are noisy.

Line R² is computed after averaging predictions within each test line: it measures how well the averaged predictions recover the line EMMean. Averaging 10-16 fly predictions per line cancels much of the within-fly noise, which is why line R² is substantially higher.

Line R² is the correct evaluation metric for this project. The goal is to predict DGRP line phenotype from FTIR spectra, not to predict the phenotype of an individual fly.

---

## Technical details

- `sklearn.model_selection.GroupKFold` was used for both outer and inner CV
- Predictions were collected across all 10 outer folds and aggregated with `pandas.groupby("line_id").mean()` before computing line-level metrics
- The `y_train` target for each fly is its line's EMMean, repeated ~16 times. This is the only possible target given the destructive assay design, since individual fly phenotypes are not available
- `sklearn.linear_model.LassoCV` and `ElasticNetCV` with `alphas=30` (integer count in sklearn 1.9.0) and `tol=0.01` were used for inner α selection

---

## Output files

- `results/DGRP/perfly_metrics.csv`: full per-method metrics table (line R², RMSE, Spearman ρ, Pearson r, per-fly R²)