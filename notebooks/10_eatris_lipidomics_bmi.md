# 10: EATRIS-Plus Lipidomics-to-BMI Analysis

**Script:** `scripts/explore_eatris_mae.R`, `scripts/prepare_eatris_lipidomics.R`, `scripts/run_eatris_lipidomics_bmi.py`.
**Status:** Complete
**Key result:** A real, non-collapsed spectral-style signal — elastic net LOO-CV R²=0.283 (positive-mode lipidomics), 0.271 (combined), 0.219 (negative-mode) predicting BMI in 125 healthy human adults. This clears a trivial Sex+Age baseline (R²=0.039) by 5-7x, so it is not a demographic proxy, but it sits well below the ceiling reported in large, externally-validated lipidomics-BMI studies, and is consistent with — not better than — what this exact cohort's own authors found when comparing single-omics layers for BMI prediction. This is the first genuinely positive, thoroughly-vetted result outside the DGRP starvation-resistance work in this project, and is reported with the same caution accordingly.

---

## Why this dataset

The third and final task from Adam's most recent email: whether lipidome features predict BMI quantitatively, this time in a human cohort rather than the DGRP fly panel. Adam identified the EATRIS-Plus multi-omics dataset (Zenodo DOI 10.5281/zenodo.17514796) as "a very promising dataset" for this — a human clinical cohort with paired lipidomics and BMI, letting the FTIR-style question (does a lipid-sensitive omics layer predict a body-composition phenotype?) be asked directly in humans rather than by analogy from flies.

## The MultiAssayExperiment structure

The Zenodo record ships two files: `mae_mae.rds` (4.3MB, the R `MultiAssayExperiment` object) and `mae_experiments.h5` (104MB, the actual assay matrices, HDF5-backed). Both are required together. `scripts/explore_eatris_mae.R` (inspection only, writes nothing) established the structure before any extraction code was written: 15 experiments across six omics platforms (acylcarnitines, amino acids, very long chain fatty acids, lipidomics positive and negative mode, a separate SPLINE-corrected positive-mode variant, proteomics, mRNA-seq, miRNA-seq/qRT-PCR, EM-seq), all keyed to the same 125 healthy adults. `colData` has 27 phenotype columns; BMI is present, complete (125/125, range 21-37), and there is zero sample-ID mismatch between `sampleMap`'s `primary`/`colname` columns across all 15 experiments.

Two lipidomics-relevant experiments were selected: "Lipidomics, positive | transformed" (196 features) and "Lipidomics, negative | transformed" (164 features), both with full 125/125 overlap against `colData` and against the BMI-complete sample set. A third lipidomics experiment, "Lipidomics-SPLINE | positive" (257 features, 124 samples, a differently-processed positive-mode variant with a distinct, messier rowData schema), was set aside as redundant with the "transformed" positive assay rather than combined with it.

**The HDF5 path-resolution problem.** Assay data in this object is `HDF5Matrix`/`DelayedArray`-backed: the seed's filepath is stored as a bare relative string (`"mae_experiments.h5"`) resolved against R's working directory *at the moment data is accessed* — not relative to the `.rds` file's own location. `explore_eatris_mae.R` worked around this with a temporary `setwd()` into the raw-data directory, which is fragile: it silently breaks if the script is ever invoked from somewhere other than the assumed location, or if anything else changes the working directory first. `scripts/prepare_eatris_lipidomics.R` fixes this properly: it resolves the raw-data directory to an absolute path via `here()` (the same package `run_survival_analysis.R` already uses elsewhere in this project) and overwrites the HDF5 seed's path directly (`path(a) <- H5_PATH`) before any data is touched, so extraction never depends on the working directory at all. Verified working when invoked from the repo root and from a subdirectory (`notebooks/`); it only fails when invoked from genuinely outside the repo tree (e.g. `/tmp`), which is expected `here()` behaviour and matches how every other script in this project is documented to be run.

## Data preparation

`scripts/prepare_eatris_lipidomics.R` extracts three lipidomics conditions, the same way `prepare_unckless_data.py` produced pooled/high-glucose/low-glucose variants of one underlying comparison:

| Condition | Features | Samples |
| --- | --- | --- |
| `EATRIS_Lipidomics_positive.tsv` | 196 | 125 |
| `EATRIS_Lipidomics_negative.tsv` | 164 | 125 |
| `EATRIS_Lipidomics_combined.tsv` | 360 (`pos_*` + `neg_*`, prefixed for clarity) | 125 |

Positive and negative mode have identical, fully-overlapping 125-sample sets, so the combined matrix drops no samples. Positive and negative feature IDs never collide (checked directly, 0 overlap), so the `pos_`/`neg_` prefix on the combined file is a clarity measure, not a necessity. `EATRIS_BMI.tsv` provides the matched phenotype vector (125 samples, BMI 125/125 complete, range 21-37), and `EATRIS_covariates.tsv` (Sex, Age, BMI) was added afterward for the confound check below. All four files use plain `sample_id` + value columns, not the DGRP/sex/value convention used elsewhere in this project, since these are individual human subjects rather than DGRP line means.

## Method

`scripts/run_eatris_lipidomics_bmi.py` is a new standalone script rather than an adaptation of the DGRP-specific ones. This dataset needs none of that machinery: no `.dat` spectral parsing (`ftir_loader.py`), no per-DGRP-line averaging, no `--sex`/`--spectral-sex` split for cross-sex matching, since EATRIS is already one row per human subject with fully complete phenotype coverage. Reusing the FTIR scripts would mean carrying fly-specific code paths through a dataset that never exercises them.

What it does mirror exactly: the same LOO-CV and per-fold `StandardScaler` discipline as `run_regularised_regression.py` and `run_pls_analysis.py` (scaler and all hyperparameter selection fit inside the training fold only), the same four methods (PLS, Ridge, LASSO, elastic net) with the same hyperparameter grids (`RidgeCV` GCV over `np.logspace(-3,6,100)`, `LassoCV`/`ElasticNetCV` with `cv=3` and a 30-point alpha grid, `l1_ratio=[0.5,0.7,0.9,0.95,1.0]`, PLS component sweep `[1,2,3,5,10,15,20]`), and the CLI/results-accumulation shape of `run_dgrpool_phenotype.py` (one file argument + a `--condition` label, appends to `results/EATRIS/lipidomics_bmi_summary.csv`).

Run once per condition, the same way the Unckless diet conditions and Huang temperature conditions were each run separately, chained into one `caffeinate`-wrapped background job (a precaution from two earlier sleep-interruption incidents on the longer FTIR runs). This dataset's far smaller feature count (164-360 vs. FTIR's 1,723 wavenumbers) meant all three conditions finished in about 32 seconds total — fast enough that the background/monitoring setup turned out to be unnecessary in practice, though harmless to have had in place.

## Results: 3 conditions x 4 methods

| Condition | Method | n | p | CV R² | RMSE | Spearman ρ |
| --- | --- | --- | --- | --- | --- | --- |
| Positive mode | PLS (n_components=2) | 125 | 196 | 0.238 | 2.924 | +0.530 |
| Positive mode | Ridge (GCV) | 125 | 196 | 0.241 | 2.918 | +0.524 |
| Positive mode | LASSO (cv=3) | 125 | 196 | 0.179 | 3.035 | +0.451 |
| **Positive mode** | **Elastic net** | **125** | **196** | **0.283** | **2.837** | **+0.545** |
| Negative mode | PLS (n_components=2) | 125 | 164 | 0.166 | 3.059 | +0.493 |
| Negative mode | Ridge (GCV) | 125 | 164 | 0.128 | 3.129 | +0.516 |
| Negative mode | LASSO (cv=3) | 125 | 164 | 0.095 | 3.187 | +0.386 |
| **Negative mode** | **Elastic net** | **125** | **164** | **0.219** | **2.960** | **+0.526** |
| Combined | PLS (n_components=2) | 125 | 360 | 0.252 | 2.898 | +0.581 |
| Combined | Ridge (GCV) | 125 | 360 | 0.238 | 2.924 | +0.551 |
| Combined | LASSO (cv=3) | 125 | 360 | 0.251 | 2.899 | +0.532 |
| **Combined** | **Elastic net** | **125** | **360** | **0.271** | **2.859** | **+0.546** |

Elastic net wins in positive mode and combined, and is competitive (close second) in negative mode. Positive mode alone (R²=0.283) slightly outperforms the 360-feature combined matrix (R²=0.271) — adding negative-mode features does not improve on positive mode alone — and both clearly outperform negative mode alone (R²=0.219). None of these results shows the mean-collapse pattern seen throughout the DGRP work: every method in every condition produces real, non-degenerate CV R² in the 0.10-0.28 range.

## Verification checks

A positive, non-collapsed R² this far into a project full of nulls is exactly the situation most likely to be over-interpreted, so three checks were run before accepting it as a real result.

### 1. Per-fold fitting discipline

Confirmed directly from the code. The core LOO-CV helper:

```python
def loo_cv(X, y, model_fn):
    """LOO-CV with per-fold StandardScaling. model_fn returns a fresh estimator."""
    y_pred = np.zeros(len(y))
    for train_idx, test_idx in loo.split(X):
        sc = StandardScaler()
        X_tr = sc.fit_transform(X[train_idx])
        X_te = sc.transform(X[test_idx])
        mdl = model_fn()
        mdl.fit(X_tr, y[train_idx])
        y_pred[test_idx] = np.asarray(mdl.predict(X_te)).ravel()
    return y_pred
```

`StandardScaler` is instantiated fresh every fold; `fit_transform` only ever sees `X[train_idx]`; `X_te` only ever gets `.transform()`. For Ridge/LASSO/elastic net, hyperparameter selection happens entirely inside `mdl.fit(X_tr, ...)` — `RidgeCV`'s GCV, and `LassoCV`/`ElasticNetCV`'s inner `cv=3` splits, never see the held-out fold. Identical discipline to `run_regularised_regression.py`.

**One caveat, inherited rather than new.** PLS's reported "optimal `n_components`" is chosen post-hoc, by comparing the already-computed LOO-CV R² across all 7 sweep values (`best_k = max(pls_scores, key=lambda k: pls_scores[k][0])`). Each individual sweep value's LOO-CV is properly nested, but selecting the best of 7 after seeing all of them is a mild form of model selection using outer-fold outcomes collectively — the same practice the original `run_pls_analysis.py` already uses for the DGRP data, not something introduced for this analysis. It does not affect the elastic net result, which has no such post-hoc step and wins outright in 2 of the 3 conditions.

### 2. Sex/Age confound check

The concern: with a real cohort (unlike the DGRP lines), BMI could correlate with a simple demographic variable that the lipidome also happens to track, making the "lipidomics" result a proxy for something trivial. `Sex`, `Age`, and `BMI` were pulled directly from `colData` (`phenotype-data/EATRIS_covariates.tsv`):

- Sex vs BMI: point-biserial r=+0.291, p=0.001 — real, but modest.
- Age vs BMI: r=−0.011, p=0.90 — no relationship in this cohort.
- A trivial two-feature baseline (Sex+Age only, `LinearRegression`, same LOO-CV/scaling discipline as the lipidomics models): **R²=+0.039** (Sex-only: 0.055; Age-only: −0.033).

The lipidomics models (R²=0.219-0.283) sit 5-7x above this baseline. The result is not a demographic proxy.

### 3. Literature context

- **Large, externally-validated benchmark** (plasma lipidomics, obesity estimation, PLOS Biology 2019, FINRISK 2012 n=1,061 training + Malmö Diet and Cancer Cardiovascular Cohort n=250 external validation): 183 lipid species reduced via LASSO to 50-75 predictors, **BMI R²=0.47**. The feature-to-sample ratio after selection (~50-75 features / ~800 training samples) is far more favourable than here (164-360 features / 125 samples, no external validation).
- **This exact cohort's own companion paper** (bioRxiv 2024.11.07.622407, the EATRIS-Plus flagship publication, same 127-person cohort and same lipidomics assays): ran a directly comparable analysis, predicting clinical variables including BMI from single-omics-layer components versus a combined multi-omics model. Their finding: for the majority of clinical variables including BMI, "the multi-omics components explained more of the variance than the individual-omics components" — no single omics layer, lipidomics included, was a strong standalone BMI predictor in their hands either, and their feature-level association ranking for BMI in this cohort was metabolites > lipids > proteins, with lipids not the top layer. This is directly consistent with what was found here: a real but non-dominant single-omics-layer signal.

**Verdict:** modest-but-real by the standards of this specific field. R²=0.283 clears a trivial demographic baseline by a wide margin and shows none of the collapse pattern seen everywhere else in this project, so it is not noise dressed up as signal. But it sits well below the ~0.47 ceiling reported in the large, externally-validated, carefully feature-selected lipidomics-BMI literature, was obtained with no external validation in a small-sample, high-dimensional (p>n) regime (164-360 features over 125 samples), and is consistent with — not better than — what this exact cohort's own authors already reported for single-omics-layer BMI prediction. This is the first genuinely positive, thoroughly-vetted result outside the DGRP starvation-resistance work in the whole project, and is reported with commensurate care: neither dismissed as another null, nor oversold as a validated biomarker signature.

## This completes Adam's three-task email

Three tasks were requested: (1) mated-female lifespan and age-specific fecundity, covered in markdown 09; (2) this lipidomics-BMI analysis; and this markdown covers the third and final one. Combined with markdown 09's conclusion (lifespan/fecundity nulls regardless of mating status, one weak 25°C candidate) and the running total through markdown 09 (30 tests: 2 signals, 26 nulls, 2 weak candidates), this analysis adds a genuinely different kind of result to the project: not a DGRP FTIR test at all, but an independent, human-cohort confirmation that lipidome-derived features can carry real (if modest) BMI-relevant signal, using an entirely different profiling technology (LC-MS lipidomics, not FTIR) and species (human, not fly). It does not speak to whether FTIR itself would show the same relationship in humans, only that the underlying biological premise — lipidome composition relates to body-composition phenotypes — holds up outside the DGRP panel.

## Output files

- `phenotype-data/raw/mae_mae.rds`, `phenotype-data/raw/mae_experiments.h5`: raw Zenodo download (DOI 10.5281/zenodo.17514796), not generated by any script here. The `.h5` file (104MB) is gitignored (over GitHub's 100MB file-size limit); see `REPRODUCE.md` for the download commands and checksums.
- `phenotype-data/EATRIS_Lipidomics_positive.tsv`, `EATRIS_Lipidomics_negative.tsv`, `EATRIS_Lipidomics_combined.tsv`, `EATRIS_BMI.tsv`, `EATRIS_covariates.tsv`: generated by `scripts/prepare_eatris_lipidomics.R` (lipidomics + BMI) and a one-off extraction for the Sex/Age confound check (covariates)
- `results/EATRIS/lipidomics_bmi_summary.csv`: appended with 12 rows (3 conditions x 4 methods) by `scripts/run_eatris_lipidomics_bmi.py`
- `scripts/explore_eatris_mae.R`: inspection-only, establishes the MAE structure, writes nothing
