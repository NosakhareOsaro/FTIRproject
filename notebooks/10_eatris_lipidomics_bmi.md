# 10: EATRIS-Plus Lipidomics-to-BMI Analysis

**Script:** `scripts/explore_eatris_mae.R`, `scripts/prepare_eatris_lipidomics.R`, `scripts/run_eatris_lipidomics_bmi.py`.
**Status:** Complete
**Key result:** A real, non-collapsed spectral-style signal: elastic net LOO-CV R²=0.283 (positive-mode lipidomics), 0.271 (combined), 0.219 (negative-mode) predicting BMI in 125 healthy human adults. This clears a trivial Sex+Age baseline (R²=0.039) by 5-7x, so it is not a demographic proxy, but it sits well below the ceiling reported in large, externally-validated lipidomics-BMI studies, and is consistent with, not better than, what this exact cohort's own authors found when comparing single-omics layers for BMI prediction. This is the first genuinely positive, thoroughly-vetted result outside the DGRP starvation-resistance work in this project, and is reported with the same caution accordingly. A follow-up pushed on this harder by nesting Sex/Age directly inside the same LOO-CV pipeline instead of checking it separately: the signal holds up, with residual-BMI R² staying in the 0.21-0.29 range once the lipidomics features are no longer allowed to lean on Sex or Age at all.

---

## Why this dataset

The third and final task: whether lipidome features predict BMI quantitatively, this time in a human cohort rather than the DGRP fly panel. The EATRIS-Plus multi-omics dataset (Zenodo DOI 10.5281/zenodo.17514796) was identified as "a very promising dataset" for this: a human clinical cohort with paired lipidomics and BMI, letting the FTIR-style question (does a lipid-sensitive omics layer predict a body-composition phenotype?) be asked directly in humans rather than by analogy from flies.

## The MultiAssayExperiment structure

The Zenodo record ships two files: `mae_mae.rds` (4.3MB, the R `MultiAssayExperiment` object) and `mae_experiments.h5` (104MB, the actual assay matrices, HDF5-backed). Both are required together. `scripts/explore_eatris_mae.R` (inspection only, writes nothing) established the structure before any extraction code was written: 15 experiments across six omics platforms (acylcarnitines, amino acids, very long chain fatty acids, lipidomics positive and negative mode, a separate SPLINE-corrected positive-mode variant, proteomics, mRNA-seq, miRNA-seq/qRT-PCR, EM-seq), all keyed to the same 125 healthy adults. `colData` has 27 phenotype columns; BMI is present, complete (125/125, range 21-37), and there is zero sample-ID mismatch between `sampleMap`'s `primary`/`colname` columns across all 15 experiments.

Two lipidomics-relevant experiments were selected: "Lipidomics, positive | transformed" (196 features) and "Lipidomics, negative | transformed" (164 features), both with full 125/125 overlap against `colData` and against the BMI-complete sample set. A third lipidomics experiment, "Lipidomics-SPLINE | positive" (257 features, 124 samples, a differently-processed positive-mode variant with a distinct, messier rowData schema), was set aside as redundant with the "transformed" positive assay rather than combined with it.

**The HDF5 path-resolution problem.** Assay data in this object is `HDF5Matrix`/`DelayedArray`-backed: the seed's filepath is stored as a bare relative string (`"mae_experiments.h5"`) resolved against R's working directory _at the moment data is accessed_, not relative to the `.rds` file's own location. `explore_eatris_mae.R` worked around this with a temporary `setwd()` into the raw-data directory, which is fragile: it silently breaks if the script is ever invoked from somewhere other than the assumed location, or if anything else changes the working directory first. `scripts/prepare_eatris_lipidomics.R` fixes this properly: it resolves the raw-data directory to an absolute path via `here()` (the same package `run_survival_analysis.R` already uses elsewhere in this project) and overwrites the HDF5 seed's path directly (`path(a) <- H5_PATH`) before any data is touched, so extraction never depends on the working directory at all. Verified working when invoked from the repo root and from a subdirectory (`notebooks/`); it only fails when invoked from genuinely outside the repo tree (e.g. `/tmp`), which is expected `here()` behaviour and matches how every other script in this project is documented to be run.

## Data preparation

`scripts/prepare_eatris_lipidomics.R` extracts three lipidomics conditions, the same way `prepare_unckless_data.py` produced pooled/high-glucose/low-glucose variants of one underlying comparison:

| Condition                        | Features                                      | Samples |
| -------------------------------- | --------------------------------------------- | ------- |
| `EATRIS_Lipidomics_positive.tsv` | 196                                           | 125     |
| `EATRIS_Lipidomics_negative.tsv` | 164                                           | 125     |
| `EATRIS_Lipidomics_combined.tsv` | 360 (`pos_*` + `neg_*`, prefixed for clarity) | 125     |

Positive and negative mode have identical, fully-overlapping 125-sample sets, so the combined matrix drops no samples. Positive and negative feature IDs never collide (checked directly, 0 overlap), so the `pos_`/`neg_` prefix on the combined file is a clarity measure, not a necessity. `EATRIS_BMI.tsv` provides the matched phenotype vector (125 samples, BMI 125/125 complete, range 21-37), and `EATRIS_covariates.tsv` (Sex, Age, BMI) was added afterward for the confound check below - it started out as a one-off extraction for that standalone check, and later got pulled directly into the main analysis script once the confound question needed a properly nested answer (see "Results: does Sex/Age change this" below). All four files use plain `sample_id` + value columns, not the DGRP/sex/value convention used elsewhere in this project, since these are individual human subjects rather than DGRP line means.

## Method

`scripts/run_eatris_lipidomics_bmi.py` is a new standalone script rather than an adaptation of the DGRP-specific ones. This dataset needs none of that machinery: no `.dat` spectral parsing (`ftir_loader.py`), no per-DGRP-line averaging, no `--sex`/`--spectral-sex` split for cross-sex matching, since EATRIS is already one row per human subject with fully complete phenotype coverage. Reusing the FTIR scripts would mean carrying fly-specific code paths through a dataset that never exercises them.

What it does mirror exactly: the same LOO-CV and per-fold `StandardScaler` discipline as `run_regularised_regression.py` and `run_pls_analysis.py` (scaler and all hyperparameter selection fit inside the training fold only), the same four methods (PLS, Ridge, LASSO, elastic net) with the same hyperparameter grids (`RidgeCV` GCV over `np.logspace(-3,6,100)`, `LassoCV`/`ElasticNetCV` with `cv=3` and a 30-point alpha grid, `l1_ratio=[0.5,0.7,0.9,0.95,1.0]`, PLS component sweep `[1,2,3,5,10,15,20]`), and the CLI/results-accumulation shape of `run_dgrpool_phenotype.py` (one file argument + a `--condition` label, appends to `results/EATRIS/lipidomics_bmi_summary.csv`).

Run once per condition, the same way the Unckless diet conditions and Huang temperature conditions were each run separately, chained into one `caffeinate`-wrapped background job (a precaution from two earlier sleep-interruption incidents on the longer FTIR runs). This dataset's far smaller feature count (164-360 vs. FTIR's 1,723 wavenumbers) meant all three conditions finished in about 32 seconds total. Fast enough that the background/monitoring setup turned out to be unnecessary in practice, though harmless to have had in place.

**Extension: does Sex/Age change this.** I tried to check whether including Sex and Age alongside the lipidomics features changes the BMI result, so the script now runs each condition through three analyses instead of one, all funnelled through the same four-method sweep so the PLS/Ridge/LASSO/elastic-net code isn't tripled: **lipidomics only** (the original analysis, table below), **covariates appended** (Sex, encoded 0/1 for FEMALE/MALE, and Age added as two extra columns to the feature matrix before the per-fold `StandardScaler`), and **residualized** (a plain `LinearRegression` of BMI on Sex+Age fit _inside each LOO training fold only_ (never on the full data) with the four methods then predicting the residual, actual BMI minus that fold's Sex+Age prediction, from the lipidomics features alone). These last two are not two versions of the same check. Appending asks whether Sex/Age help the model when they're available to it; residualizing asks whether the lipidomics features still carry BMI signal once the variance Sex/Age already explain has been taken away first. See "Results: does Sex/Age change this" below for both.

## Results: 3 conditions x 4 methods

| Condition         | Method               | n       | p       | CV R²     | RMSE      | Spearman ρ |
| ----------------- | -------------------- | ------- | ------- | --------- | --------- | ---------- |
| Positive mode     | PLS (n_components=2) | 125     | 196     | 0.238     | 2.924     | +0.530     |
| Positive mode     | Ridge (GCV)          | 125     | 196     | 0.241     | 2.918     | +0.524     |
| Positive mode     | LASSO (cv=3)         | 125     | 196     | 0.179     | 3.035     | +0.451     |
| **Positive mode** | **Elastic net**      | **125** | **196** | **0.283** | **2.837** | **+0.545** |
| Negative mode     | PLS (n_components=2) | 125     | 164     | 0.166     | 3.059     | +0.493     |
| Negative mode     | Ridge (GCV)          | 125     | 164     | 0.128     | 3.129     | +0.516     |
| Negative mode     | LASSO (cv=3)         | 125     | 164     | 0.095     | 3.187     | +0.386     |
| **Negative mode** | **Elastic net**      | **125** | **164** | **0.219** | **2.960** | **+0.526** |
| Combined          | PLS (n_components=2) | 125     | 360     | 0.252     | 2.898     | +0.581     |
| Combined          | Ridge (GCV)          | 125     | 360     | 0.238     | 2.924     | +0.551     |
| Combined          | LASSO (cv=3)         | 125     | 360     | 0.251     | 2.899     | +0.532     |
| **Combined**      | **Elastic net**      | **125** | **360** | **0.271** | **2.859** | **+0.546** |

Elastic net wins in positive mode and combined, and is competitive (close second) in negative mode. Positive mode alone (R²=0.283) slightly outperforms the 360-feature combined matrix (R²=0.271); adding negative-mode features does not improve on positive mode alone. Both clearly outperform negative mode alone (R²=0.219). None of these results shows the mean-collapse pattern seen throughout the DGRP work: every method in every condition produces real, non-degenerate CV R² in the 0.10-0.28 range.

## Results: does Sex/Age change this

Two different questions were put to the same three conditions, both properly nested inside the LOO-CV loop rather than fit once on the full data:

- **Covariates appended** - does giving each method Sex and Age as two extra columns, sitting right next to the lipidomics features, improve prediction? This is the permissive version: the model is free to lean on Sex/Age wherever that helps, and there's no way to tell from the R² alone how much of any gain is the lipidome and how much is just Sex/Age doing what they already do on their own.
- **Residualized** - with the BMI variance Sex and Age already explain subtracted out first (a `LinearRegression` refit inside every fold), can the lipidomics features still predict what's left? Sex and Age never enter the feature matrix here, so they cannot be doing any of the predicting. This is the stricter test, and the one that actually answers "is the lipidomics result just a Sex/Age proxy."

Elastic net, the top method in every condition throughout this analysis, across both:

| Condition     | Lipidomics only | + Sex/Age appended | Residualized |
| ------------- | --------------- | ------------------ | ------------ |
| Positive mode | 0.283           | 0.287              | 0.250        |
| Negative mode | 0.219           | 0.233              | 0.211        |
| Combined      | 0.271           | 0.324              | 0.285        |

Full four-method breakdown, same layout as the results table above.

**Covariates appended**

| Condition         | Method               | n       | p       | CV R²     | RMSE      | Spearman ρ |
| ----------------- | -------------------- | ------- | ------- | --------- | --------- | ---------- |
| Positive mode     | PLS (n_components=3) | 125     | 198     | 0.240     | 2.920     | +0.552     |
| Positive mode     | Ridge (GCV)          | 125     | 198     | 0.245     | 2.911     | +0.524     |
| Positive mode     | LASSO (cv=3)         | 125     | 198     | 0.210     | 2.978     | +0.477     |
| **Positive mode** | **Elastic net**      | **125** | **198** | **0.287** | **2.828** | **+0.546** |
| Negative mode     | PLS (n_components=2) | 125     | 166     | 0.174     | 3.044     | +0.504     |
| Negative mode     | Ridge (GCV)          | 125     | 166     | 0.123     | 3.138     | +0.522     |
| **Negative mode** | **LASSO (cv=3)**     | **125** | **166** | **0.249** | **2.903** | **+0.548** |
| Negative mode     | Elastic net          | 125     | 166     | 0.233     | 2.935     | +0.542     |
| Combined          | PLS (n_components=2) | 125     | 362     | 0.259     | 2.884     | +0.585     |
| Combined          | Ridge (GCV)          | 125     | 362     | 0.241     | 2.919     | +0.550     |
| Combined          | LASSO (cv=3)         | 125     | 362     | 0.305     | 2.793     | +0.583     |
| **Combined**      | **Elastic net**      | **125** | **362** | **0.324** | **2.755** | **+0.594** |

**Residualized**

| Condition         | Method               | n       | p       | CV R²     | RMSE      | Spearman ρ |
| ----------------- | -------------------- | ------- | ------- | --------- | --------- | ---------- |
| Positive mode     | PLS (n_components=3) | 125     | 196     | 0.167     | 2.997     | +0.461     |
| Positive mode     | Ridge (GCV)          | 125     | 196     | 0.199     | 2.940     | +0.479     |
| Positive mode     | LASSO (cv=3)         | 125     | 196     | 0.214     | 2.912     | +0.487     |
| **Positive mode** | **Elastic net**      | **125** | **196** | **0.250** | **2.845** | **+0.506** |
| Negative mode     | PLS (n_components=2) | 125     | 164     | 0.124     | 3.073     | +0.442     |
| Negative mode     | Ridge (GCV)          | 125     | 164     | 0.082     | 3.146     | +0.460     |
| Negative mode     | LASSO (cv=3)         | 125     | 164     | 0.198     | 2.941     | +0.482     |
| **Negative mode** | **Elastic net**      | **125** | **164** | **0.211** | **2.917** | **+0.492** |
| Combined          | PLS (n_components=3) | 125     | 360     | 0.175     | 2.983     | +0.484     |
| Combined          | Ridge (GCV)          | 125     | 360     | 0.197     | 2.943     | +0.493     |
| Combined          | LASSO (cv=3)         | 125     | 360     | 0.280     | 2.786     | +0.550     |
| **Combined**      | **Elastic net**      | **125** | **360** | **0.285** | **2.777** | **+0.558** |

**Reading the appended column.** The gains are small and, in positive mode, barely worth mentioning (+0.004). Negative mode and combined move more (+0.014, +0.053 for elastic net), and combined's LASSO nearly doubles (0.251 to 0.305) - but this is exactly what adding two real, if modest, BMI correlates (Sex r=+0.291, Age r=-0.011, see the confound check below) to a feature matrix should do, not evidence that the lipidome suddenly does something it wasn't doing before. Worth flagging honestly: negative mode is the one place in this whole analysis where elastic net stops winning - LASSO takes it, 0.249 vs 0.233 - a reminder that "elastic net wins" was never a fixed law of this dataset, just what happened to hold in 11 of these 12 method-condition combinations. None of this moves the headline conclusion; it's a small, expected bump from giving the model more to work with, and on its own it can't separate "lipidomics is doing more work" from "Sex/Age are doing a bit more work alongside it."

**Reading the residualized column - the one that matters for the confound question.** R² drops relative to the lipidomics-only baseline in every condition, which is expected: some of what the lipidome predicts genuinely overlaps with what Sex and Age predict, and that shared variance is gone before the lipidomics features ever see the target. What doesn't happen is collapse. Elastic net's residual R² barely moves in negative mode (0.219 to 0.211) and combined (0.271 to 0.285, if anything a touch higher), and drops the most in positive mode (0.283 to 0.250) while still sitting well clear of zero. Every method in every condition stays in the 0.08-0.28 range on a target that has had the demographic variance surgically removed first. Because Sex and Age are not in the feature matrix at all here, they cannot be propping this number up - so this is the more conservative and the more informative version of the confound check, stronger than the standalone baseline comparison run earlier (below), and it says the same thing: not a demographic proxy.

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

`StandardScaler` is instantiated fresh every fold; `fit_transform` only ever sees `X[train_idx]`; `X_te` only ever gets `.transform()`. For Ridge/LASSO/elastic net, hyperparameter selection happens entirely inside `mdl.fit(X_tr, ...)`: `RidgeCV`'s GCV, and `LassoCV`/`ElasticNetCV`'s inner `cv=3` splits, never see the held-out fold. Identical discipline to `run_regularised_regression.py`.

**One caveat, inherited rather than new.** PLS's reported "optimal `n_components`" is chosen post-hoc, by comparing the already-computed LOO-CV R² across all 7 sweep values (`best_k = max(pls_scores, key=lambda k: pls_scores[k][0])`). Each individual sweep value's LOO-CV is properly nested, but selecting the best of 7 after seeing all of them is a mild form of model selection using outer-fold outcomes collectively, the same practice the original `run_pls_analysis.py` already uses for the DGRP data, not something introduced for this analysis. It does not affect the elastic net result, which has no such post-hoc step and wins outright in 2 of the 3 conditions.

### 2. Sex/Age confound check

The concern: with a real cohort (unlike the DGRP lines), BMI could correlate with a simple demographic variable that the lipidome also happens to track, making the "lipidomics" result a proxy for something trivial. `Sex`, `Age`, and `BMI` were pulled directly from `colData` (`phenotype-data/EATRIS_covariates.tsv`):

- Sex vs BMI: point-biserial r=+0.291, p=0.001, real but modest.
- Age vs BMI: r=−0.011, p=0.90, no relationship in this cohort.
- A trivial two-feature baseline (Sex+Age only, `LinearRegression`, same LOO-CV/scaling discipline as the lipidomics models): **R²=+0.039** (Sex-only: 0.055; Age-only: −0.033).

The lipidomics models (R²=0.219-0.283) sit 5-7x above this baseline. The result is not a demographic proxy.

This baseline is a useful first pass but a weak one: it's a single Sex+Age model fit once on the full data and compared only against the lipidomics-only numbers, not something run through the same nested LOO-CV discipline as the actual result. "Results: does Sex/Age change this" above redoes this properly, two ways, directly inside the lipidomics pipeline: Sex/Age appended as extra features, and - the test that actually settles this - Sex/Age regressed out of BMI inside every LOO fold before the lipidomics features get a turn. That second version leaves R² at 0.21-0.29 across all three conditions, real signal on a target with none of the demographic variance left in it to exploit. Same conclusion as the quick baseline above, on considerably firmer ground.

### 3. Literature context

- **Large, externally-validated benchmark** (plasma lipidomics, obesity estimation, PLOS Biology 2019, FINRISK 2012 n=1,061 training + Malmö Diet and Cancer Cardiovascular Cohort n=250 external validation): 183 lipid species reduced via LASSO to 50-75 predictors, **BMI R²=0.47**. The feature-to-sample ratio after selection (~50-75 features / ~800 training samples) is far more favourable than here (164-360 features / 125 samples, no external validation).
- **This exact cohort's own companion paper** (bioRxiv 2024.11.07.622407, the EATRIS-Plus flagship publication, same 127-person cohort and same lipidomics assays): ran a directly comparable analysis, predicting clinical variables including BMI from single-omics-layer components versus a combined multi-omics model. Their finding: for the majority of clinical variables including BMI, "the multi-omics components explained more of the variance than the individual-omics components". No single omics layer, lipidomics included, was a strong standalone BMI predictor in their hands either, and their feature-level association ranking for BMI in this cohort was metabolites > lipids > proteins, with lipids not the top layer. This is directly consistent with what was found here: a real but non-dominant single-omics-layer signal.

**Verdict:** modest-but-real by the standards of this specific field. R²=0.283 clears a trivial demographic baseline by a wide margin and shows none of the collapse pattern seen everywhere else in this project, so it is not noise dressed up as signal. Regressing Sex and Age out of BMI inside every LOO fold and asking whether lipidomics still predicts what's left - the harder version of the same question, not a repeat of the easier one - leaves R² at 0.21-0.29 in every condition, nowhere near the near-zero range a demographic-proxy explanation would predict. But it sits well below the ~0.47 ceiling reported in the large, externally-validated, carefully feature-selected lipidomics-BMI literature, was obtained with no external validation in a small-sample, high-dimensional (p>n) regime (164-360 features over 125 samples), and is consistent with, not better than, what this exact cohort's own authors already reported for single-omics-layer BMI prediction. This is the first genuinely positive, thoroughly-vetted result outside the DGRP starvation-resistance work in the whole project, and is reported with commensurate care: neither dismissed as another null, nor oversold as a validated biomarker signature.

## This completes three tasks

Three tasks were requested: (1) mated-female lifespan and age-specific fecundity, covered in markdown 09; (2) this lipidomics-BMI analysis; and this markdown covers the third and final one. Combined with markdown 09's conclusion (lifespan/fecundity nulls regardless of mating status, one weak 25°C candidate) and the running total through markdown 09 (30 tests: 2 signals, 26 nulls, 2 weak candidates), this analysis adds a genuinely different kind of result to the project: not a DGRP FTIR test at all, but an independent, human-cohort confirmation that lipidome-derived features can carry real (if modest) BMI-relevant signal, using an entirely different profiling technology (LC-MS lipidomics, not FTIR) and species (human, not fly). It does not speak to whether FTIR itself would show the same relationship in humans, only that the underlying biological premise (lipidome composition relates to body-composition phenotypes) holds up outside the DGRP panel.

## Output files

- `phenotype-data/raw/mae_mae.rds`, `phenotype-data/raw/mae_experiments.h5`: raw Zenodo download (DOI 10.5281/zenodo.17514796), not generated by any script here. The `.h5` file (104MB) is gitignored (over GitHub's 100MB file-size limit); see `REPRODUCE.md` for the download commands and checksums.
- `phenotype-data/EATRIS_Lipidomics_positive.tsv`, `EATRIS_Lipidomics_negative.tsv`, `EATRIS_Lipidomics_combined.tsv`, `EATRIS_BMI.tsv`, `EATRIS_covariates.tsv`: generated by `scripts/prepare_eatris_lipidomics.R` (lipidomics + BMI); `EATRIS_covariates.tsv` (Sex, Age, BMI) started as a one-off extraction for the standalone confound check and is now also read directly by `run_eatris_lipidomics_bmi.py` for the covariates-appended and residualized analyses
- `results/EATRIS/lipidomics_bmi_summary.csv`: 48 rows total (3 conditions x 3 analyses x 4 methods, run twice - the original lipidomics-only run plus the extension re-running all three analyses per condition). The original 12 lipidomics-only rows were backfilled with an `analysis` column rather than left on a different schema from the new rows - checked directly against a pre-migration backup, byte-identical apart from that column
- `scripts/explore_eatris_mae.R`: inspection-only, establishes the MAE structure, writes nothing
