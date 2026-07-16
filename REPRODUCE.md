# Reproducing this project

> **Maintenance note:** this file must be updated whenever a new analysis
> script is added or a new phenotype is run through the pipeline. Add the
> exact command, expected outputs, and approximate runtime as a new step
> at the end of the relevant section. This file should always reflect
> every command that has actually been run and committed to this
> repository, in the order it was run. (Section 7 has already been
> revised once, when `prepare_unckless_data.py` was extended to cover
> the diet-specific columns as well as pooled-diet: an example of this
> maintenance practice in action, not just a description of it.)

This document lists every command needed to reproduce the results in this
repository from a fresh clone, in the order the analyses were originally
run. All commands assume your working directory is the repository root.

Timing notes are approximate. Steps marked "(timed)" were measured on the
machine used to develop this project (Apple Silicon Mac, one job running
at a time); everything else is an estimate based on the size of the
computation and has not been benchmarked precisely. Actual times will vary
with hardware and system load.

---

## 1. Environment setup

### Python

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install pandas==3.0.3 numpy==2.4.6 scikit-learn==1.9.0 scipy==1.17.1 \
            matplotlib==3.11.0 seaborn==0.13.2 xgboost==3.2.0 openpyxl==3.1.5
```

Developed and tested with Python 3.14. `openpyxl` is only needed for the
Unckless raw-data conversion step (Section 6); everything else only needs
the other packages.

All commands below invoke the interpreter directly as `.venv/bin/python`,
so you do not need to keep the venv activated in your shell.

### R

```bash
Rscript -e 'install.packages(c("survival","coxme","car","dplyr","emmeans",
  "ggplot2","effectsize","ggh4x","knitr","ciTools","here","survminer","rms","see"))'
```

Developed and tested with R 4.5.1.

---

## 2. Baseline reproduction

### 2.1 R survival analysis

```bash
Rscript scripts/run_survival_analysis.R
```

Fits a parametric survival model (`rms::psm`) to
`Survival-data/DGRP-starvationresistance.csv` and extracts one EMMean per
DGRP line.

**Output:** `Emmeans.csv`, `sensitive_df_20pct_emmean.csv`,
`resistant_df_80pct_emmean.csv` (repo root).
**Time:** not timed; a parametric survival fit and `emmeans` call over 108
lines, expect well under a minute.

### 2.2 Baseline SVM reproduction

Requires the R step above to have produced
`sensitive_df_20pct_emmean.csv` and `resistant_df_80pct_emmean.csv`.

```bash
.venv/bin/python scripts/run_dgrp_baseline.py
```

Runs XGBoost feature selection followed by 5 classifiers
(LR, SVM, kNN, RF, XGBoost) with 20-fold `GridSearchCV` on
`FTIR-data/DGRPFTIR.dat`.

**Output:** `results/DGRP/DGRP_XGBoost_CV_values.csv`,
`results/DGRP/DGRP_XGBoost_WNS_list.csv`, plus confusion-matrix, spectrum,
tSNE, and accuracy-boxplot PDFs and supporting CSVs in `results/DGRP/`.
**Time:** not timed; GridSearchCV over 5 classifiers with 20-fold CV, expect
several minutes.

Then, to print the SVM reclassification rates cleanly from the saved CV
values:

```bash
.venv/bin/python scripts/check_svm_rates.py
```

Prints to console only, no file output. Expected: ~88.5% resistant,
~83.9% sensitive.

### 2.3 FTIR data loader self-check

```bash
.venv/bin/python scripts/ftir_loader.py
```

Loads all six `.dat` files, validates the wavenumber axis, and prints
per-file metadata summaries.

**Output:** console only, no files written.
**Time:** not timed; just file I/O and validation, expect a few seconds.

### 2.4 Morgante external validation

Requires `Emmeans.csv` from Section 2.1.

```bash
.venv/bin/python scripts/check_morgante_overlap.py
```

Checks line-ID overlap between our EMMeans and the Morgante et al. 2015
starvation resistance data, and computes the cross-lab correlation.

**Output:** `results/DGRP/emmeans_vs_morgante_correlation.pdf`.
**Time:** not timed; a single correlation and scatter plot over ~104 lines,
expect a few seconds.

---

## 3. Method comparison on line-mean spectra

### 3.1 PCA compression analysis

```bash
.venv/bin/python scripts/run_compression_analysis.py
```

Compares PCA-then-average against average-then-PCA collapse orders and
plots the explained-variance curve.

**Output:** `results/DGRP/pca_explained_variance.pdf`,
`results/DGRP/pca_coloured_by_emmean.pdf`,
`results/DGRP/pca_linemeans_coloured_by_emmean.pdf`.
**Time:** not timed; a handful of PCA fits on ≤108 or ≤1,772 rows, expect
under a minute.

### 3.2 PLS regression analysis

```bash
.venv/bin/python scripts/run_pls_analysis.py
```

LOO-CV over 108 DGRP line-mean spectra, sweeping `n_components`, plus a
PCA+Ridge comparison on the same folds.

**Output:** `results/DGRP/pls_component1_vs_emmean.pdf`,
`results/DGRP/pls_loading1_vs_wavenumber.pdf`.
**Time:** not timed; 108 LOO folds with a lightweight PLS fit per fold,
expect a few minutes.

### 3.3 Regularised regression (Ridge, LASSO, elastic net)

```bash
.venv/bin/python scripts/run_regularised_regression.py
```

LOO-CV over 108 DGRP line-mean spectra for all three regularised methods,
each with its own inner hyperparameter search.

**Output:** `results/DGRP/regularised_coefficients_vs_wavenumber.pdf`.
**Time:** not timed; elastic net and LASSO each run an inner alpha/l1-ratio
search per fold across 1,723 features, so this is one of the slower
line-mean steps, plausibly comparable in order of magnitude to the
DGRPool per-phenotype runs in Section 6 (tens of minutes), though over
fewer folds (108 vs typically ≤108 lines there too).

---

## 4. Per-fly evaluation and random forest

### 4.1 Per-fly GroupKFold pipeline

```bash
.venv/bin/python scripts/run_perfly_pipeline.py
```

Trains PLS, Ridge, LASSO, and elastic net on ~1,772 individual fly spectra
with `GroupKFold(10)` (line-stratified), plus an inner `GroupKFold(5)` for
PLS component selection. Averages per-fly predictions to the line level
for evaluation.

**Output:** `results/DGRP/perfly_metrics.csv`.
**Time:** not timed; this trains on individual fly spectra (~1,772 rows,
16x more than the line-mean steps) across 10 outer folds with nested CV for
PLS: likely the most computationally expensive step in the project.
Expect this to take longer than any single line-mean script above.

### 4.2 Random forest

```bash
.venv/bin/python scripts/run_random_forest.py
```

LOO-CV over 108 DGRP line-mean spectra with `GridSearchCV` over
`n_estimators` and `max_features` inside each fold.

**Output:** `results/DGRP/rf_feature_importance_vs_wavenumber.pdf`.
**Time:** not timed; 108 folds each with a small grid search over tree
hyperparameters, expect several minutes.

---

## 5. Fecundity cross-phenotype (dedicated script)

```bash
.venv/bin/python scripts/run_fecundity_enet.py
```

Elastic net LOO-CV predicting lifetime fecundity
(`phenotype-data/S18_LifeFecundity_mean.tsv`) from the same spectra and
hyperparameters as Section 3.3, restricted to the 96 lines with a fecundity
value.

**Output:** console only, no files written (the null result does not
produce an informative plot).
**Time:** not timed; same per-fold cost profile as Section 3.3 but on 96
lines instead of 108, expect a broadly similar order of magnitude.

Note: this predates the general-purpose runner in Section 6 and is not run
through it: it is a standalone script with the phenotype hardcoded.

---

## 6. General-purpose DGRPool phenotype runner

`scripts/run_dgrpool_phenotype.py` runs elastic net LOO-CV against any
DGRPool-format phenotype TSV (columns `DGRP`, `sex`, `value`). It takes
`--sex` (filters the phenotype file, default `F`) and `--spectral-sex`
(filters the FTIR spectra, default `F`, since `DGRPFTIR.dat` is
female-only). The two only need to differ for cross-sex comparisons; a
console warning prints automatically whenever they do.

Each run appends one row to `results/DGRP/dgrpool_phenotype_summary.csv`.

**Time (timed):** every phenotype below took roughly 30–90 minutes in this
environment. All are the same computation (elastic net LOO-CV with a
nested alpha/l1-ratio search over 1,723 wavenumbers), so runtime mostly
tracks the number of overlapping lines and system load rather than the
phenotype itself.

### 6.1 Smoke test (internal EMMeans, validates the script is correct)

```bash
.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/S00_EMMeans_starvation.tsv \
  --study "Internal PSM model" --phenotype "Starvation resistance (EMMeans)"
```

Should reproduce CV R² ≈ 0.673 (matching Section 3.3's elastic net result);
if it doesn't, the general-purpose script has a bug.

### 6.2 Morgante starvation resistance (cross-lab validation)

```bash
.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/S24_StarvationRes_summary_mean.tsv \
  --study "Morgante 2015" --phenotype "Starvation resistance"
```

### 6.3 Lifespan

```bash
.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/S_Lifespan_mean.tsv \
  --study "Ivanov 2015" --phenotype "Lifespan"
```

### 6.4 Chill coma recovery

```bash
.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/S24_ChillComaRec_mean.tsv \
  --study "Morgante 2015" --phenotype "Chill coma recovery"
```

### 6.5 Cuticle hydrocarbon n-C25

```bash
.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/S_CuticHC_nC25_mean.tsv \
  --study "Dembeck 2015" --phenotype "Cuticle HC n-C25"
```

Sections 6.1–6.5 all use the default `--sex F` / `--spectral-sex F`
(same-sex comparison); the flags can be omitted since F is the default
for both, but are shown here for clarity.

---

## 7. Unckless 2015 nutritional indices (male, cross-sex comparison)

### 7.1 Obtain the raw supplementary file

This project does not download this file automatically. Obtain
Table S2 from Unckless RL, Rottschaefer SM, Lazzaro BP (2015). "A
Genome-Wide Association Study for Nutritional Indices in *Drosophila*."
*G3: Genes|Genomes|Genetics*, 5(3), 417–425.
https://doi.org/10.1534/g3.114.016477, and place it at:

```
phenotype-data/raw/016477_tables2.xlsx
```

### 7.2 Convert to per-measure TSVs

```bash
.venv/bin/python scripts/prepare_unckless_data.py
```

Extracts three source columns per metabolic measure (glucose, glycerol,
glycogen, triglyceride, protein, mean weight): `_pooled`, `_high_glucose`,
and `_low_glucose` (the two single-diet conditions), for 18 output files
total. Drops rows with missing (`.`) values independently per column,
reformats DGRP line IDs to match our spectral convention, and hardcodes
`sex=M` (all Unckless measures were assayed in pools of 10 adult males per
line, per the paper's Materials and Methods).

**Output:** for each of the six measures, `Unckless_<Measure>_pooled.tsv`,
`Unckless_<Measure>_highglucose.tsv`, and `Unckless_<Measure>_lowglucose.tsv`
(18 files total, all in `phenotype-data/`). Row counts differ by diet
condition, not by measure, because a different set of DGRP lines failed
the assay under each diet: 145 valid rows for pooled (153 source rows,
minus 7 missing values and 1 duplicate line), 147 for high-glucose (minus
5 missing values and the same duplicate line), 150 for low-glucose (minus
2 missing values and the same duplicate line). The missing/duplicate lines
are dropped identically across all six measures within a given diet
condition.
**Time:** a few seconds.

### 7.3 Run each pooled-diet measure through the general-purpose runner

Because the FTIR spectra used here (`DGRPFTIR.dat`) are female-only, these
are necessarily cross-sex comparisons (male metabolic phenotype vs. female
spectra, matched by DGRP line/genotype, not by matched sex). `--sex M`
is required; `--spectral-sex` is left at its default (`F`), so the console
will print the cross-sex warning for every run below.

```bash
.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_Glucose_pooled.tsv \
  --sex M --study "Unckless 2015" --phenotype "Glucose (male, pooled diet)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_Glycerol_pooled.tsv \
  --sex M --study "Unckless 2015" --phenotype "Glycerol (male, pooled diet)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_Glycogen_pooled.tsv \
  --sex M --study "Unckless 2015" --phenotype "Glycogen (male, pooled diet)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_Triglyceride_pooled.tsv \
  --sex M --study "Unckless 2015" --phenotype "Triglyceride (male, pooled diet)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_Protein_pooled.tsv \
  --sex M --study "Unckless 2015" --phenotype "Protein (male, pooled diet)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_MeanWeight_pooled.tsv \
  --sex M --study "Unckless 2015" --phenotype "MeanWeight (male, pooled diet)"
```

Each of the six overlaps to 77 lines (108 female spectral lines ∩ 145 male
phenotype lines).

### 7.4 Run each diet-specific measure through the general-purpose runner

Same cross-sex setup as Section 7.3 (`--sex M`, `--spectral-sex` left at
its default `F`), but against the single-diet columns extracted in
Section 7.2. Run in the order below: triglyceride and glycerol first
(both diets), as the most theoretically informative given the pooled
results, then the remaining four measures, both diets each.

```bash
.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_Triglyceride_highglucose.tsv \
  --sex M --study "Unckless 2015" --phenotype "Triglyceride (male, high-glucose diet)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_Triglyceride_lowglucose.tsv \
  --sex M --study "Unckless 2015" --phenotype "Triglyceride (male, low-glucose diet)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_Glycerol_highglucose.tsv \
  --sex M --study "Unckless 2015" --phenotype "Glycerol (male, high-glucose diet)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_Glycerol_lowglucose.tsv \
  --sex M --study "Unckless 2015" --phenotype "Glycerol (male, low-glucose diet)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_Glucose_highglucose.tsv \
  --sex M --study "Unckless 2015" --phenotype "Glucose (male, high-glucose diet)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_Glucose_lowglucose.tsv \
  --sex M --study "Unckless 2015" --phenotype "Glucose (male, low-glucose diet)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_Glycogen_highglucose.tsv \
  --sex M --study "Unckless 2015" --phenotype "Glycogen (male, high-glucose diet)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_Glycogen_lowglucose.tsv \
  --sex M --study "Unckless 2015" --phenotype "Glycogen (male, low-glucose diet)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_Protein_highglucose.tsv \
  --sex M --study "Unckless 2015" --phenotype "Protein (male, high-glucose diet)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_Protein_lowglucose.tsv \
  --sex M --study "Unckless 2015" --phenotype "Protein (male, low-glucose diet)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_MeanWeight_highglucose.tsv \
  --sex M --study "Unckless 2015" --phenotype "MeanWeight (male, high-glucose diet)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Unckless_MeanWeight_lowglucose.tsv \
  --sex M --study "Unckless 2015" --phenotype "MeanWeight (male, low-glucose diet)"
```

Each of the six high-glucose files overlaps to 77 lines (108 female
spectral lines ∩ 147 male phenotype lines); each of the six low-glucose
files overlaps to 80 lines (108 ∩ 150). The overlap counts differ from
the pooled-diet run (77 lines) because a different set of DGRP lines is
present in each diet condition's file, not because of the measure.

---

## Where results accumulate

- `results/DGRP/dgrpool_phenotype_summary.csv`: one row per run through
  Section 6's general-purpose runner (Sections 6 and 7).
- `results/DGRP/perfly_metrics.csv`: per-fly pipeline metrics (Section 4.1).
- All plots land in `results/DGRP/`.
- `Emmeans.csv` and the sensitive/resistant CSVs land in the repo root
  (Section 2.1).
