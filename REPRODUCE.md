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
Unckless raw-data conversion step (Section 7); everything else only needs
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

### 5.1 Obtain the raw DGRPool phenotype files

This section and Section 6 both use phenotype files downloaded by hand
from DGRPool (dgrpool.epfl.ch) - these are browser downloads from the
site's phenotype browser, not a single scriptable URL the way the
Unckless/EATRIS sources are, so there's no `curl` command to give here.
Place the files at:

```
phenotype-data/S18_LifeFecundity_mean.tsv      (Study 18: dgrpool.epfl.ch/studies/18)
phenotype-data/S24_StarvationRes_summary_mean.tsv  (Study 24, Morgante et al. 2015: dgrpool.epfl.ch/phenotypes/2798)
phenotype-data/S_Lifespan_mean.tsv             (Ivanov et al. 2015; study number wasn't noted at download time)
phenotype-data/S24_ChillComaRec_mean.tsv       (Study 24, Morgante et al. 2015)
phenotype-data/S_CuticHC_nC25_mean.tsv         (Dembeck et al. 2015; study number wasn't noted at download time)
```

The exact phenotype-browser URL wasn't recorded for the last three at the
time they were downloaded (2026-07-03) - only the study/paper. If these
ever need re-downloading, `phenotype-data/README.md`'s per-file notes are
the fullest record of what's in each: line counts, sex composition, and
which of our 108 FTIR lines are present.

`phenotype-data/S00_EMMeans_starvation.tsv` is different from the other
four: it isn't a DGRPool download at all, it's our own `Emmeans.csv`
(Section 2.1) reformatted into the same DGRP/sex/value shape. Regenerate
it with:

```bash
.venv/bin/python scripts/prepare_smoke_test_phenotype.py
```

**Output:** `phenotype-data/S00_EMMeans_starvation.tsv` (108 rows).
**Time:** instant, it's a column rename and a hardcoded `sex="F"`.

### 5.2 Run the fecundity analysis

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

The RESULTS block also prints the prediction-SD/true-SD ratio on every
run (not just when it drops below the 0.2 collapse threshold, which was
the original behaviour). This lets you distinguish a genuine null where
predictions collapse to the training mean from a null where predictions
vary but simply aren't accurate, without needing a separate diagnostic
script.

Each run appends one row to `results/DGRP/dgrpool_phenotype_summary.csv`.

**Time (timed):** every phenotype below took roughly 30–90 minutes in this
environment. All are the same computation (elastic net LOO-CV with a
nested alpha/l1-ratio search over 1,723 wavenumbers), so runtime mostly
tracks the number of overlapping lines and system load rather than the
phenotype itself.

Sections 6.1-6.5 all need the phenotype files obtained/generated in
Section 5.1 above.

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
Genome-Wide Association Study for Nutritional Indices in _Drosophila_."
_G3: Genes|Genomes|Genetics_, 5(3), 417–425.
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

### 7.4 Run each diet-specific measure through the general-purpose runner (verified)

All 12 commands below have been run to completion; each produced exactly
one row in `results/DGRP/dgrpool_phenotype_summary.csv`, in the same
order as documented here, with no errors. See
`phenotype-data/README.md` for the full 18-row results table (6 measures
× pooled/high-glucose/low-glucose) and the diagnostic used to check
whether the one non-null result (Protein, low-glucose) is a genuine
signal or a collapse artifact.

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

## 8. Lifespan/fecundity extension (Durham 2014 + Huang 2020, same-sex comparison)

Follow-up requested by Adam (email, 2026-08-01): test whether FTIR
predicts lifespan in mated females (§6.3's Ivanov result used virgin
females), and test age-specific fecundity at the youngest and oldest
measured timepoints.

### 8.1 Obtain the raw source files

This project does not download these files automatically. Place them at:

```
phenotype-data/raw/S12_lsm_Lifespan_original.tsv
phenotype-data/raw/S12_lsm_Week1_Fecundity_original.tsv
phenotype-data/raw/S12_lsm_Week7_Fecundity_original.tsv
phenotype-data/raw/S40_Lifespan_18C_original.tsv
phenotype-data/raw/S40_Lifespan_25C_original.tsv
phenotype-data/raw/S40_Lifespan_28C_original.tsv
```

The three `S12_lsm_*` files are Durham, Magwire, Stone & Leips (2014,
DGRPool study 12), mated female lifespan and age-specific fecundity
(least-squares means, not raw values). The three `S40_*` files are Huang
et al. (2020, DGRPool study 40), individual-fly mated female lifespan at
18°C, 25°C, and 28°C.

### 8.2 Convert to DGRPool-format TSVs

```bash
.venv/bin/python scripts/prepare_lifespan_fecundity_data.py
```

The three Durham files are already in clean DGRP/sex/value format and are
copied through unchanged to `Durham_Lifespan_mated.tsv`,
`Durham_Week1_Fecundity.tsv`, and `Durham_Week7_Fecundity.tsv` (189, 189,
and 135 rows respectively). One negative value in the Week 7 file
(DGRP_042, −0.165) is expected: these are least-squares means corrected
for body size and block, not raw egg counts.

The three Huang files are individual fly-level raw data (23,000+ rows
each, both sexes). The script filters to `sex == 'F'`, groups by DGRP,
and averages to produce `Huang_Lifespan_18C_female.tsv`,
`Huang_Lifespan_25C_female.tsv`, and `Huang_Lifespan_28C_female.tsv`
(183, 186, and 177 female lines respectively).

Prints a six-row summary table (filename, n, min, max, mean) for review
before proceeding.
**Time:** a few seconds.

### 8.3 Run each phenotype through the general-purpose runner

All six use default `--sex F`/`--spectral-sex F` (same-sex comparison; no
`--spectral-sex` flag needed, unlike Section 7's cross-sex Unckless work).

```bash
.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Durham_Lifespan_mated.tsv \
  --study "Durham 2014" --phenotype "Lifespan (mated, Durham 2014)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Durham_Week1_Fecundity.tsv \
  --study "Durham 2014" --phenotype "Fecundity Week 1 (Durham 2014)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Durham_Week7_Fecundity.tsv \
  --study "Durham 2014" --phenotype "Fecundity Week 7 (Durham 2014)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Huang_Lifespan_18C_female.tsv \
  --study "Huang 2020" --phenotype "Lifespan 18C (mated, Huang 2020)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Huang_Lifespan_25C_female.tsv \
  --study "Huang 2020" --phenotype "Lifespan 25C (mated, Huang 2020)"

.venv/bin/python scripts/run_dgrpool_phenotype.py \
  phenotype-data/Huang_Lifespan_28C_female.tsv \
  --study "Huang 2020" --phenotype "Lifespan 28C (mated, Huang 2020)"
```

Overlap with our 108 FTIR lines: 96 (Durham lifespan and week 1
fecundity), 76 (Durham week 7 fecundity), 99 (Huang 18°C), 102 (Huang
25°C), 99 (Huang 28°C). See `phenotype-data/README.md` and
`PROJECT_NOTES.md` §7e for full results and interpretation: 5 of 6 null,
one weak unconfirmed candidate (Huang 25°C, R²=+0.056).

**Timing note (observed on this machine):** each of these six runs took
noticeably longer than the Unckless runs in Section 7 (roughly 40–100
minutes wall-clock each), because these phenotype files have more
overlapping DGRP lines (76–102 vs. 77–91 for Unckless), and LOO-CV cost
scales with the number of lines. If running overnight or unattended,
wrap each command in `caffeinate -i` (macOS) or the equivalent to prevent
system sleep from pausing the computation mid-run.

---

## 9. EATRIS-Plus lipidomics-to-BMI (human cohort, R + Python)

Third task from Adam's most recent email: whether lipidome features predict
BMI quantitatively, this time in the EATRIS-Plus human clinical cohort
(125 healthy adults) rather than the DGRP fly panel.

### 9.1 R packages

```bash
Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")
BiocManager::install(c("MultiAssayExperiment","HDF5Array","rhdf5","SummarizedExperiment"), update=FALSE, ask=FALSE)'
```

`here` is also required (already installed for Section 2.1's survival analysis).

### 9.2 Obtain the raw MultiAssayExperiment object

This project does not download these files automatically. Obtain both files
from Zenodo (DOI 10.5281/zenodo.17514796) and place them at:

```
phenotype-data/raw/mae_mae.rds
phenotype-data/raw/mae_experiments.h5
```

```bash
cd phenotype-data/raw
curl -sL -o mae_mae.rds "https://zenodo.org/records/17514796/files/mae_mae.rds?download=1"
curl -sL -o mae_experiments.h5 "https://zenodo.org/records/17514796/files/mae_experiments.h5?download=1"
```

Verify checksums:

```bash
md5 mae_mae.rds mae_experiments.h5
# mae_mae.rds:        329e82b9341071e2a91c20e537903505
# mae_experiments.h5: 9938d9731fe58ca67326109a58a47a82
```

`mae_experiments.h5` (104MB) is gitignored, since it's over GitHub's
100MB file-size limit, so it has to be downloaded separately even after
cloning this repo. `mae_mae.rds` (4.3MB) is small enough to be tracked
normally.

### 9.3 Inspect the structure (no files written)

```bash
Rscript scripts/explore_eatris_mae.R
```

Reports all 15 experiments and dimensions, the full `colData` column list
and BMI completeness, and sample ID overlap between the lipidomics assays
and the phenotype table. See `notebooks/10_eatris_lipidomics_bmi.md` for
the full findings.

### 9.4 Extract the lipidomics matrices and BMI phenotype

```bash
Rscript scripts/prepare_eatris_lipidomics.R
```

Extracts positive-mode (196 features), negative-mode (164 features), and a
combined (360 features, `pos_`/`neg_` prefixed) lipidomics matrix, the BMI
phenotype vector, and a small Sex/Age/BMI covariates file, all matched to
the same 125 samples. Resolves the HDF5-backed assay data via an absolute
path computed with `here()` rather than a working-directory-dependent
`setwd()`, so it can be run from anywhere inside the repo. Writes:

```
phenotype-data/EATRIS_Lipidomics_positive.tsv
phenotype-data/EATRIS_Lipidomics_negative.tsv
phenotype-data/EATRIS_Lipidomics_combined.tsv
phenotype-data/EATRIS_BMI.tsv
phenotype-data/EATRIS_covariates.tsv
```

`EATRIS_covariates.tsv` (Sex, Age, BMI) isn't used by the regression
pipeline itself - it only feeds the demographic-confound check in
`notebooks/10_eatris_lipidomics_bmi.md` (is the lipidomics signal just
recovering Sex or Age?). It was originally pulled out with a quick
interactive `Rscript -e` one-liner and only folded into this script
afterward, once it was clear the confound check was worth keeping around
rather than a throwaway sanity check.

Prints a summary (dimensions, sample counts, NA counts) for all five
outputs. **Time:** a few seconds.

### 9.5 Run the method comparison on all three conditions

```bash
.venv/bin/python scripts/run_eatris_lipidomics_bmi.py \
  phenotype-data/EATRIS_Lipidomics_positive.tsv --condition "Positive mode"

.venv/bin/python scripts/run_eatris_lipidomics_bmi.py \
  phenotype-data/EATRIS_Lipidomics_negative.tsv --condition "Negative mode"

.venv/bin/python scripts/run_eatris_lipidomics_bmi.py \
  phenotype-data/EATRIS_Lipidomics_combined.tsv --condition "Combined"
```

Each runs PLS, Ridge, LASSO, and elastic net with LOO-CV over all 125
samples, and appends 4 rows to
`results/EATRIS/lipidomics_bmi_summary.csv`. **Time (observed on this
machine):** all three conditions complete in well under a minute combined,
far faster than the DGRP FTIR runs, since these matrices have 164-360
features versus FTIR's 1,723 wavenumbers.

Full results, the per-fold discipline verification, the Sex/Age confound
check, and literature context are in
`notebooks/10_eatris_lipidomics_bmi.md`.

---

## Where results accumulate

- `results/DGRP/dgrpool_phenotype_summary.csv`: one row per run through
  Section 6's general-purpose runner (Sections 6, 7, and 8).
- `results/DGRP/perfly_metrics.csv`: per-fly pipeline metrics (Section 4.1).
- `results/EATRIS/lipidomics_bmi_summary.csv`: one row per method per
  condition run through Section 9's runner.
- All plots land in `results/DGRP/`.
- `Emmeans.csv` and the sensitive/resistant CSVs land in the repo root
  (Section 2.1).
