# 08: Unckless Metabolic Cross-Phenotype Analysis

**Script:** `scripts/prepare_unckless_data.py`, `scripts/run_dgrpool_phenotype.py`
**Status:** Complete
**Key result:** 17 of 18 elastic net LOO-CV tests against Unckless et al. 2015 metabolic measures show no spectral signal, including glycogen, the strongest published correlate of starvation resistance in their own data (r=0.307). One weak result (protein, low-glucose diet, R²=0.066) is not a mean-collapse artefact, but it is uncorrected for multiple comparisons across 18 tests and should not be read as a finding. Combined with the five phenotype targets from markdown 07, the FTIR signal looks specific to this lab's own starvation assay rather than to metabolic status in general.

---

## What Adam asked for

At the meeting on 15 July 2026, Adam set out the next task directly. Markdown 07 had shown that our FTIR spectra predict starvation resistance measured in our own lab, but not the same phenotype measured by another lab (Morgante 2015), and not any of the other traits tested so far (lifespan, chill coma recovery, cuticle hydrocarbons, fecundity). The open question was whether that pattern means the FTIR chemotype only tracks our specific starvation assay, or whether it has no predictive power for anything outside our own lab's data at all. As Adam put it: "the key thing would be to assess whether the FTIR spectra from our lab have any predictive power for other traits at all." He sent over a dataset from Unckless et al. 2015 to test this directly.

## Why the Unckless dataset

Unckless, Rottschaefer and Lazzaro (2015) measured six metabolic indices across the DGRP panel: glucose, glycerol, glycogen, triglyceride, protein, and wet weight, each under two diet conditions (high-glucose and low-glucose) plus a pooled measure. This is a good test case for two reasons. First, it targets metabolic status directly rather than a downstream behavioural outcome like starvation survival, so it separates two possible explanations for the starvation result: does FTIR detect lipid and energy reserves generally, or does it only track the specific assay we run in this lab? Second, the paper reports its own correlations between these metabolic measures and starvation resistance (measured by Mackay et al. 2012), which gives an independent, published benchmark to compare against our own null results.

## The male-only data problem

Reading the paper's Materials and Methods showed that all six measures were assayed in pools of 10 adult males per DGRP line. Our FTIR spectral data (`DGRPFTIR.dat`) is female-only. There is no way around this with the data available: any comparison against the Unckless measures is necessarily a cross-sex comparison, matching DGRP lines by genotype rather than by sex.

This is a defensible design given the DGRP's use as an inbred panel: each line is effectively a fixed genotype, and a genotype-level correlation between a male metabolic measure and a female spectral signature is a real biological question, just not the same question as a same-sex comparison would be. It does mean any result has to be read with that caveat attached, which is why it is stated explicitly in every run below.

## Fixing the pipeline for cross-sex comparisons

`scripts/run_dgrpool_phenotype.py` originally took a single `--sex` flag and used it to filter both the phenotype file and the spectral data. That works fine when both are the same sex, which every phenotype tested up to markdown 07 was. It breaks for the Unckless data: passing `--sex M` would have filtered `DGRPFTIR.dat` down to its male spectra too, and since that file contains none, the spectral side would come out empty and the run would abort with zero overlapping lines before it ever compared anything.

The fix splits this into two flags. `--sex` filters the phenotype file. A new `--spectral-sex` flag, defaulting to `F`, filters the spectral file separately. Every phenotype tested before this markdown used `F` for both, so nothing about those results changes. For the Unckless runs, `--sex M --spectral-sex F` (the second flag left at its default) does the cross-sex lookup correctly. The script also now prints an explicit warning whenever the two flags differ, so a cross-sex run can never happen silently.

## Pipeline audit

Finding one real bug in a script that had already produced five published-looking results was reason enough to check the rest of the pipeline before trusting any new numbers from it. Four things were checked.

**Sex contamination in earlier phenotype files.** Every phenotype TSV used in markdowns 06 and 07 was checked for its raw sex composition, independent of how the pipeline had filtered it. Three of the six files (Morgante starvation, chill coma, cuticle HC) do contain both sexes at the file level, as already documented in `phenotype-data/README.md`. All five prior runs used the default `--sex F`, and the script filters the phenotype file down to that sex before building the analysis matrix, so the reported results used only the female rows. No contamination in the actual numbers.

**DGRP line ID normalisation.** Four scripts each convert DGRP line identifiers to the same format independently: `check_morgante_overlap.py`, `run_fecundity_enet.py`, `run_dgrpool_phenotype.py`, and `prepare_unckless_data.py`. Tested against a set of sample IDs (`DGRP_021`, `DGRP229`, `DGRP_1`, and so on), all four produce identical output. Two of the four use byte-identical code; the third strips the underscore and prefix in a different order but ends up equivalent; the fourth works from an already-integer column rather than a string, since that is the format the Unckless source file uses, so it does not need the same parsing at all.

**StandardScaler discipline.** Every script that runs LOO-CV or GroupKFold (`run_pls_analysis.py`, `run_regularised_regression.py`, `run_perfly_pipeline.py`, `run_random_forest.py`, `run_dgrpool_phenotype.py`, `run_fecundity_enet.py`) instantiates a fresh scaler inside the fold loop, fits it on the training indices only, and applies it to the test indices without refitting. No leakage across any of the six.

**Cross-sex warning.** Confirmed working as described above: the header block prints a warning whenever `--sex` and `--spectral-sex` differ, and it fired correctly on every Unckless run.

## Data preparation

`scripts/prepare_unckless_data.py` reads the paper's Table S2 supplementary file directly (`phenotype-data/raw/016477_tables2.xlsx`), rather than going through DGRPool, since this dataset does not have a DGRPool entry. For each of the six measures it extracts three columns: pooled diet, high-glucose diet, and low-glucose diet, giving 18 output TSVs. Missing values, marked with a literal "." in the source file, are converted to NaN and dropped; one duplicate DGRP line (287) with identical values in both rows is dropped as well. The set of missing lines differs by diet condition rather than by measure, since a failed assay for one diet knocks out that line for all six measures under that diet: 145 of 153 source lines are valid for pooled, 147 for high-glucose, 150 for low-glucose.

## Results: 18 elastic net LOO-CV runs

Same model and hyperparameters as every other DGRPool phenotype run in this project: `ElasticNetCV(cv=3, l1_ratio=[0.5,0.7,0.9,0.95,1.0], alphas=30, max_iter=5000, tol=0.01)`, leave-one-line-out cross-validation, `StandardScaler` fitted inside each fold. The "Unckless r" column is their own published phenotypic correlation between that measure and starvation resistance (Table 1 of Unckless et al. 2015, correlated against starvation resistance as measured by Mackay et al. 2012), included for comparison since it is an independent check on which of these measures is actually expected to relate to starvation biology at all.

| Measure | Diet | n lines | Our CV R² | Unckless r with starvation resistance | Result |
|---|---|---|---|---|---|
| Glucose | Pooled | 77 | -0.026 | 0.246 (P<0.01) | Null, predictions collapse to the mean |
| Glucose | High-glucose | 77 | -0.026 | 0.246 (P<0.01) | Null, predictions collapse to the mean |
| Glucose | Low-glucose | 80 | -0.048 | 0.246 (P<0.01) | Null, predictions collapse to the mean |
| Glycerol | Pooled | 77 | +0.012 | 0.079 (ns) | Null, no collapse but no predictive power |
| Glycerol | High-glucose | 77 | +0.012 | 0.079 (ns) | Null, no collapse but no predictive power |
| Glycerol | Low-glucose | 80 | +0.007 | 0.079 (ns) | Null, predictions collapse to the mean |
| Glycogen | Pooled | 77 | -0.070 | 0.307 (P<0.001) | Null, no collapse but no predictive power |
| Glycogen | High-glucose | 77 | -0.045 | 0.307 (P<0.001) | Null, no collapse but no predictive power |
| Glycogen | Low-glucose | 80 | -0.025 | 0.307 (P<0.001) | Null, predictions collapse to the mean |
| Triglyceride | Pooled | 77 | -0.051 | -0.071 (ns) | Null, predictions collapse to the mean |
| Triglyceride | High-glucose | 77 | -0.006 | -0.071 (ns) | Null, predictions collapse to the mean |
| Triglyceride | Low-glucose | 80 | -0.030 | -0.071 (ns) | Null, predictions collapse to the mean |
| Protein | Pooled | 77 | -0.074 | -0.113 (ns) | Null, predictions collapse to the mean |
| Protein | High-glucose | 77 | -0.026 | -0.113 (ns) | Null, predictions collapse to the mean |
| **Protein** | **Low-glucose** | **80** | **+0.066** | -0.113 (ns) | See below |
| MeanWeight | Pooled | 77 | -0.026 | 0.241 (P<0.01), wet weight | Null, predictions collapse to the mean |
| MeanWeight | High-glucose | 77 | -0.026 | 0.241 (P<0.01), wet weight | Null, predictions collapse to the mean |
| MeanWeight | Low-glucose | 80 | -0.032 | 0.241 (P<0.01), wet weight | Null, predictions collapse to the mean |

"Collapse" is the same diagnostic used for the fecundity null in markdown 06: if the standard deviation of the LOO-CV predictions falls below 20% of the standard deviation of the true values, the model has effectively given up and is predicting close to the training mean rather than finding any structure. `run_dgrpool_phenotype.py` now prints this ratio on every run, not only when it drops below the threshold, so the distinction between a collapsed null and a non-collapsed null with genuinely weak predictive power can always be checked directly.

## The protein, low-glucose result

Protein under the low-glucose diet is the only one of these 18 tests that did not fit the null pattern. Because it stood out, it got the same scrutiny as the fecundity null in markdown 06 rather than being taken at face value in either direction.

The prediction spread rules out the obvious false-positive explanation. Predictions had a standard deviation of 0.0070 against a true value standard deviation of 0.0180, a ratio of 0.388, well above the 0.2 threshold that flags a collapsed model elsewhere in this project. So this is not a case of the model predicting a near-constant value and getting a small positive R² by chance the way a near-collapsed model sometimes can. The model produced real, varying predictions.

That said, R²=0.066 is still weak on its own terms, nowhere near the 0.673 seen for starvation resistance in our own data, and this project has now run roughly two dozen independent phenotype and diet combinations through this pipeline with no correction for multiple comparisons. A single weak positive result out of that many tests is exactly what would be expected by chance alone, collapsed or not. The honest reading is that this is not a mean-collapse artefact, but it is also not a validated finding: a candidate worth a follow-up test if an independent protein dataset becomes available, not something to report as a result in its own right.

## Interpretation

Seventeen of eighteen tests here show no spectral signal, and the one exception does not survive scrutiny as anything more than a candidate for further checking. What makes this result more informative than a simple failure to replicate is glycogen: it is the strongest published correlate of starvation resistance in the Unckless data itself (r=0.307, P<0.001, stronger than glucose or wet weight), and the FTIR spectra show no ability to predict it under any of the three diet conditions. If the FTIR chemotype were picking up a general lipid or energy-storage signature that happens to correlate with starvation resistance, glycogen is exactly the measure it should have picked up. It did not.

Combined with the five phenotype targets from markdown 07, where only starvation resistance measured in our own lab showed a real signal (and starvation resistance measured by Morgante 2015 showed only a weak, cross-lab-noise-consistent one), the pattern across the whole project now points somewhere more specific than "FTIR detects lipid content." The signal looks tied to this lab's own starvation assay rather than to metabolic or energy-storage status as a general biological category, even when tested against a metabolic measure with a documented, independent link to starvation biology.

## Running total across the project

| Phenotype | Source | Result |
|---|---|---|
| Starvation resistance | Own lab EMMeans | Strong positive, R²=0.673 |
| Starvation resistance | Morgante 2015 | Weak positive, R²=0.041, consistent with cross-lab measurement noise |
| Fecundity, lifespan, chill coma recovery, cuticle HC n-C25 | DGRPool (4 tests) | Null |
| Six Unckless metabolic measures, pooled and high/low-glucose diets | Unckless 2015 (17 of 18 tests) | Null |
| Protein, low-glucose diet | Unckless 2015 (1 of 18 tests) | Weak, not a collapse artefact, unconfirmed |

24 distinct phenotype and diet combinations tested in total: 2 involving starvation resistance itself, 21 genuine nulls, and 1 weak candidate that has not been corrected for multiple comparisons and should not be treated as a finding.

## Output files

- `phenotype-data/raw/016477_tables2.xlsx`: raw supplementary Table S2 from the paper, not generated by any script here
- `phenotype-data/Unckless_<Measure>_<diet>.tsv`: 18 files, generated by `scripts/prepare_unckless_data.py`
- `results/DGRP/dgrpool_phenotype_summary.csv`: appended with 18 new rows, one per run
- `phenotype-data/README.md`: full results table and methodology notes
- `REPRODUCE.md` Section 7: exact commands for regenerating every file and result in this markdown
