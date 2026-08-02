# 09: Lifespan and Fecundity Extension (Durham 2014, Huang 2020)

**Script:** `scripts/prepare_lifespan_fecundity_data.py`, `scripts/run_dgrpool_phenotype.py`.
**Status:** Complete
**Key result:** Five of six new elastic net LOO-CV tests show no spectral signal. Combined with the earlier Ivanov 2015 lifespan null (markdown 07), lifespan now shows no reliable FTIR signal across 5 independent tests spanning both virgin and mated females. The mating-status hypothesis that motivated this extension does not explain the pattern. The two Durham fecundity results (week 1, week 7) are also both null, ruling out early- and late-life fecundity as alternative explanations. One weak candidate survives scrutiny (Huang lifespan at 25°C, R²=0.056), and it is tied to the standard DGRP rearing temperature rather than to mating status.

---

## Why these two datasets

The lifespan null reported in markdown 07 (Ivanov et al. 2015, R²=−0.052) left an obvious confound unaddressed: Ivanov's females were unmated (virgin). If the FTIR chemotype tracks something about reproductive physiology, or if mating itself changes the lipid/metabolic state that the spectra are thought to reflect, a virgin-only lifespan test could miss a real signal that would show up in mated females. This was the exact check that was requested: test whether FTIR predicts lifespan using studies where females were mated, and separately, test age-specific fecundity using the youngest and oldest measured timepoints from Durham et al. 2014, motivated by the well-established relationship between lipid content and ageing in _Drosophila_. If the FTIR signal reflects lipid reserves, it might be expected to track fecundity most strongly early in life (when reserves are being allocated to egg production) or late in life (when reserves and fecundity both decline).

Two DGRPool sources cover this. Durham, Magwire, Stone & Leips (2014, study 12) measured mated female lifespan and fecundity at multiple ages; week 1 and week 7 are the youngest and oldest timepoints in their design. Huang et al. (2020, study 40) measured mated female lifespan at three rearing temperatures (18°C, 25°C, 28°C), giving three further mated-lifespan tests independent of the Durham cohort.

## The Huang raw-data structure

The three Durham files (`S12_lsm_*_original.tsv`) arrived already in DGRPool line/sex/value format, one row per line, matching every other phenotype file used in this project. The three Huang files (`S40_Lifespan_*_original.tsv`) did not: each is individual fly-level raw data, 23,000+ rows, both sexes present, matching the paper's reported design of roughly 72 flies per sex/line/temperature. Confirmed directly before writing any processing code: `S40_Lifespan_18C_original.tsv` has 11,683 female and 11,345 male rows across 183 unique DGRP lines; `_25C` has 12,155 female and 11,905 male rows across 186 lines; `_28C` has 11,657 female and 11,494 male rows across 177 lines. No missing values in any of the three. These needed filtering to females and averaging per line before they were usable as a DGRPool-format phenotype target, unlike every file used in markdowns 06–08, which arrived pre-aggregated.

## Data preparation

`scripts/prepare_lifespan_fecundity_data.py` handles the two sources differently, matching the structural difference above. The three Durham files are copied through unchanged into `Durham_Lifespan_mated.tsv`, `Durham_Week1_Fecundity.tsv`, and `Durham_Week7_Fecundity.tsv` (189, 189, and 135 rows respectively). No filtering or aggregation needed. The three Huang files are filtered to `sex == 'F'`, grouped by DGRP line, and averaged into `Huang_Lifespan_18C_female.tsv`, `Huang_Lifespan_25C_female.tsv`, and `Huang_Lifespan_28C_female.tsv` (183, 186, and 177 female line-means respectively).

One value needs an explicit note: DGRP_042 in the Durham Week 7 file is negative (−0.165). This is expected and correct, not a data error. Durham's reported values are least-squares means corrected for body size and block, not raw egg counts, so a line with fecundity well below the panel's LS-mean baseline can legitimately land below zero after correction. The Week 1 file has one mildly negative value for the same reason (−0.145), for the same reason.

## Results: 6 elastic net LOO-CV runs

Same model and hyperparameters as every other DGRPool phenotype run in this project: `ElasticNetCV(cv=3, l1_ratio=[0.5,0.7,0.9,0.95,1.0], alphas=30, max_iter=5000, tol=0.01)`, leave-one-line-out cross-validation, `StandardScaler` fitted inside each fold. All six use default `--sex F`/`--spectral-sex F`: these are same-sex comparisons, unlike the cross-sex Unckless work in markdown 08, since both the FTIR spectra and every phenotype tested here are female.

| Phenotype                 | Study          | n lines | CV R²      | Spearman ρ | SD ratio  | Result                                    |
| ------------------------- | -------------- | ------- | ---------- | ---------- | --------- | ----------------------------------------- |
| Lifespan (mated)          | Durham 2014    | 96      | −0.085     | −0.970     | 0.103     | Null, predictions collapse to the mean    |
| Fecundity Week 1          | Durham 2014    | 96      | −0.058     | −0.686     | 0.049     | Null, predictions collapse to the mean    |
| Fecundity Week 7          | Durham 2014    | 76      | −0.006     | +0.219     | 0.334     | Null, no collapse but no predictive power |
| Lifespan 18°C (mated)     | Huang 2020     | 99      | −0.023     | −0.999     | 0.014     | Null, predictions collapse to the mean    |
| **Lifespan 25°C (mated)** | **Huang 2020** | **102** | **+0.056** | **+0.263** | **0.294** | See below                                 |
| Lifespan 28°C (mated)     | Huang 2020     | 99      | −0.034     | −1.000     | 0.035     | Null, predictions collapse to the mean    |

"SD ratio" is the same collapse diagnostic used since markdown 06: the standard deviation of the LOO-CV predictions divided by the standard deviation of the true values. Below 0.2, the model has effectively given up and is predicting close to the training mean; `run_dgrpool_phenotype.py` prints this ratio on every run. Four of the six new tests show extreme Spearman ρ (|ρ| > 0.9) despite negative or near-zero R². This is the LOO mean-shift artefact documented since markdown 06, not a genuine rank correlation, and should not be read alongside the R² as if it were independent evidence.

## The Huang 25°C result

Lifespan at 25°C is the only one of these six tests that did not fit the null pattern, so it gets the same scrutiny as the fecundity null in markdown 06 and the protein/low-glucose result in markdown 08, rather than being taken at face value in either direction.

The prediction spread rules out the obvious false-positive explanation: predictions had a standard deviation of 3.065 against a true value standard deviation of 10.411, a ratio of 0.294, above the 0.2 threshold that flags a collapsed model. This is not a case of a near-constant model producing a small positive R² by chance the way a collapsed model sometimes can. The model produced real, varying predictions.

That said, R²=0.056 is weak on its own terms, and this project has now run 30 independent phenotype/diet/temperature tests through this pipeline (arithmetic below) with no correction for multiple comparisons. A single weak positive result out of that many tests is close to what would be expected by chance alone, non-collapsed or not. It is worth noting what is specific about this one condition: 25°C is the standard DGRP rearing temperature, the condition under which the DGRP lines are normally maintained and under which most other DGRP phenotyping (including our own FTIR spectra) is collected. That is a more defensible observation than any claim tied to mating status, since 18°C and 28°C are both also mated-female Huang cohorts and both are null. Mating status cannot be what distinguishes the 25°C result from the other two Huang temperatures. The honest reading is that this is not a mean-collapse artefact, but it is also not a validated finding: a candidate worth a follow-up test if an independent 25°C lifespan dataset becomes available, not something to report as a result in its own right.

## Interpretation

**Lifespan across mating status.** Including the Ivanov 2015 result from markdown 07 (virgin females, R²=−0.052), five independent lifespan tests have now been run against the FTIR spectra: one virgin cohort (Ivanov) and four mated cohorts (Durham, and Huang at 18°C, 25°C, and 28°C). Four of the five are null. This should be framed precisely: it does not show that mating status is irrelevant to FTIR spectra or to lifespan in general, only that **the specific hypothesis motivating this extension (that Ivanov's null was an artefact of using virgin females) does not explain the pattern of nulls**. Three of the four mated cohorts are null under the same diagnostic used throughout this project, so mated status alone does not unlock a signal that virgin status was suppressing. The one exception (Huang 25°C) is better explained by rearing temperature, a factor that varies across the mated cohorts themselves, than by mating status, which does not.

**Age-specific fecundity.** The two Durham fecundity results, week 1 (youngest measured timepoint) and week 7 (oldest), are both null. This was tested specifically because a lipid-content explanation for the FTIR signal might predict a fecundity relationship strongest at one age extreme or the other, as reserves are allocated to early reproduction or depleted with age. Neither timepoint shows a spectral relationship, ruling out early- or late-life fecundity as an alternative phenotype where the hypothesised lipid signal might have been easier to detect than in the lifetime-fecundity null already reported in markdown 06.

Combined with every prior null in markdowns 07 and 08, the pattern remains what it has been since the Unckless glycogen result: the FTIR signal detected for starvation resistance in this lab does not generalise to lifespan (virgin or mated), fecundity (lifetime, week-1, or week-7), chill coma recovery, cuticle hydrocarbons, or six independent metabolic measures under three diet conditions. Nothing in this extension changes that picture; it closes off two specific alternative explanations (virgin-only lifespan measurement, and fecundity measured only at one life stage) that could otherwise have been raised as reasons the earlier nulls were incomplete.

## Running total across the project

Verified directly against `results/DGRP/dgrpool_phenotype_summary.csv` (30 logged rows) plus the one test that used a dedicated script instead of the general-purpose runner (lifetime fecundity, markdown 06, `run_fecundity_enet.py`, not logged to that CSV) and one smoke-test duplicate in the CSV (the internal EMMeans validation run twice on different dates, same R²=0.673 both times, counted once here since it is the same test, not two).

**Category A, genuine signal (2):**

1. Starvation resistance, own lab (elastic net on raw spectra): R²=+0.673
2. Starvation resistance, Morgante 2015 (cross-lab): R²=+0.041

**Category B, null (26):**

| Count | Source                                                                                       |
| ----- | -------------------------------------------------------------------------------------------- |
| 1     | Lifetime fecundity, DGRPool study 18 (markdown 06)                                           |
| 1     | Lifespan, Ivanov 2015, virgin (markdown 07)                                                  |
| 1     | Chill coma recovery, Morgante 2015 (markdown 07)                                             |
| 1     | Cuticle HC n-C25, Dembeck 2015 (markdown 07)                                                 |
| 17    | Unckless 2015 metabolic measures, 17 of 18 measure×diet combinations (markdown 08)           |
| 5     | This markdown: Durham lifespan (mated), Durham week 1, Durham week 7, Huang 18°C, Huang 28°C |

Running sum: 1 + 1 + 1 + 1 + 17 + 5 = **26**

**Category C, flagged candidates (2):**

1. Unckless protein, low-glucose diet: R²=+0.066 (markdown 08)
2. Huang lifespan, 25°C (mated): R²=+0.056 (this markdown)

**Grand total: 2 + 26 + 2 = 30 independent phenotype/diet/temperature tests.**

Cross-check against the CSV directly: 30 logged rows = 2 rows for the smoke-test duplicate (both Category A test #1, the same result run twice) + 1 row for Category A test #2 (Morgante) + 25 of the 26 Category B nulls (all except lifetime fecundity, which was never logged to this CSV because it ran through a separate script) + 2 Category C rows = 2 + 1 + 25 + 2 = 30. ✓

30 tests run in total: 2 real signals, 26 genuine nulls, and 2 weak candidates that have not been corrected for multiple comparisons and should not be treated as findings.

## Output files

- `phenotype-data/raw/S12_lsm_Lifespan_original.tsv`, `S12_lsm_Week1_Fecundity_original.tsv`, `S12_lsm_Week7_Fecundity_original.tsv`: raw Durham et al. 2014 DGRPool downloads, not generated by any script here
- `phenotype-data/raw/S40_Lifespan_18C_original.tsv`, `S40_Lifespan_25C_original.tsv`, `S40_Lifespan_28C_original.tsv`: raw Huang et al. 2020 DGRPool downloads, individual fly-level, not generated by any script here
- `phenotype-data/Durham_Lifespan_mated.tsv`, `Durham_Week1_Fecundity.tsv`, `Durham_Week7_Fecundity.tsv`: generated by `scripts/prepare_lifespan_fecundity_data.py`
- `phenotype-data/Huang_Lifespan_18C_female.tsv`, `Huang_Lifespan_25C_female.tsv`, `Huang_Lifespan_28C_female.tsv`: generated by `scripts/prepare_lifespan_fecundity_data.py`
- `results/DGRP/dgrpool_phenotype_summary.csv`: appended with 6 new rows, one per run
- `phenotype-data/README.md`: full results table and methodology notes
- `PROJECT_NOTES.md` §7e: results table and interpretation
- `REPRODUCE.md` Section 8: exact commands for regenerating every file and result in this markdown
