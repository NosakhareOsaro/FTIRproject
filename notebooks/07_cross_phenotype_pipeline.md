# 07: Cross-Phenotype Pipeline: General-Purpose DGRPool Runner

**Script:** `scripts/run_dgrpool_phenotype.py`
**Status:** Complete
**Key result:** Across five phenotype targets, only starvation resistance shows a spectral signal (CV R² = 0.673 on our own EMMeans; +0.041 on an independent cross-lab replicate). Lifespan, chill coma recovery, and cuticle HC n-C25 are all null (CV R² between −0.05 and −0.06). The FTIR chemotype's predictive signal is highly specific, so the lipid-metabolism hypothesis from markdown 06 needs refining rather than abandoning.

---

## What the pipeline does

Markdown 06 tested one additional phenotype (lifetime fecundity) using a bespoke script (`run_fecundity_enet.py`) that duplicated most of the starvation-resistance pipeline with the phenotype hardcoded. To test more phenotypes without copy-pasting the same ~150 lines each time, `scripts/run_dgrpool_phenotype.py` generalises this into a single command-line tool.

Given any DGRPool-format phenotype TSV (columns `DGRP`, `sex`, `value`), the script:

1. Filters to the requested sex (default `F`, matching all analyses so far)
2. Normalises DGRP line IDs, stripping underscores and leading zeros (`DGRP_021` → `DGRP21`) so DGRPool's zero-padded IDs match our spectral line IDs
3. Loads `DGRPFTIR.dat` via `load_ftir()`, filters to the same sex, and averages spectra within each DGRP line (the "Order B" collapse established in markdown 02)
4. Inner-joins the spectral lines against the phenotype lines and reports the overlap, including which spectral lines have no phenotype value
5. Runs elastic net with leave-one-line-out cross-validation, identical hyperparameters to `run_regularised_regression.py`: `ElasticNetCV(cv=3, l1_ratio=[0.5,0.7,0.9,0.95,1.0], alphas=30, max_iter=5000, tol=0.01)`, `StandardScaler` fitted inside each training fold
6. Reports n lines, CV R², RMSE, and Spearman ρ, with an explicit warning if predictions collapse toward the training mean (prediction SD < 20% of true SD), since this is the signature of a null result and Spearman ρ becomes numerically unreliable in that regime (documented first in markdown 06)
7. Appends one row to `results/DGRP/dgrpool_phenotype_summary.csv`, so results accumulate automatically across runs rather than needing to be manually copied into a table

### Validation

Before trusting the script on any new phenotype, it needed to reproduce a known result. `Emmeans.csv` was reformatted into a mock DGRPool TSV (`phenotype-data/S00_EMMeans_starvation.tsv`, same `DGRP`/`sex`/`value` columns, `sex` set to `F` for all 108 rows) and run through the pipeline as a smoke test. It reproduced **CV R² = 0.673**, matching `run_regularised_regression.py` almost exactly (small differences are due to a slightly wider `l1_ratio` grid). This confirms the general-purpose script is a faithful reimplementation, not just superficially similar code.

---

## Complete cross-phenotype results

| Phenotype                                   | Study              | n lines | CV R²      | RMSE   | Spearman ρ        | Signal? |
| ------------------------------------------- | ------------------ | ------- | ---------- | ------ | ----------------- | ------- |
| Starvation resistance (EMMeans, smoke test) | Internal PSM model | 108     | **+0.673** | 0.4213 | +0.817            | Yes     |
| Starvation resistance                       | Morgante 2015      | 104     | +0.041     | 12.952 | +0.150            | Weak    |
| Lifespan                                    | Ivanov 2015        | 104     | −0.052     | 9.917  | −0.012            | No      |
| Chill coma recovery                         | Morgante 2015      | 95      | −0.060     | 5.332  | −0.305 (artefact) | No      |
| Cuticle HC n-C25                            | Dembeck 2015       | 91      | −0.052     | 0.0258 | −0.947 (artefact) | No      |

The two Spearman ρ values marked "artefact" are the same LOO mean-shift artefact documented in markdown 06: when predictions collapse to a near-constant value, holding out a high-true-value line slightly shifts the training mean, producing a spurious monotone trend across folds that inflates |ρ| without reflecting genuine rank agreement. R² is the reliable metric in these cases, not ρ.

---

## Interpreting each result honestly

A null CV R² is not automatically informative: it could mean "no biological signal," "wrong evaluation," "data quality problem," or "measurement mismatch." Each result was checked individually rather than assumed to mean the same thing.

### Starvation resistance (Morgante 2015): weak signal, not absence of signal

R² = +0.041 is much lower than our own EMMeans (R² = 0.673), but this is not evidence that FTIR fails to predict starvation resistance in general. Morgante's measurement is an **independent replicate of the same phenotype**, collected in a different lab, at a different time, using a different assay protocol. The correlation between our EMMeans and Morgante's means, on the 104 overlapping lines, is only Pearson r = 0.428 (p = 5.88×10⁻⁶; `scripts/check_morgante_overlap.py`). Two independent measurements of the _same_ underlying trait only agree with each other at r ≈ 0.43, and this is the ceiling on how well any model trained on one could ever predict the other, regardless of the predictor's quality. A weak but positive and significant correlation between the two starvation measures is expected; the modest R² = 0.041 is consistent with real signal attenuated by cross-lab phenotype noise, not with the FTIR chemotype failing.

### Lifespan: genuine null

R² = −0.052, predictions collapse to the training mean (documented in the run output). Lifespan is a whole-organism integrative trait shaped by many physiological systems (metabolic rate, immune function, oxidative stress resistance, reproductive trade-offs), of which lipid storage is only one contributor among many. There is no strong a priori reason the FTIR lipid signature should dominate a trait this diffuse, and the result is consistent with that: no detectable signal.

### Chill coma recovery: genuine null, different stress axis

R² = −0.060, predictions collapse to the mean. Chill coma recovery measures how quickly a fly regains motor function after cold-induced paralysis. This is a stress-resistance phenotype, like starvation resistance, but the underlying physiology is different: cold tolerance in _Drosophila_ is linked to membrane lipid composition (fluidity, unsaturation) and ion homeostasis during chilling, rather than to bulk triglyceride energy reserves. The FTIR C-H stretching signal that predicts starvation resistance reflects total lipid quantity, not membrane lipid composition or saturation state, so there is no strong mechanistic reason to expect the same spectral region to carry a chill coma signal. The null result is informative: it shows the FTIR signal is specific to energy-reserve lipid content, not to "any lipid-related trait."

### Cuticle HC n-C25: genuine null, but likely a measurement mismatch rather than a true absence of signal

R² = −0.052, and this null was checked most carefully because it is the most surprising: cuticular hydrocarbons are laid down directly on the body surface, the same location FTIR measures. If any phenotype should be predictable from an ATR-FTIR surface measurement, cuticle chemistry should be it.

**Distribution check (91 lines):** mean 0.0716, SD 0.0252, strongly right-skewed (skewness 2.50), with one clear outlier (DGRP303 = 0.220, nearly double the next-highest line). Re-running the LOO-CV independently confirmed predictions genuinely collapse to a near-constant value (prediction SD = 0.0015 vs true SD = 0.0252, a ratio of 5.8%), so this is a real null, not a coding error or a skew/outlier artefact inflating the reported R².

**Why the null is not surprising on reflection: the measurement mismatch.** FTIR's C-H stretching region (~2,900-3,000 cm⁻¹, the same band driving the starvation prediction in markdown 03) is sensitive to _total alkyl C-H content_: it sums the absorbance contribution of every hydrocarbon chain on the cuticle surface, dominated by whichever hydrocarbon species are most abundant. n-C25 (pentacosane) is a single hydrocarbon species, and Dembeck's phenotype is its **relative abundance as a proportion of total CHC profile**, not its absolute quantity. A DGRP line could have a high n-C25 _proportion_ while having a low total CHC quantity, or vice versa: FTIR cannot distinguish "more of everything" from "more of this one component, less of the others," because a proportion is a compositional measure and FTIR gives a bulk absorbance measure. Predicting a single component's relative share of a multi-component mixture from a bulk spectral signal is a fundamentally harder, and possibly different-in-kind, problem than predicting bulk quantity, which is what the starvation resistance prediction actually relies on (bulk lipid content driving the C-H peak).

This means the cuticle HC null result should **not** be read as "FTIR cannot predict cuticle chemistry." It should be read as "FTIR did not predict this specific compositional measure of cuticle chemistry," a narrower and more defensible claim. A phenotype more aligned with what FTIR actually measures (total CHC abundance, or overall cuticle lipid load) is a better test of whether FTIR carries cuticle-surface lipid information.

---

## What the overall pattern means

Markdown 06 posed the hypothesis that "FTIR predictability tracks the degree to which a phenotype is determined by surface lipid chemistry." The results in this markdown complicate that hypothesis rather than confirming it: fecundity, lifespan, chill coma recovery, and cuticle HC n-C25 are all lipid-adjacent or surface-adjacent phenotypes in some loose sense, and all four show no signal.

The more precise conclusion is: **the FTIR signal detected so far is specific to starvation resistance (or its measurement, the EMMean), not to "lipid" as a broad category.** The mechanistic story from markdown 03, that the C-H stretching region reflects bulk triglyceride/energy-reserve lipid, is still the best available explanation for _why_ starvation resistance specifically is predictable, since starvation survival time is mechanistically almost synonymous with fat reserve size. But the corollary "any lipid-related trait should be predictable" does not hold, because:

- Lifespan and chill coma recovery are lipid-_adjacent_ at most, not lipid-_determined_: many other physiological systems dominate
- Cuticle HC n-C25 is a compositional (relative) measure, not a bulk quantity measure, and FTIR is fundamentally a bulk measurement technique

The refined hypothesis is that **FTIR predicts phenotypes that are directly and dominantly determined by bulk energy-reserve lipid quantity**, a narrower and more testable claim than "lipid metabolism in general." Testing this properly requires phenotypes chosen to match what FTIR actually measures, not just phenotypes that happen to involve lipids in some way.

---
