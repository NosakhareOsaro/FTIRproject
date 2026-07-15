# 06: Cross-Phenotype Analysis: Lifetime Fecundity

**Script:** `scripts/run_fecundity_enet.py`
**Status:** Complete
**Key result:** Elastic net CV R² = −0.109 on lifetime fecundity, with no spectral signal detected. The FTIR chemotype does not predict fecundity, suggesting the spectral signal is specific to starvation resistance rather than a generic predictor of life-history variation.

---

## What this step does and why

The starvation resistance analysis (markdowns 02-05) established that FTIR spectra of DGRP flies contain a reliable signal for predicting starvation resistance, driven primarily by the C-H stretching region (lipid content). A natural question follows: is this a general property of the spectra, or is it specific to starvation resistance?

If the FTIR chemotype is a broad readout of overall physiological state, it might predict many life-history traits. If it is specific to lipid-related traits, it would predict starvation resistance (which depends directly on fat reserves) but not other traits with different biochemical bases.

This step tests the pipeline on a second phenotype, **lifetime fecundity**, to begin answering this question. The same elastic net model, the same cross-validation design, and the same spectra were used. Only the phenotype target changed.

---

## The phenotype: lifetime fecundity

**Source:** DGRPool Study 18, Durham et al. 2014
**File:** `phenotype-data/S18_LifeFecundity_mean.tsv`
**Measure:** Mean number of eggs laid per female over a lifetime (not corrected for survival)
**Coverage:** 189 DGRP lines, all female, no missing values

Fecundity was chosen as the second phenotype for several reasons. It is a core life-history trait measured on the same DGRP lines. It has good female coverage (189 lines) and strong overlap with the spectral dataset. It is biologically distinct from starvation resistance: fecundity depends on reproductive investment and resource allocation to egg production, which is a different biochemical process than storing lipid reserves for starvation survival.

---

## Line overlap

After normalising DGRP line IDs to remove underscores (DGRPool uses `DGRP_100`; our spectral data uses `DGRP100`), the overlap was computed:

- Our spectral lines: 108
- DGRPool fecundity lines: 189
- Overlap: **96 lines** (12 of our 108 lines absent from the fecundity dataset)
- All 96 overlapping lines have valid (non-NA) fecundity values

The 12 missing lines (DGRP100, DGRP31, DGRP32, DGRP320, DGRP321, DGRP348, DGRP354, DGRP390, DGRP405, DGRP406, DGRP513, DGRP859) were simply not measured in Durham et al. 2014. They were dropped, with no imputation performed.

---

## Model and evaluation design

Identical to the starvation resistance analysis in markdown 04:

- **Method:** Elastic net (`ElasticNetCV`)
- **CV:** Leave-one-line-out (LOO) over the 96 overlapping lines
- **Hyperparameters:** cv=3, l1_ratio ∈ [0.5, 0.7, 0.9, 0.95, 1.0], 30 α candidates, max_iter=5000, tol=0.01
- **Scaling:** StandardScaler fitted inside each fold on the 95 training lines
- **Target:** Female lifetime fecundity mean (eggs per female)

The only changes from the starvation analysis are the phenotype target and the reduced sample size (96 lines instead of 108).

---

## Result

| Phenotype | n lines | Elastic net CV R² | RMSE |
|---|---|---|---|
| Starvation resistance | 108 | +0.673 | 0.4244 |
| **Lifetime fecundity** | **96** | **−0.109** | **20.913** |

CV R² = −0.109 means the elastic net model predicts fecundity **worse than simply predicting the mean fecundity value for every line**. A model that always outputs the mean fecundity would achieve R² = 0 by definition. Negative R² indicates the predictions are actively worse than this naive baseline.

---

## Diagnosing the result: why R² is negative

When a model with R² = 0 or below, it is worth checking that the result reflects genuine absence of signal rather than a bug. Three checks were run:

**Check 1: Prediction spread:**
The standard deviation of LOO-CV predictions was 2.8 eggs/female. The standard deviation of true fecundity values was 19.9 eggs/female. The model is essentially predicting the training mean for every test line, with only tiny variation between predictions. This is the expected behaviour of a regularised model when there is no signal: the penalty drives all coefficients to near zero and the model degenerates to predicting the mean.

**Check 2: Spearman ρ artefact:**
The raw LOO-CV output showed Spearman ρ ≈ −1, which looks alarming but is a known numerical artefact. When predictions are near-constant, a subtle bias emerges: when a high-fecundity line is held out, its inclusion in the training set would have slightly raised the training mean. Removing it slightly lowers the training mean, which slightly lowers the prediction for that line. This creates an artificial monotone negative trend across the 96 folds: high-fecundity lines get slightly lower predictions than low-fecundity lines, purely because of the mean-shift effect. With near-constant predictions, this tiny effect dominates the rank correlation and drives ρ toward −1. Spearman ρ is not reported for this result. R² = −0.109 is the correct metric.

**Check 3: Cross-phenotype correlation:**
The Pearson correlation between starvation EMMeans and fecundity means for the 96 overlapping lines is ρ = −0.04 (p = 0.68), effectively zero. Starvation resistance and lifetime fecundity are not correlated with each other in this dataset. This confirms the null result is not caused by the two phenotypes measuring the same thing in opposite directions.

All three checks are consistent: there is genuinely no detectable spectral signal for lifetime fecundity.

---

## Interpretation

The contrast between the two phenotypes is stark and interpretable:

| Phenotype | CV R² | Spectral signal? |
|---|---|---|
| Starvation resistance | +0.673 | Yes, strong |
| Lifetime fecundity | −0.109 | No |

**Why starvation resistance but not fecundity?**

The FTIR signal is primarily driven by the C-H stretching region (~2,900-3,000 cm⁻¹), which reflects lipid content on the cuticle surface (markdown 03). Starvation resistance in Drosophila is directly linked to lipid storage: flies with larger fat reserves survive longer without food. Lines that store more lipid have higher FTIR absorption in this region and higher starvation resistance. The spectral signal and the phenotype share a common biochemical basis.

Lifetime fecundity depends on different biochemistry. Reproductive investment in Drosophila involves resource allocation to ovary development, vitellogenin production, and egg maturation, processes driven by protein synthesis and carbohydrate metabolism more than by bulk lipid content. The FTIR measurement at the cuticle surface does not directly reflect the biochemical variation underlying fecundity.

This result is consistent with Adam's observation after seeing the starvation results: the FTIR prediction approach may be most useful for traits where lipid metabolism is involved. The fecundity null result supports this hypothesis and narrows the scope of what the chemotype actually measures.

**This is a meaningful negative result, not a failure.** A method that predicts everything would be uninformative about biology. The specificity of the FTIR signal to starvation-relevant lipid variation, and its absence for fecundity, is itself a finding about the biochemical content of the measurement.

---

## Next steps for cross-phenotype analysis

Fecundity is one data point. To build a broader picture of which phenotypes FTIR can and cannot predict, the same pipeline should be run on additional DGRPool traits. The most informative next phenotypes to test are those with different biochemical bases:

- **Lifespan** (Ivanov 2015, 197 female lines): is longevity predicted by the lipid chemotype?
- **Chill coma recovery** (Morgante 2015, 174 female lines): a stress resistance trait with different biochemistry to starvation
- **Cuticle hydrocarbons** (Dembeck 2015, ~169 female lines): a direct biochemical measure from the cuticle surface; the most mechanistically direct prediction the FTIR chemotype should make

The hypothesis emerging from the starvation and fecundity results is that FTIR predictability tracks the degree to which a phenotype is determined by surface lipid chemistry. Testing this hypothesis across multiple phenotypes is the next phase of the project.

---

## Output files

No figures were generated for this analysis: the result (predictions collapsing to the training mean) does not produce an informative scatter plot. The result is reported numerically.

- `phenotype-data/S18_LifeFecundity_mean.tsv`: raw DGRPool phenotype file
- `phenotype-data/README.md`: provenance record updated with this file