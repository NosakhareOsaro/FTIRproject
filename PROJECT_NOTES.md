# PROJECT_NOTES.md

Context file for working on this project locally.
Drop this in the repo root. It summarises the agreed plan, data status, and
evaluation design so any assistant starting fresh has the full picture.

**Author:** Nosakhare Odionfo Osaro (MSc Bioinformatics, University of Glasgow)
**Supervisor:** Dr Adam Dobson
**Last updated:** 18 July 2026

---

## 1. What this project is

MSc dissertation. Working title: _Machine Learning Signatures of Microbiome
Perturbation_ (the official title is broader than the actual work, see below).

**Actual focus:** thorough exploration of multivariate methods for
phenotype prediction from FTIR spectra in _Drosophila melanogaster_. The
core question is: can we generate continuous vectors from spectra that
reliably predict phenotype in individuals not seen during training? This
builds directly on a lab pre-print (Ibrahim et al., bioRxiv 2026,
doi:10.64898/2026.03.22.713522), nicknamed "chemotyping", which uses ATR-FTIR + ML
to classify biological variation and predict starvation resistance.

**Framing:** the project is now "explore the method space thoroughly"
rather than "test a specific hypothesis." The FTIR/DGRP work is the primary
case study. Methods span dimensionality reduction (PCA, tSNE), methods that
do both reduction and prediction (PLS-DA, sPLS-DA), and direct prediction
(LASSO, ridge regression, elastic net, random forest). Regularised/frequentist
methods (PLS-DA, elastic net, LASSO) serve as baselines within the broader
comparison: the frequentist angle is not abandoned, just embedded in a wider
sweep. Common evaluation yardstick: held-out prediction performance (R² for
continuous targets, accuracy for discrete). A set of external clinical cohorts
(NoMIC, EATRIS-Plus, ROP, IBD, IMPACC, CMAISE) are referenced as motivation
only, not worked through, unless time allows after the FTIR work.

---

## 2. The gap this project fills

The pre-print has two limitations that this project addresses:

1. **No frequentist baseline.** It compares five ML classifiers (LR, SVM, kNN,
   RF, XGBoost) only to each other. It does NOT include PLS-DA or PLS
   regression: the field-standard chemometric baseline. Given the strong
   colinearity in spectral data, PLS is the natural comparator and its absence
   is a real gap. (Supervisor explicitly endorsed adding it.)

2. **Continuous outcomes collapsed to categories.** Every continuous or ordinal
   outcome is discretised before ML: starvation resistance (continuous EMMeans
   → 20th/80th percentile binary), dietary restriction (5/10/20% yeast →
   3-class), age (9d/4wk/6wk → 3-class). The paper's own thesis, "chemical
   state encodes biological state," is fundamentally a continuous claim, so
   regression is a more principled test of it than classification.

---

## 3. The three contributions (agreed with supervisor, updated 24 Jun 2026)

1. **Systematic method comparison** across the full space of multivariate
   approaches: dimensionality reduction only (PCA, tSNE), methods that reduce
   and predict (PLS-DA, sPLS-DA), and direct prediction (LASSO, ridge,
   elastic net, random forest, SVR). Both the dimensionality-reduction and
   direct-prediction families must be tested; PLS-DA, elastic net, and LASSO
   serve as the regularised/frequentist baselines within the comparison.
2. **Build a regression pipeline that predicts continuously.** Starvation
   resistance (108 DGRP EMMeans) is the primary target. Use **SVR** as the
   like-for-like counterpart to the pre-print's best classifier (SVM).
   Two evaluation settings: (a) per-fly spectra with line-stratified CV,
   (b) collapsed mean spectrum per DGRP line.
3. **Extend to other DGRPool phenotypes** (lifespan, fecundity, age at
   measurement). The DGRP panel has extensive public phenotype data
   (https://dgrpool.epfl.ch/), so the same spectra can be tested against many
   continuous targets. Starvation resistance remains the positive control
   (signal is confirmed to be there).

---

## 4. CRITICAL: evaluation design (do not get this wrong)

Both the FTIR assay and the starvation assay are **destructive**: the same
fly cannot be measured for both. So FTIR spectra and phenotype can only be
linked at the **line (genotype) level**, not the individual fly level.

Consequences for the pipeline:

- **Training unit:** individual fly (~16 flies per line × 108 lines ≈ 1,700
  spectra). This preserves within-line spectral variation as signal and gives
  a much larger training set than 108 line-means.
- **Cross-validation:** MUST be **line-stratified** (group k-fold by DGRP line).
  No fly from a given line may appear in both train and test folds. Random
  fly-level CV would leak the line-level target across folds and inflate
  performance; this is the single most important correctness requirement.
- **Evaluation:** collapse per-fly predictions to the line level (average
  predictions within each test-fold line → one prediction per line), then
  compute metrics against that line's true EMMean.
- **Metrics:** R², RMSE, Spearman rho, Pearson r. Consider weighting by the
  EMMean standard errors (already available in the R output) since lines are
  estimated with varying precision.

This is a genuine methodological subtlety. It needs to be implemented
explicitly AND understood well enough to defend in the viva.

---

## 5. Data

Repo: github.com/r-ib-code/FTIRproject (Rita Ibrahim).
As of 15 June 2026 all spectral files contain real data (earlier versions had
empty placeholders; Rita re-uploaded). Sample counts verified against the
pre-print figures: they match.

| File                         | Spectra | Experiment                           |
| ---------------------------- | ------- | ------------------------------------ |
| `FTIR-data/DGRPFTIR.dat`     | 1,772   | DGRP starvation (primary)            |
| `FTIR-data/SexGenoFTIR.dat`  | 1,684   | Sex + genotype (944F / 740M)         |
| `FTIR-data/Diet2FTIR.dat`    | 1,187   | Dietary restriction (5/10/20% yeast) |
| `FTIR-data/AgeFTIR.dat`      | 2,046   | Ageing (9d / 4wk / 6wk)              |
| `FTIR-data/Diet-AgeFTIR.dat` | 1,870   | Age × diet                           |
| `FTIR-data/Diet1FTIR.dat`    | 177     | High-fat diet (control / high-lipid) |

**File format:** tab-separated. First 4-5 columns are metadata
(`Genot.`, `Sex`, `Age`, `Diet`, sometimes `StoTime`), then ~1,723 wavenumber
columns running from 3900 down to ~456 cm⁻¹. Column count varies (1727 vs 1728) depending on whether `StoTime` is present: the loader must handle both.

**Metadata coding is inconsistent and needs normalising** in the loader:

- Diet: `D05`/`D10`/`D20` (yeast %), `LI`/`SY` (lipid vs standard yeast)
- Age: `09D`, `04W`, `06W`, `03W`, `UNK`
- Sex: `F`, `M`, and `S` (`S` in Diet1 = confirmed female by data owner); some `UNK`/`SY` placeholders

Phenotype data:

- `Survival-data/DGRP-starvationresistance.csv`: individual fly survival
  events (2,887 rows) → fed into the R survival model.
- `Survival-data/Diet2-lifespan.csv`, `Fecundity-data/Diet2-fecundity.csv`:
  other phenotypes.

---

## 6. Existing code in the repo

- `scripts/DGRP_survival_analysis.R`: fits a parametric survival model
  (rms::psm, logistic dist) across the 108 DGRP lines, extracts EMMeans per
  line (the continuous starvation-resistance coefficient → `Emmeans.csv`),
  then bins at 20th/80th percentiles into sensitive/resistant. Also runs tSNE.
  **The EMMeans this produces are the regression target.**
- `scripts/FTIR-script.ipynb`: one big `MIRSPIPELINE()` function (by Rita
  Ibrahim & Mario Gonzalez Jimenez). Handles all experiments. XGBoost for
  feature (wavenumber) selection, then trains 5 classifiers with GridSearchCV,
  evaluates with stratified k=20 CV. All CLASSIFIERS, no regressors.
- Upstream QC uses the `bad-blood` package
  (github.com/magonji/bad-blood): discards low-intensity spectra, atmospheric
  interference, etc. CHECK whether the .dat files are already QC'd or whether
  this needs running first.

---

## 7. Work order (updated 25 Jun 2026)

Steps 1–4 complete. Steps 5–6 are the active frontier.

1. **[DONE] Reproduce Rita's classification baseline** on `DGRPFTIR.dat`.
   Recovered ~88.5% resistant / ~83.9% sensitive SVM reclassification
   (paper: 88% / 84%). Gate passed.
2. **[DONE] Write a clean data loader** (`scripts/ftir_loader.py`): reads
   any of the 6 .dat files, normalises messy metadata, validates wavenumber
   axis.
3. **[DONE] Reproduce the R survival analysis** to regenerate `Emmeans.csv`
   (+ SEs). EMMeans match Morgante 2015 (r=0.428, 104/108 lines overlap;
   corrected 2026-07-03 from an earlier 93-line count that undercounted
   zero-padded Morgante IDs like `DGRP_021`).
4. **[DONE] Build the method-comparison pipeline** on line-mean spectra
   (108 lines × 1,723 wavenumbers, LOO-CV). See §7a for full results.
5. **[IN PROGRESS] Extend to other DGRPool phenotypes**: assess whether
   spectral signal generalises beyond starvation resistance. Fecundity,
   lifespan, chill coma recovery, and cuticle HC n-C25 tested (§7b, §7c): no
   signal on any. Six Unckless et al. 2015 metabolic measures across three
   diet conditions tested next at Adam's request (§7d): 17 of 18 null, one
   weak unconfirmed candidate (protein, low-glucose diet).
6. **Meeting with Vinny Davies** (mathematician, potential collaborator):
   Monday 29 June 2026, 2 pm. No fixed agenda; may shape the mathematical
   approach, particularly around how to compare dimensionality-reduction
   methods with direct prediction methods, and whether Bayesian/GP approaches
   are worth adding given small N.

---

## 7a. Results so far (25 Jun 2026)

### PCA compression (`scripts/run_compression_analysis.py`)

- **Individual-fly PCA (Order A):** PC1 vs EMMean r = −0.05, p = 0.61:
  within-line noise swamps the genotype signal when PCA is fitted on all
  ~1,772 individual fly spectra.
- **Line-mean PCA (Order B):** averaging ~16 spectra per line first cancels
  within-line noise; PC1 of the 108-line-mean matrix vs EMMean r = +0.685,
  p = 2.92e−16. **Order B is the correct collapse order.**
- **Dimensionality:** only 4 PCs reach 95% explained variance: FTIR spectra
  of _Drosophila_ are very low-dimensional (PC1 = 53.6%, PC2 = 32.2%).

### Method comparison: LOO-CV on 108 DGRP line-mean spectra

All methods use StandardScaler fitted inside each training fold. α/hyperparameters
selected by inner CV within the training fold only (test line never seen).

| Method                           | CV R²     | RMSE       | Spearman ρ | Script                                                          |
| -------------------------------- | --------- | ---------- | ---------- | --------------------------------------------------------------- |
| PCA + Ridge (4 PCs)              | 0.553     | 0.4927     | +0.743     | `run_compression_analysis.py` + `run_regularised_regression.py` |
| PLS (10 components, optimal)     | 0.623     | 0.4524     | +0.801     | `run_pls_analysis.py`                                           |
| Random forest (max_features=0.3) | 0.540     | 0.4998     | +0.729     | `run_random_forest.py`                                          |
| Ridge (raw 1,723 wn)             | 0.635     | 0.4451     | +0.809     | `run_regularised_regression.py`                                 |
| LASSO                            | 0.669     | 0.4266     | +0.813     | `run_regularised_regression.py`                                 |
| **Elastic net (best)**           | **0.673** | **0.4244** | **+0.816** | `run_regularised_regression.py`                                 |

**Key finding:** Direct sparse regression (elastic net, LASSO) on raw spectra
outperforms all other methods. Random forest (R²=0.540) performs similarly to
PCA+Ridge: the non-linear ensemble offers no advantage over linear methods
here, consistent with spectral data being highly collinear and the signal being
largely captured by a few broad spectral regions rather than non-linear
interactions. The L1 penalty is a more efficient compression for this phenotype
prediction task than either PLS latent structure or RF tree splits.

**PLS model selection:** performance peaks at n=10 components (CV R²=0.623) and
declines at n=15 (0.550) and n=20 (0.515), classic overfitting with n=108
lines.

**PLS loading vector:** prediction driven primarily by the C–H stretching region
(~2,900–3,000 cm⁻¹), consistent with lipid content as the main spectral
correlate of starvation resistance.

**RF hyperparameter selection:** modal best params across 108 LOO folds:
`max_features=0.3, n_estimators=100`. Inner 3-fold GridSearchCV over
{n_estimators: [100,300], max_features: ["sqrt","log2",0.1,0.3]}.

**Still to do on line-mean setting:** SVR, tSNE (visualisation). Random forest done (see table above).

### Per-fly pipeline: GroupKFold(10) (`scripts/run_perfly_pipeline.py`)

Training unit: individual fly spectra (~1,772 females). Outer CV:
GroupKFold(n_splits=10) with DGRP line as the group (no fly from a given line
appears in both train and test). After all folds, predictions are averaged
within each test-fold line → one predicted value per line, then compared to
EMMeans. This is the dissertation-grade evaluation.

| Method                             | Line R² | RMSE   | Spearman ρ | Pearson r |
| ---------------------------------- | ------- | ------ | ---------- | --------- |
| PLS (n_comp by inner GroupKFold-5) | 0.534   | 0.5029 | +0.794     | +0.802    |
| Ridge                              | 0.515   | 0.5130 | +0.767     | +0.764    |
| LASSO                              | 0.517   | 0.5121 | +0.781     | +0.768    |
| Elastic net                        | 0.518   | 0.5113 | +0.783     | +0.771    |

**Key finding:** PLS wins the per-fly setting (ranking inverted vs line-mean
LOO-CV, where sparse models dominated). The L1 sparse models drop ~0.15 R²
relative to their line-mean performance; they appear to overfit to
individual-fly spectral variation that does not generalise across lines.
All four methods are tightly clustered (line R² 0.515–0.534), suggesting
within-line noise is the dominant bottleneck rather than model class.

Per-fly R² (before line-averaging) is much lower (0.31–0.36); this is
expected: within-line spectral noise is real and cannot be learned from
line-level targets. The line R² values above are the correct dissertation
metric. Results also saved to `results/DGRP/perfly_metrics.csv`.

---

## 7b. Cross-phenotype results (25 Jun 2026)

### Lifetime fecundity (`scripts/run_fecundity_enet.py`)

Phenotype: DGRPool study 18, female lifetime fecundity means (eggs/female).
96 of 108 spectral lines have a fecundity value; 12 dropped (no DGRPool entry).
Same elastic net LOO-CV as starvation resistance (identical hyperparameters).

| Phenotype             | n lines | Elastic net CV R² | RMSE   |
| --------------------- | ------- | ----------------- | ------ |
| Starvation resistance | 108     | 0.673             | 0.4244 |
| Lifetime fecundity    | 96      | −0.109            | 20.913 |

**Finding: no spectral signal for lifetime fecundity.** The model predicts
approximately the training mean for every test point (prediction SD ≈ 2.8 vs
true SD ≈ 19.9). Elastic net selects maximum regularisation, driving all
coefficients near zero. R² = −0.109 (worse than a mean baseline).

**Methodological note on Spearman ρ:** the raw LOO output shows ρ ≈ −1, which
is a numerical artefact: with near-constant predictions, the LOO mean-shift
effect dominates (holding out a high-fecundity line slightly lowers the training
mean, so the model predicts slightly lower for high-true lines → artificial
monotone negative trend). ρ is not reported for this result; R² is the correct
metric. The starvation/fecundity cross-phenotype correlation is only ρ = −0.04
(p = 0.68), confirming no biological confound.

**Interpretation:** the FTIR spectral signal is specific to starvation
resistance (and its biochemical correlates, primarily lipid content), and does
not generalise to lifetime fecundity. This is a meaningful negative result for
the dissertation: it rules out the FTIR chemotype being a generic indicator of
all life-history variation.

---

## 7c. General-purpose cross-phenotype pipeline (3–4 Jul 2026)

`scripts/run_dgrpool_phenotype.py` generalises the fecundity script into a
command-line tool: any DGRPool-format phenotype TSV (`DGRP`, `sex`, `value`)
can be run through the same elastic net LOO-CV pipeline used for starvation
resistance. Each run appends one row to
`results/DGRP/dgrpool_phenotype_summary.csv`, so results accumulate
automatically across phenotypes.

**Validation:** run first against `S00_EMMeans_starvation.tsv` (our own
EMMeans reformatted as a mock DGRPool TSV) to confirm the script reproduces
the known result (R²≈0.673) before trusting it on real external phenotypes.

| Phenotype                                   | Study              | n lines | Elastic net CV R² | RMSE   | Spearman ρ        |
| ------------------------------------------- | ------------------ | ------- | ----------------- | ------ | ----------------- |
| Starvation resistance (EMMeans, smoke test) | Internal PSM model | 108     | +0.673            | 0.4213 | +0.817            |
| Starvation resistance                       | Morgante 2015      | 104     | +0.041            | 12.952 | +0.150            |
| Lifespan                                    | Ivanov 2015        | 104     | −0.052            | 9.917  | −0.012            |
| Chill coma recovery                         | Morgante 2015      | 95      | −0.060            | 5.332  | −0.305 (artefact) |
| Cuticle HC n-C25                            | Dembeck 2015       | 91      | −0.052            | 0.0258 | −0.947 (artefact) |

**Smoke test passed:** the S00 row exactly reproduces the
`run_regularised_regression.py` result, confirming the general-purpose script
is correct.

**Morgante starvation resistance (cross-lab validation):** R²=+0.041 is far
weaker than our own EMMeans (R²=0.673), but this is expected: Morgante's
measurement is an independent, noisy replicate of the same phenotype in a
different lab, not the same target the model was fitted against. The
line-level correlation between our EMMeans and Morgante's means is only
Pearson r=0.428 (p=5.88e-06; see `scripts/check_morgante_overlap.py`), so a
model trained on one is not expected to predict the other well even before
FTIR enters the picture. (Overlap was corrected 2026-07-03 from 93→104 lines:
the original script under-normalised zero-padded Morgante IDs.)

**Lifespan, chill coma recovery, cuticle HC n-C25: no spectral signal.** All
three show negative CV R² with predictions collapsing to the training mean
(see per-run notes in `phenotype-data/README.md`). The cuticle HC null result
is the most notable: cuticular hydrocarbons are the most mechanistically
direct cuticle-surface measurement tested against the FTIR chemotype, yet no
signal is detected. Combined with the fecundity, lifespan, and chill coma
nulls, this sharpens the picture from markdown 06: the FTIR signal so far
appears specific to starvation resistance itself (or its EMMean
representation), not to lipid content or cuticle chemistry as a general
category. This narrows rather than confirms the original lipid-metabolism
hypothesis and is worth discussing directly with Adam.

---

## 7d. Unckless metabolic cross-phenotype (15–18 Jul 2026)

Task 4 from the 15 July meeting with Adam: markdowns 06/07 showed our FTIR
spectra predict starvation resistance measured in our own lab, but not the
same phenotype measured by another lab, and not any other trait tested so
far. Adam asked whether the FTIR chemotype has any predictive power at all
outside our own starvation assay, and sent over Unckless, Rottschaefer &
Lazzaro (2015, G3, doi:10.1534/g3.114.016477), which measured six metabolic
indices in the DGRP (glucose, glycerol, glycogen, triglyceride, protein, wet
weight), each under pooled, high-glucose, and low-glucose diet conditions.
See `notebooks/08_unckless_metabolic_crossphenotype.md` for the full
writeup; summarised here.

**Data problem found before any analysis:** the paper's Materials and
Methods states all six measures were assayed in pools of 10 adult males per
line. `DGRPFTIR.dat` is female-only. Every comparison here is necessarily
cross-sex, matched by DGRP line/genotype, not by sex.

**Bug found and fixed:** `run_dgrpool_phenotype.py` used a single `--sex`
flag to filter both the phenotype file and the spectral data. Passing
`--sex M` would have filtered the female-only spectra down to zero rows and
aborted the run. Fixed by splitting into `--sex` (phenotype) and
`--spectral-sex` (spectra, default `F`), with a console warning whenever the
two differ. Default behaviour for every phenotype tested before this is
unchanged (both were always `F`).

**Pipeline audit run after the fix**, before trusting any new numbers from
the script: (1) raw sex composition of all six previously-used phenotype
files checked directly, confirmed no contamination despite three of the six
containing both sexes at the file level; (2) DGRP line ID normalisation
compared across `check_morgante_overlap.py`, `run_fecundity_enet.py`,
`run_dgrpool_phenotype.py`, and `prepare_unckless_data.py`, confirmed
identical output on sample IDs; (3) StandardScaler fit-on-training-fold-only
discipline confirmed in all six scripts that run LOO-CV or GroupKFold; (4)
the new cross-sex warning confirmed firing correctly on every Unckless run.

`scripts/prepare_unckless_data.py` converts the paper's raw Table S2
(`phenotype-data/raw/016477_tables2.xlsx`, not a DGRPool download) into 18
TSVs: 6 measures × 3 diet conditions. 145 valid lines for pooled, 147 for
high-glucose, 150 for low-glucose (missing values and one duplicate line
dropped identically per diet condition, differs by diet not by measure).

| Measure      | Diet         | n lines | CV R²  | Unckless r vs starvation |
| ------------ | ------------ | ------- | ------ | ------------------------ |
| Glucose      | all 3        | 77/77/80 | −0.026/−0.026/−0.048 | 0.246 (P<0.01) |
| Glycerol     | all 3        | 77/77/80 | +0.012/+0.012/+0.007 | 0.079 (ns) |
| Glycogen     | all 3        | 77/77/80 | −0.070/−0.045/−0.025 | 0.307 (P<0.001), strongest |
| Triglyceride | all 3        | 77/77/80 | −0.051/−0.006/−0.030 | −0.071 (ns) |
| Protein      | pooled/high  | 77/77   | −0.074/−0.026        | −0.113 (ns) |
| **Protein**  | **low**      | **80**  | **+0.066**            | −0.113 (ns) |
| MeanWeight   | all 3        | 77/77/80 | −0.026/−0.026/−0.032 | 0.241 (P<0.01), wet weight |

17 of 18 tests are genuine nulls. `run_dgrpool_phenotype.py` now prints the
prediction-SD/true-SD ratio on every run (not only when it drops below the
0.2 collapse threshold, the original behaviour), which lets a collapsed null
be distinguished from a null where predictions vary but simply aren't
accurate. Both patterns show up here (see notebook 08 for the full
per-measure breakdown).

**Protein, low-glucose diet (R²=+0.066) checked separately.** Prediction
SD/true SD ratio = 0.388, well above the 0.2 collapse threshold, so this is
not a mean-collapse artefact the way most of the nulls are. It is still
weak in absolute terms and uncorrected for multiple comparisons across the
roughly two dozen phenotype/diet tests run in this project so far. Treated
as a candidate for follow-up, not a finding.

**Interpretation:** glycogen is the strongest published correlate of
starvation resistance in the Unckless data itself (r=0.307, P<0.001,
stronger than glucose or wet weight), and the FTIR spectra show no ability
to predict it under any diet condition. If the FTIR chemotype tracked lipid
or energy-storage status generally, glycogen is exactly what it should have
picked up. Combined with the five phenotype targets in §7c, the signal now
looks specific to this lab's own starvation assay rather than to metabolic
status as a general category.

**Running total across the whole project:** starvation resistance (own lab,
R²=0.673, strong positive), starvation resistance (Morgante, R²=0.041, weak
positive, consistent with cross-lab measurement noise), 21 other independent
phenotype/diet tests showing no signal (fecundity, lifespan, chill coma,
cuticle HC, and 17 of 18 Unckless measures/diets), and 1 weak, unconfirmed
candidate (Unckless protein, low-glucose diet).

---

## 8. Environment

Python: pandas, numpy, scikit-learn, xgboost, seaborn, matplotlib.
(PLSRegression and PLS-DA are in scikit-learn: sklearn.cross_decomposition.)
R 4.4.2: rms, survival, emmeans, MASS, car, Rtsne, ggplot2.

---

## 9. Working norms (from the lab)

- FAIR principles: everything reproducible, version-controlled, documented.
- Keep code in git; commit incrementally with clear messages.
- This is assessed work. Every part of the pipeline must be understood and
  defensible by the author.

---

## 10. Key reference

Ibrahim R, González-Jiménez M, ..., Wynne K, Dobson AJ. "Barcoding biology:
Chemotype predicts variation in genotype, physiology, and stress response."
bioRxiv 2026. doi:10.64898/2026.03.22.713522.
