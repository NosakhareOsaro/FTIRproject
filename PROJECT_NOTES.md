# PROJECT_NOTES.md

Context file for working on this project locally (e.g. with Claude Code).
Drop this in the repo root. It summarises the agreed plan, data status, and
evaluation design so any assistant starting fresh has the full picture.

**Author:** Nosakhare Odionfo Osaro (MSc Bioinformatics, University of Glasgow)
**Supervisor:** Dr Adam Dobson
**Last updated:** 24 June 2026

---

## 1. What this project is

MSc dissertation. Working title: _Machine Learning Signatures of Microbiome
Perturbation_ (the official title is broader than the actual work — see below).

**Actual focus:** comprehensive exploration of multivariate methods for
phenotype prediction from FTIR spectra in _Drosophila melanogaster_. The
core question is: can we generate continuous vectors from spectra that
reliably predict phenotype in individuals not seen during training? This
builds directly on a lab pre-print (Ibrahim et al., bioRxiv 2026,
doi:10.64898/2026.03.22.713522) — "chemotyping" — which uses ATR-FTIR + ML
to classify biological variation and predict starvation resistance.

**Framing:** the project is now "explore the method space comprehensively"
rather than "test a specific hypothesis." The FTIR/DGRP work is the primary
case study. Methods span dimensionality reduction (PCA, tSNE), methods that
do both reduction and prediction (PLS-DA, sPLS-DA), and direct prediction
(LASSO, ridge regression, elastic net, random forest). Regularised/frequentist
methods (PLS-DA, elastic net, LASSO) serve as baselines within the broader
comparison — the frequentist angle is not abandoned, just embedded in a wider
sweep. Common evaluation yardstick: held-out prediction performance (R² for
continuous targets, accuracy for discrete). A set of external clinical cohorts
(NoMIC, EATRIS-Plus, ROP, IBD, IMPACC, CMAISE) are referenced as motivation
only, not worked through — unless time allows after the FTIR work.

---

## 2. The gap this project fills

The pre-print has two limitations that this project addresses:

1. **No frequentist baseline.** It compares five ML classifiers (LR, SVM, kNN,
   RF, XGBoost) only to each other. It does NOT include PLS-DA or PLS
   regression — the field-standard chemometric baseline. Given the strong
   colinearity in spectral data, PLS is the natural comparator and its absence
   is a real gap. (Supervisor explicitly endorsed adding it.)

2. **Continuous outcomes collapsed to categories.** Every continuous or ordinal
   outcome is discretised before ML: starvation resistance (continuous EMMeans
   → 20th/80th percentile binary), dietary restriction (5/10/20% yeast →
   3-class), age (9d/4wk/6wk → 3-class). The paper's own thesis — "chemical
   state encodes biological state" — is fundamentally a continuous claim, so
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
   measurement). The DGRP panel has extensive public phenotype data —
   https://dgrpool.epfl.ch/ — so the same spectra can be tested against many
   continuous targets. Starvation resistance remains the positive control
   (signal is confirmed to be there).

---

## 4. CRITICAL: evaluation design (do not get this wrong)

Both the FTIR assay and the starvation assay are **destructive** — the same
fly cannot be measured for both. So FTIR spectra and phenotype can only be
linked at the **line (genotype) level**, not the individual fly level.

Consequences for the pipeline:

- **Training unit:** individual fly (~16 flies per line × 108 lines ≈ 1,700
  spectra). This preserves within-line spectral variation as signal and gives
  a much larger training set than 108 line-means.
- **Cross-validation:** MUST be **line-stratified** (group k-fold by DGRP line).
  No fly from a given line may appear in both train and test folds. Random
  fly-level CV would leak the line-level target across folds and inflate
  performance — this is the single most important correctness requirement.
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
pre-print figures — they match.

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
columns running from 3900 down to ~456 cm⁻¹. Column count varies (1727 vs 1728) depending on whether `StoTime` is present — the loader must handle both.

**Metadata coding is inconsistent and needs normalising** in the loader:

- Diet: `D05`/`D10`/`D20` (yeast %), `LI`/`SY` (lipid vs standard yeast)
- Age: `09D`, `04W`, `06W`, `03W`, `UNK`
- Sex: `F`, `M`, and `S` (`S` in Diet1 = confirmed female by data owner); some `UNK`/`SY` placeholders

Phenotype data:

- `Survival-data/DGRP-starvationresistance.csv` — individual fly survival
  events (2,887 rows) → fed into the R survival model.
- `Survival-data/Diet2-lifespan.csv`, `Fecundity-data/Diet2-fecundity.csv` —
  other phenotypes.

---

## 6. Existing code in the repo

- `scripts/DGRP_survival_analysis.R` — fits a parametric survival model
  (rms::psm, logistic dist) across the 108 DGRP lines, extracts EMMeans per
  line (the continuous starvation-resistance coefficient → `Emmeans.csv`),
  then bins at 20th/80th percentiles into sensitive/resistant. Also runs tSNE.
  **The EMMeans this produces are the regression target.**
- `scripts/FTIR-script.ipynb` — one big `MIRSPIPELINE()` function (by Rita
  Ibrahim & Mario Gonzalez Jimenez). Handles all experiments. XGBoost for
  feature (wavenumber) selection, then trains 5 classifiers with GridSearchCV,
  evaluates with stratified k=20 CV. All CLASSIFIERS — no regressors.
- Upstream QC uses the `bad-blood` package
  (github.com/magonji/bad-blood) — discards low-intensity spectra, atmospheric
  interference, etc. CHECK whether the .dat files are already QC'd or whether
  this needs running first.

---

## 7. Work order (updated 24 Jun 2026)

Steps 1–3 are done or in progress; steps 4–6 are the active frontier.

1. **[DONE] Reproduce Rita's classification baseline** on `DGRPFTIR.dat`.
   Recovered ~88.5% resistant / ~83.9% sensitive SVM reclassification
   (paper: 88% / 84%). Gate passed.
2. **[DONE] Write a clean data loader** (`scripts/ftir_loader.py`) — reads
   any of the 6 .dat files, normalises messy metadata, validates wavenumber
   axis.
3. **[DONE] Reproduce the R survival analysis** to regenerate `Emmeans.csv`
   (+ SEs). EMMeans match Morgante 2015 (r=0.46, 93/108 lines overlap).
4. **Build the method-comparison pipeline.** Implement and evaluate all methods
   from §3 on `DGRPFTIR.dat` with starvation EMMeans as target. Two evaluation
   settings must both be reported:
   - **(a) Per-fly with line-stratified CV** (see §4 for correctness
     requirements — this is the most important).
   - **(b) Line-mean spectra** (collapse first, then fit and evaluate on 108
     line means; simpler but lower-N).
   Evaluation metric: held-out R², RMSE, Spearman ρ for continuous targets;
   accuracy for any discrete targets. PLS-DA, elastic net, and LASSO are the
   regularised baselines. SVR is the like-for-like ML comparator to the
   pre-print's SVM classifier.
5. **Extend to other DGRPool phenotypes** — produce a shortlist of
   well-measured continuous phenotypes, run the same pipeline, assess whether
   spectral signal generalises beyond starvation resistance.
6. **Meeting with Vinny Davies** (mathematician, potential collaborator) —
   Monday 29 June 2026, 2 pm. No fixed agenda; may shape the mathematical
   approach, particularly around how to compare dimensionality-reduction
   methods with direct prediction methods, and whether Bayesian/GP approaches
   are worth adding given small N.

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
  defensible by the author — use AI assistance to accelerate and learn, not to
  outsource understanding. Be transparent with the supervisor about tool use,
  and check programme rules on AI use in assessed work.

---

## 10. Key reference

Ibrahim R, González-Jiménez M, ..., Wynne K, Dobson AJ. "Barcoding biology:
Chemotype predicts variation in genotype, physiology, and stress response."
bioRxiv 2026. doi:10.64898/2026.03.22.713522.
