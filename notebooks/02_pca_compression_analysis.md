# 02 — PCA Compression Analysis

**Script:** `scripts/run_compression_analysis.py`
**Status:** Complete
**Key result:** Line-mean PCA PC1 correlates r = +0.685 (p = 2.92×10⁻¹⁶) with starvation EMMeans. Individual-fly PCA gives r = −0.05 (p = 0.61) — no signal.

---

## What this step does and why

The pre-print established that FTIR spectra contain enough information to classify flies as starvation-resistant or sensitive. The next question is whether this information is structured enough to predict a continuous starvation resistance score — not just a binary label.

Before building any prediction model, it helps to first ask a simpler question: is there any visible spectral structure that aligns with starvation resistance, without using the phenotype during compression? Principal Component Analysis (PCA) answers this. PCA finds the directions of maximum variance in the spectral data, completely unsupervised — it has no knowledge of which flies come from resistant or sensitive lines. If PC1 (the direction of greatest spectral variation) correlates with starvation resistance, that is strong evidence the spectra encode the biological signal.

This step also resolves a practical question that shapes every downstream analysis: at what level should the data be collapsed? Should spectra be averaged to DGRP line means before analysis, or should individual fly spectra be used throughout?

---

## The data

`DGRPFTIR.dat` contains 1,772 female spectra from 108 DGRP lines, approximately 16 flies per line. Each spectrum has 1,723 absorbance values at wavenumbers from 3,900 to 456 cm⁻¹. The starvation phenotype target is the EMMean from `Emmeans.csv` — one continuous score per DGRP line (not per fly).

Because the starvation assay and the FTIR measurement are both destructive, they were run on sibling flies from the same lines, not on the same individual flies. The spectra and phenotypes can only be linked at the line level, not the fly level.

---

## Step 1 — How many principal components are needed?

Before examining the correlation with phenotype, PCA was run on all 1,772 individual fly spectra (after StandardScaling each wavenumber to mean 0, standard deviation 1).

**Result: only 4 PCs capture 95.3% of all spectral variance.**

- PC1: 53.6% of variance
- PC2: 32.2% of variance
- PC3, PC4: remainder to 95.3%

FTIR spectra are highly redundant — adjacent wavenumbers absorb similar compounds and are strongly correlated. Despite having 1,723 measurements per fly, the underlying spectral variation is almost entirely captured by 4 dimensions. This is important because it means any compression method (PCA, PLS, or regularised regression) has a much simpler problem than the raw dimensionality suggests.

---

## Step 2 — The collapse-order question

This is the most important methodological finding of this step.

Starvation EMMeans exist at the line level, not the fly level. To compare PC scores with EMMeans, the fly-level PC scores have to be collapsed to line means at some point. There are two natural orderings:

**Order A — PCA first, then average:**
Run PCA on all 1,772 individual fly spectra. Project each fly into PC space. Average the PC scores within each DGRP line to get one line-level score.

**Order B — Average first, then PCA:**
Average the raw spectra within each DGRP line first (108 lines × 1,723 wavenumbers). Run PCA on the 108 line-mean spectra.

The question is whether the order matters. It does — dramatically.

| Order | Pearson r (PC1 vs EMMean) | p-value | Signal? |
|---|---|---|---|
| A — PCA first, average after | −0.050 | 0.606 | No |
| B — Average first, PCA after | **+0.685** | **2.92×10⁻¹⁶** | **Yes** |

**Why the difference is so large:**

In Order A, PCA is run on 1,772 individual fly spectra. The dominant source of variation across those spectra is within-line fly-to-fly noise — random measurement differences between sibling flies from the same genotype. This noise is large relative to the between-line genetic signal. PC1 captures this noise, not the genetic signal, so it shows no correlation with starvation resistance.

In Order B, averaging ~16 spectra per line before PCA applies the law of large numbers: random within-line noise cancels out, and what remains is predominantly between-line genetic variation. PCA on this averaged matrix finds spectral directions that separate genotypes from each other, and PC1 turns out to align strongly with starvation resistance.

**What this means for the biology:** The starvation resistance signal in FTIR spectra is a property of the genotype, not something reliably readable from an individual fly. A single fly's spectrum contains substantial measurement noise that drowns out the genetic signal. Averaging across flies from the same inbred line reveals the genotypic chemical fingerprint.

**What this means for the pipeline:** Order B is the correct collapse order. All subsequent line-mean analyses use this approach: average spectra within each DGRP line first, then fit models on the 108 line-mean spectra.

---

## Step 3 — Visualising the signal

Two scatter plots were produced to illustrate the contrast between Order A and Order B.

**`results/DGRP/pca_coloured_by_emmean.pdf` (Order A — individual flies):**
PC1 vs PC2, coloured by each fly's line EMMean. The colours are randomly scattered across the plot — no visible gradient. Confirms r = −0.05.

**`results/DGRP/pca_linemeans_coloured_by_emmean.pdf` (Order B — line means):**
PC1 vs PC2 for the 108 line-mean spectra, coloured by EMMean. A clear colour gradient runs left to right along PC1: dark purple (sensitive, low EMMean) on the left, orange and yellow (resistant, high EMMean) on the right. r = +0.685, p = 2.92×10⁻¹⁶ is annotated directly on the figure.

This plot is the clearest single visualisation of the main finding: the FTIR spectra of Drosophila genotypes encode starvation resistance information that becomes visible only after within-line averaging.

---

## Technical details

- **StandardScaler** was applied before PCA in all cases, fitted on the full dataset (no CV at this stage — this step is exploratory visualisation, not predictive modelling)
- **n_components = 100** was fitted; the 95% threshold was reached at PC4
- **The sign of PC1** is arbitrary in PCA — it can be positive or negative depending on initialisation. The magnitude of the correlation (0.685) is what matters, not the sign
- All code is in `scripts/run_compression_analysis.py`, using `sklearn.decomposition.PCA`

---

## Output files

- `results/DGRP/pca_explained_variance.pdf` — cumulative variance curve, 95% threshold line at PC4
- `results/DGRP/pca_coloured_by_emmean.pdf` — individual-fly PC scatter (no signal, Order A)
- `results/DGRP/pca_linemeans_coloured_by_emmean.pdf` — line-mean PC scatter (r = +0.685, Order B)