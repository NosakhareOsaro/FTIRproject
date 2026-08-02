#### prepare_eatris_lipidomics.R
####
#### Extracts three lipidomics feature matrices - positive mode, negative
#### mode, and a combined positive+negative matrix - plus the BMI phenotype
#### vector and a small Sex/Age/BMI covariates file (used only for the
#### demographic-confound check, not for the regression pipeline itself),
#### from the EATRIS-Plus MultiAssayExperiment (Zenodo DOI
#### 10.5281/zenodo.17514796, 125 healthy human adults). Writes plain TSV
#### files so the Python pipeline can read them with no R/Bioconductor
#### dependency downstream.
####
#### Treated as three separate test conditions, the same way
#### prepare_unckless_data.py produced pooled/high-glucose/low-glucose
#### variants of the same underlying comparison:
####   positive only   (196 features)
####   negative only   (164 features)
####   combined        (360 features, inner-joined on the 125 shared samples;
####                    feature columns prefixed pos_/neg_ to keep the two
####                    modes' features unambiguous even though the raw IDs
####                    never actually collide - checked directly, 0 overlap)
####
#### Source files (not committed to git, download from Zenodo separately,
#### see scripts/explore_eatris_mae.R for the inspection that established
#### this structure):
####   phenotype-data/raw/mae_mae.rds
####   phenotype-data/raw/mae_experiments.h5
####
#### Both files are required. Assay data is HDF5-backed and lazily read from
#### the .h5 file at access time, via a path stored *inside* the .rds object
#### and resolved relative to R's working directory when the data is
#### touched - not relative to the .rds file's own location.
#### explore_eatris_mae.R worked around this with a temporary setwd(), which
#### is fragile: it silently breaks if the script is invoked from somewhere
#### unexpected, or if anything else changes the working directory
#### mid-session. Here the HDF5 seed's stored path is instead overwritten
#### with an absolute path resolved via here() *before* any data is read, so
#### this script never depends on the working directory at all.

suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(HDF5Array)
  library(here)
})

RAW_DIR <- here("phenotype-data", "raw")
OUT_DIR <- here("phenotype-data")
RDS_PATH <- file.path(RAW_DIR, "mae_mae.rds")
H5_PATH <- normalizePath(file.path(RAW_DIR, "mae_experiments.h5"), mustWork = TRUE)

stopifnot(file.exists(RDS_PATH))

mae <- readRDS(RDS_PATH)

#' Read a SummarizedExperiment's assay as a dense samples x features matrix.
#' Rewrites the HDF5 seed's path to the absolute H5_PATH (computed above via
#' here()) before touching any data, so this never depends on getwd().
extract_samples_by_features <- function(se) {
  a <- assay(se, withDimnames = FALSE)
  path(a) <- H5_PATH
  m <- as.matrix(a)
  dimnames(m) <- dimnames(se)
  t(m)
}

exp_list <- experiments(mae)
pos_se <- exp_list[["Lipidomics, positive | transformed"]]
neg_se <- exp_list[["Lipidomics, negative | transformed"]]

pos_mat <- extract_samples_by_features(pos_se)  # samples x 196 features
neg_mat <- extract_samples_by_features(neg_se)  # samples x 164 features

stopifnot(
  "positive and negative feature IDs must not collide" =
    length(intersect(colnames(pos_mat), colnames(neg_mat))) == 0
)

shared_samples <- intersect(rownames(pos_mat), rownames(neg_mat))
stopifnot(
  "expected identical 125-sample sets for positive and negative assays" =
    length(shared_samples) == nrow(pos_mat) && length(shared_samples) == nrow(neg_mat)
)
pos_mat <- pos_mat[shared_samples, , drop = FALSE]
neg_mat <- neg_mat[shared_samples, , drop = FALSE]

pos_for_combo <- pos_mat
colnames(pos_for_combo) <- paste0("pos_", colnames(pos_for_combo))
neg_for_combo <- neg_mat
colnames(neg_for_combo) <- paste0("neg_", colnames(neg_for_combo))
combined_mat <- cbind(pos_for_combo, neg_for_combo)

# ---- BMI phenotype, matched to the same sample set ------------------------
cd <- colData(mae)
bmi_df <- data.frame(
  sample_id = rownames(cd),
  BMI = as.numeric(cd$BMI)
)
bmi_df <- bmi_df[bmi_df$sample_id %in% shared_samples, ]
bmi_df <- bmi_df[match(shared_samples, bmi_df$sample_id), ]

# ---- Sex/Age covariates, for the demographic-confound check only ----------
# Not used by the regression pipeline - written so the check in notebook 10
# (does the lipidomics signal just recover Sex/Age?) has a reproducible
# source file instead of a one-off interactive extraction.
cov_df <- data.frame(
  sample_id = rownames(cd),
  Sex = cd$Sex,
  Age = as.numeric(cd$Age),
  BMI = as.numeric(cd$BMI)
)
cov_df <- cov_df[cov_df$sample_id %in% shared_samples, ]
cov_df <- cov_df[match(shared_samples, cov_df$sample_id), ]

# ---- Write outputs ----------------------------------------------------------
write_matrix_tsv <- function(mat, path) {
  out <- data.frame(sample_id = rownames(mat), mat, check.names = FALSE)
  write.table(out, path, sep = "\t", row.names = FALSE, quote = FALSE)
}

write_matrix_tsv(pos_mat, file.path(OUT_DIR, "EATRIS_Lipidomics_positive.tsv"))
write_matrix_tsv(neg_mat, file.path(OUT_DIR, "EATRIS_Lipidomics_negative.tsv"))
write_matrix_tsv(combined_mat, file.path(OUT_DIR, "EATRIS_Lipidomics_combined.tsv"))
write.table(bmi_df, file.path(OUT_DIR, "EATRIS_BMI.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cov_df, file.path(OUT_DIR, "EATRIS_covariates.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---- Summary ----------------------------------------------------------------
cat("========================================================================\n")
cat("Output summary\n")
cat("========================================================================\n")
summarise_one <- function(name, mat) {
  cat(sprintf("%-32s %4d samples x %5d features   NAs: %d\n",
              name, nrow(mat), ncol(mat), sum(is.na(mat))))
}
summarise_one("EATRIS_Lipidomics_positive.tsv", pos_mat)
summarise_one("EATRIS_Lipidomics_negative.tsv", neg_mat)
summarise_one("EATRIS_Lipidomics_combined.tsv", combined_mat)
cat(sprintf("%-32s %4d samples,   BMI non-missing: %d, range %g - %g\n",
            "EATRIS_BMI.tsv", nrow(bmi_df), sum(!is.na(bmi_df$BMI)),
            min(bmi_df$BMI, na.rm = TRUE), max(bmi_df$BMI, na.rm = TRUE)))
cat(sprintf("%-32s %4d samples,   Sex/Age/BMI, for the confound check only\n",
            "EATRIS_covariates.tsv", nrow(cov_df)))
cat("\nSample ID order identical across all five files:",
    identical(rownames(pos_mat), rownames(neg_mat)) &&
      identical(rownames(pos_mat), rownames(combined_mat)) &&
      identical(rownames(pos_mat), bmi_df$sample_id) &&
      identical(rownames(pos_mat), cov_df$sample_id), "\n")
