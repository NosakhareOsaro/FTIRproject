#### explore_eatris_mae.R
####
#### Inspection-only script for the EATRIS-Plus multi-omics MultiAssayExperiment
#### (Zenodo DOI 10.5281/zenodo.17514796, 125 healthy adults). Reports the
#### structure needed to decide an extraction plan: does NOT write any
#### derived files.
####
#### Source files (not committed, download separately):
####   phenotype-data/raw/mae_mae.rds          - MultiAssayExperiment object
####   phenotype-data/raw/mae_experiments.h5   - HDF5-backed assay data
#### Both files must be kept together in the same directory: assay() access
#### on any experiment lazily reads from mae_experiments.h5 via a path stored
#### relative to R's working directory at access time (not relative to the
#### .rds file's own location), so this script temporarily changes into that
#### directory before touching any assay data.

suppressPackageStartupMessages(library(MultiAssayExperiment))

RAW_DIR <- file.path("phenotype-data", "raw")
RDS_FILE <- "mae_mae.rds"
H5_FILE <- "mae_experiments.h5"

stopifnot(file.exists(file.path(RAW_DIR, RDS_FILE)))
stopifnot(file.exists(file.path(RAW_DIR, H5_FILE)))

old_wd <- getwd()
setwd(RAW_DIR)
on.exit(setwd(old_wd), add = TRUE)

mae <- readRDS(RDS_FILE)

cat("========================================================================\n")
cat("1. Experiments/assays in the MultiAssayExperiment object\n")
cat("========================================================================\n")
print(mae)

cat("\n")
cat("========================================================================\n")
cat("2. Dimensions of each assay (features x samples)\n")
cat("========================================================================\n")
exp_list <- experiments(mae)
dims <- t(sapply(names(exp_list), function(nm) dim(exp_list[[nm]])))
colnames(dims) <- c("n_features", "n_samples")
print(dims)

cat("\n")
cat("========================================================================\n")
cat("3. Phenotype / colData table\n")
cat("========================================================================\n")
cd <- colData(mae)
cat("colData dimensions (samples x columns):", nrow(cd), "x", ncol(cd), "\n\n")
cat("Full column list:\n")
print(colnames(cd))

cat("\nBMI present:", "BMI" %in% colnames(cd), "\n")
if ("BMI" %in% colnames(cd)) {
  cat("BMI class:", class(cd$BMI), "\n")
  cat("BMI non-missing:", sum(!is.na(cd$BMI)), "of", nrow(cd), "samples\n")
  cat("BMI range:", paste(range(cd$BMI, na.rm = TRUE), collapse = " - "), "\n")
  cat("BMI mean (SD):", round(mean(cd$BMI, na.rm = TRUE), 2),
      "(", round(sd(cd$BMI, na.rm = TRUE), 2), ")\n")
}
if ("BMI.group" %in% colnames(cd)) {
  cat("\nBMI.group table:\n")
  print(table(cd$BMI.group, useNA = "ifany"))
}
if ("Sex" %in% colnames(cd)) {
  cat("\nSex table:\n")
  print(table(cd$Sex, useNA = "ifany"))
}
if ("Age" %in% colnames(cd)) {
  cat("\nAge non-missing:", sum(!is.na(cd$Age)), "of", nrow(cd),
      "; range:", paste(range(cd$Age, na.rm = TRUE), collapse = " - "), "\n")
}

cat("\n")
cat("========================================================================\n")
cat("4. Sample ID overlap: lipidomics assays vs. phenotype table\n")
cat("========================================================================\n")

pheno_ids <- rownames(cd)
bmi_complete_ids <- rownames(cd)[!is.na(cd$BMI)]

lipid_assay_names <- grep("Lipidomics", names(exp_list), value = TRUE)
cat("Lipidomics assay(s) found:", paste(lipid_assay_names, collapse = "; "), "\n\n")

for (nm in lipid_assay_names) {
  ex <- exp_list[[nm]]
  ids <- colnames(ex)
  cat("---", nm, "---\n")
  cat("  n samples in assay          :", length(ids), "\n")
  cat("  n overlapping with colData  :", length(intersect(ids, pheno_ids)), "\n")
  cat("  n overlapping with BMI-complete samples:",
      length(intersect(ids, bmi_complete_ids)), "\n")
  missing_pheno <- setdiff(ids, pheno_ids)
  if (length(missing_pheno) > 0) {
    cat("  assay samples with NO phenotype row:", paste(missing_pheno, collapse = ", "), "\n")
  }
}

if (length(lipid_assay_names) > 1) {
  cat("\nPairwise overlap between lipidomics assays:\n")
  for (i in seq_along(lipid_assay_names)) {
    for (j in seq_along(lipid_assay_names)) {
      if (i < j) {
        a <- colnames(exp_list[[lipid_assay_names[i]]])
        b <- colnames(exp_list[[lipid_assay_names[j]]])
        cat("  ", lipid_assay_names[i], "&", lipid_assay_names[j], ":",
            length(intersect(a, b)), "shared /", length(union(a, b)), "total\n")
        only_a <- setdiff(a, b)
        only_b <- setdiff(b, a)
        if (length(only_a) > 0) cat("    only in", lipid_assay_names[i], ":", paste(only_a, collapse = ", "), "\n")
        if (length(only_b) > 0) cat("    only in", lipid_assay_names[j], ":", paste(only_b, collapse = ", "), "\n")
      }
    }
  }
}

cat("\nsampleMap sanity check: any assay where primary != colname (ID remapping)?\n")
sm <- sampleMap(mae)
mismatches <- sum(as.character(sm$primary) != as.character(sm$colname))
cat("  Mismatches across all 15 experiments:", mismatches, "\n")

cat("\n")
cat("========================================================================\n")
cat("5. rowData preview for each lipidomics assay (feature ID / annotation columns)\n")
cat("========================================================================\n")
for (nm in lipid_assay_names) {
  ex <- exp_list[[nm]]
  cat("---", nm, "---\n")
  cat("  rowData columns:", paste(colnames(rowData(ex)), collapse = ", "), "\n")
  cat("  first 3 feature IDs:", paste(head(rownames(ex), 3), collapse = ", "), "\n\n")
}

cat("Done. No files written by this script.\n")
