
suppressPackageStartupMessages({
  library(SIAMCAT)
})

input_path <- "results/taxonomy/PRJEB27564/abundance_matrices_131"

output_path <- "results/machine_learning/PRJEB27564/input"

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

cat("SIAMCAT version:",
    as.character(packageVersion("SIAMCAT")), "\n")

# Load genus relative abundance data
feature_data <- read.delim(
  file.path(
    input_path,
    "genus_relative_abundance_assigned_reads.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Load sample metadata
meta <- read.delim(
  file.path(
    input_path,
    "matched_sample_metadata.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Validate sample IDs
feature_data$SampleID <- as.character(feature_data$SampleID)
meta$Run <- as.character(meta$Run)

if (anyDuplicated(feature_data$SampleID)) {
  stop("Duplicate sample IDs in feature table.")
}

if (anyDuplicated(meta$Run)) {
  stop("Duplicate sample IDs in metadata.")
}

# Convert feature table to a numeric matrix
sample_ids <- feature_data$SampleID

feature_data$SampleID <- NULL

feat <- t(as.matrix(feature_data))

storage.mode(feat) <- "double"

# Convert percentages to proportions
feat <- feat / 100

rownames(feat) <- colnames(feature_data)
colnames(feat) <- sample_ids

# Match metadata order to feature matrix
meta <- meta[match(colnames(feat), meta$Run), , drop = FALSE]

if (anyNA(meta$Run)) {
  stop("Some feature samples are missing from metadata.")
}

rownames(meta) <- meta$Run

# Validate disease labels
if (!all(meta$parkinson %in% c("yes", "no"))) {
  stop("Unexpected Parkinson's disease labels.")
}

meta$parkinson <- as.character(meta$parkinson)

# Validate feature matrix
if (anyNA(feat) || any(!is.finite(feat))) {
  stop("Missing or non-finite feature values detected.")
}

if (any(feat < 0)) {
  stop("Negative relative abundances detected.")
}

if (anyDuplicated(rownames(feat))) {
  stop("Duplicate genus feature names.")
}

if (!identical(colnames(feat), rownames(meta))) {
  stop("Feature and metadata sample order does not match.")
}

if (any(abs(colSums(feat) - 1) > 1e-4)) {
  stop("Relative abundance columns do not sum to 1.")
}

# Check cohort
cat("\nSamples:", ncol(feat), "\n")
cat("Genus features:", nrow(feat), "\n")

cat("\nDisease group distribution:\n")
print(table(meta$parkinson))

# Construct SIAMCAT object
cat("\nCreating SIAMCAT object...\n")

sc <- siamcat(
  feat = feat,
  meta = meta,
  label = "parkinson",
  case = "yes",
  verbose = 2
)

# Save input matrices and SIAMCAT object
saveRDS(
  feat,
  file.path(output_path, "genus_features_131.rds")
)

saveRDS(
  meta,
  file.path(output_path, "metadata_131.rds")
)

saveRDS(
  sc,
  file.path(output_path, "siamcat_genus_131.rds")
)

write.table(
  feat,
  file = file.path(output_path, "genus_features_131.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

cat("\nSIAMCAT input preparation completed.\n")
cat("Samples:", ncol(feat), "\n")
cat("Features:", nrow(feat), "\n")
cat("Positive class: PD (yes)\n")
cat("Control class: no\n")

cat("\nOutputs saved in:", output_path, "\n")
