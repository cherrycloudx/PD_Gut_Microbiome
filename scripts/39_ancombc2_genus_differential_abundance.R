suppressPackageStartupMessages({
  library(ANCOMBC)
})

set.seed(56987)

count_file <- "results/taxonomy/PRJEB27564/abundance_matrices_131/genus_counts.tsv"

metadata_file <- "results/taxonomy/PRJEB27564/abundance_matrices_131/matched_sample_metadata.tsv"

out <- "results/statistics/PRJEB27564/ANCOMBC2_genus"

dir.create(
  out,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Loading genus count data...\n")

counts <- read.delim(
  count_file,
  check.names = FALSE
)

metadata <- read.delim(
  metadata_file,
  check.names = FALSE
)

if (!all(c("SampleID") %in% names(counts))) {
  stop("SampleID missing from count data.")
}

if (!all(c("Run", "parkinson") %in% names(metadata))) {
  stop("Required metadata columns missing.")
}

if (anyDuplicated(counts$SampleID) ||
    anyDuplicated(metadata$Run)) {
  stop("Duplicate sample identifiers found.")
}

if (!setequal(counts$SampleID, metadata$Run)) {
  stop("Count and metadata sample IDs do not match.")
}

metadata <- metadata[
  match(counts$SampleID, metadata$Run),
]

if (!identical(counts$SampleID, metadata$Run)) {
  stop("Sample alignment failed.")
}

group <- tolower(trimws(as.character(metadata$parkinson)))

if (!all(group %in% c("no", "yes"))) {
  stop("Unexpected Parkinson status values.")
}

metadata$parkinson <- factor(
  group,
  levels = c("no", "yes")
)

cat("\nSample groups:\n")
print(table(metadata$parkinson))

if (sum(group == "yes") != 67 ||
    sum(group == "no") != 64) {
  stop("Unexpected group sizes.")
}

cat("\nPreparing taxa-by-sample count matrix...\n")

genus_names <- setdiff(
  names(counts),
  "SampleID"
)

count_matrix <- t(
  as.matrix(counts[, genus_names, drop = FALSE])
)

storage.mode(count_matrix) <- "numeric"

rownames(count_matrix) <- genus_names
colnames(count_matrix) <- counts$SampleID

if (any(!is.finite(count_matrix)) ||
    any(count_matrix < 0) ||
    any(count_matrix != round(count_matrix))) {
  stop("Count matrix contains invalid count values.")
}

cat("Count matrix dimensions:\n")
print(dim(count_matrix))

cat("\nLibrary size summary:\n")
print(summary(colSums(count_matrix)))

if (any(colSums(count_matrix) < 1000)) {
  stop("At least one sample falls below the specified library-size threshold.")
}

meta_data <- data.frame(
  parkinson = factor(group, levels = c("no", "yes")),
  ancombc_dummy = rep(1L, length(group)),
  row.names = metadata$Run,
  check.names = FALSE
)
stopifnot(
  ncol(meta_data) == 2,
  identical(rownames(meta_data), colnames(count_matrix)),
  identical(as.character(meta_data$parkinson), group)
)

cat("\nValidated metadata columns:\n")
print(names(meta_data))

cat("\nParkinson status:\n")
print(table(meta_data$parkinson))

if (!identical(colnames(count_matrix), rownames(meta_data))) {
  stop("Metadata row names do not match count matrix sample names.")
}

cat("\nSaving validated input files...\n")

write.table(
  count_matrix,
  file = file.path(out, "ANCOMBC2_input_genus_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  meta_data,
  file = file.path(out, "ANCOMBC2_input_metadata.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

cat("\nRunning ANCOM-BC2...\n")

fit <- ancombc2(
  data = count_matrix,
  taxa_are_rows = TRUE,
  meta_data = meta_data,
  fix_formula = "parkinson",
  rand_formula = NULL,
  p_adj_method = "BH",
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 1000,
  s0_perc = 0.05,
  group = "parkinson",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 1,
  verbose = TRUE,
  global = FALSE,
  pairwise = FALSE,
  dunnet = FALSE,
  trend = FALSE
)

cat("\nSaving complete ANCOM-BC2 output...\n")

saveRDS(
  fit,
  file = file.path(out, "ANCOMBC2_full_output.rds")
)

if (is.null(fit$res)) {
  stop("ANCOM-BC2 did not return primary results.")
}

res <- as.data.frame(fit$res)

write.table(
  res,
  file = file.path(out, "ANCOMBC2_all_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (!is.null(fit$zero_ind)) {
  write.table(
    as.data.frame(fit$zero_ind),
    file = file.path(out, "ANCOMBC2_structural_zeros.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

cat("\n====================================\n")
cat("ANCOM-BC2 RESULTS\n")
cat("====================================\n")

cat("Genera in primary results:", nrow(res), "\n")

cat("\nAvailable result columns:\n")
print(names(res))

cat("\nFirst 15 result rows:\n")
print(head(res, 15))

cat("\nResults directory:", out, "\n")

cat("\nANCOM-BC2 analysis completed.\n")
