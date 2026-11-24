suppressPackageStartupMessages({
  library(ANCOMBC)
})

set.seed(56987)

count_file <- "results/taxonomy/PRJEB27564/abundance_matrices_131/genus_counts.tsv"

metadata_file <- "results/taxonomy/PRJEB27564/abundance_matrices_131/matched_sample_metadata.tsv"

out <- "results/statistics/PRJEB27564/ANCOMBC2_genus_nonconservative"

dir.create(
  out,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("ANCOMBC version:", as.character(packageVersion("ANCOMBC")), "\n")

cat("\nLoading count data and metadata...\n")

counts <- read.delim(count_file, check.names = FALSE)
metadata <- read.delim(metadata_file, check.names = FALSE)

if (!setequal(counts$SampleID, metadata$Run)) {
  stop("Sample identifiers do not match.")
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

if (sum(group == "no") != 64 ||
    sum(group == "yes") != 67) {
  stop("Unexpected group sizes.")
}

metadata$parkinson <- factor(
  group,
  levels = c("no", "yes")
)

cat("\nSample groups:\n")
print(table(metadata$parkinson))

genus_names <- setdiff(names(counts), "SampleID")

count_matrix <- t(
  as.matrix(counts[, genus_names, drop = FALSE])
)

storage.mode(count_matrix) <- "numeric"

rownames(count_matrix) <- genus_names
colnames(count_matrix) <- counts$SampleID

if (any(!is.finite(count_matrix)) ||
    any(count_matrix < 0) ||
    any(count_matrix != round(count_matrix))) {
  stop("Invalid count values.")
}

meta_data <- data.frame(
  parkinson = factor(group, levels = c("no", "yes")),
  ancombc_dummy = rep(1L, length(group)),
  row.names = metadata$Run,
  check.names = FALSE
)

stopifnot(
  identical(colnames(count_matrix), rownames(meta_data))
)

cat("\nCount matrix dimensions:\n")
print(dim(count_matrix))

cat("\nRunning secondary ANCOM-BC2 sensitivity analysis...\n")
cat("conservative = FALSE\n")
cat("pseudo_sens = TRUE\n")

fit <- ancombc2(
  data = count_matrix,
  taxa_are_rows = TRUE,
  meta_data = meta_data,
  fix_formula = "parkinson",
  rand_formula = NULL,
  p_adj_method = "BH",
  pseudo_sens = TRUE,
  conservative = FALSE,
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

cat("\nSaving secondary analysis results...\n")

saveRDS(
  fit,
  file = file.path(out, "ANCOMBC2_nonconservative_full_output.rds")
)

res <- as.data.frame(fit$res)

write.table(
  res,
  file = file.path(out, "ANCOMBC2_nonconservative_all_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (!is.null(fit$ss_tab)) {
  write.table(
    as.data.frame(fit$ss_tab),
    file = file.path(out, "ANCOMBC2_nonconservative_sensitivity.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

cat("\n====================================\n")
cat("NON-CONSERVATIVE ANCOM-BC2 SUMMARY\n")
cat("====================================\n")

cat("Genera tested:", nrow(res), "\n")

cat(
  "BH-significant:",
  sum(res$diff_parkinsonyes %in% TRUE),
  "\n"
)

cat(
  "Passed sensitivity:",
  sum(res$passed_ss_parkinsonyes %in% TRUE),
  "\n"
)

cat(
  "Significant and robust:",
  sum(res$diff_robust_parkinsonyes %in% TRUE),
  "\n"
)

cat("\nTop 20 results by adjusted p-value:\n")

res <- res[order(res$q_parkinsonyes), ]

print(
  head(
    res[, c(
      "taxon",
      "lfc_parkinsonyes",
      "p_parkinsonyes",
      "q_parkinsonyes",
      "diff_parkinsonyes",
      "passed_ss_parkinsonyes",
      "diff_robust_parkinsonyes"
    )],
    20
  ),
  row.names = FALSE
)

cat("\nSecondary analysis completed.\n")
cat("Results saved in:", out, "\n")
