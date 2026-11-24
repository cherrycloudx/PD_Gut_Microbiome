root <- "results/machine_learning/PRJEB27564/experiments"

out <- file.path(root, "final_comparison")

dir.create(out, recursive = TRUE, showWarnings = FALSE)

cat("Searching for individual model summaries...\n")

files <- list.files(
  root,
  pattern = "^(ridge_model_summary|model_summary)\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No model summary files found.")
}

cat("Summary files found:", length(files), "\n")

results <- lapply(files, function(f) {
  x <- read.delim(
    f,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  x$Source_file <- substring(f, nchar(root) + 2)

  x
})

combined <- do.call(rbind, results)

expected_models <- expand.grid(
  Prevalence = c("5%", "10%", "20%", "30%"),
  Normalization = c("log.std", "log.clr"),
  Model = c("Ridge", "LASSO", "ElasticNet", "RandomForest"),
  stringsAsFactors = FALSE
)

make_key <- function(x) {
  paste(x$Prevalence, x$Normalization, x$Model, sep = "|")
}

observed_keys <- make_key(combined)
expected_keys <- make_key(expected_models)

duplicates <- unique(observed_keys[duplicated(observed_keys)])
missing <- expected_models[!expected_keys %in% observed_keys, ]

if (length(duplicates) > 0) {
  cat("\nDuplicate configurations found:\n")
  print(duplicates)
  stop("Resolve duplicates before creating the final table.")
}

if (nrow(missing) > 0) {
  cat("\nMissing configurations:\n")
  print(missing)
  stop("Some of the 32 expected configurations are missing.")
}

if (nrow(combined) != 32) {
  stop("Expected exactly 32 model summaries.")
}

cat("\nAll 32 unique configurations verified.\n")

required_columns <- c(
  "Model",
  "Prevalence",
  "Normalization",
  "Samples",
  "Genera",
  "Folds",
  "Repetitions",
  "Mean_repeat_AUROC",
  "SD_repeat_AUROC",
  "Mean_curve_AUROC",
  "Mean_repeat_AUPRC",
  "SD_repeat_AUPRC",
  "Mean_curve_AUPRC",
  "Source_file"
)

if (!all(required_columns %in% colnames(combined))) {
  stop("One or more expected summary columns are missing.")
}

combined$Prevalence_order <- match(
  combined$Prevalence,
  c("5%", "10%", "20%", "30%")
)

combined$Normalization_order <- match(
  combined$Normalization,
  c("log.std", "log.clr")
)

combined$Model_order <- match(
  combined$Model,
  c("Ridge", "LASSO", "ElasticNet", "RandomForest")
)

combined$Rank_AUROC <- rank(
  -combined$Mean_curve_AUROC,
  ties.method = "min"
)

combined$Rank_AUPRC <- rank(
  -combined$Mean_curve_AUPRC,
  ties.method = "min"
)

combined <- combined[
  order(
    combined$Prevalence_order,
    combined$Normalization_order,
    combined$Model_order
  ),
]

combined$Prevalence_order <- NULL
combined$Normalization_order <- NULL
combined$Model_order <- NULL

master_file <- file.path(
  out,
  "master_ML_model_comparison_32.tsv"
)

write.table(
  combined,
  master_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

by_auroc <- combined[
  order(-combined$Mean_curve_AUROC),
]

by_auprc <- combined[
  order(-combined$Mean_curve_AUPRC),
]

write.table(
  by_auroc,
  file.path(out, "ML_ranked_by_AUROC.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  by_auprc,
  file.path(out, "ML_ranked_by_AUPRC.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n====================================\n")
cat("MASTER ML RESULTS\n")
cat("====================================\n")

print(combined[, c(
  "Model",
  "Prevalence",
  "Normalization",
  "Genera",
  "Mean_curve_AUROC",
  "Mean_curve_AUPRC",
  "Rank_AUROC",
  "Rank_AUPRC"
)])

cat("\n====================================\n")
cat("TOP MODEL BY AUROC\n")
cat("====================================\n")

print(by_auroc[1, required_columns])

cat("\n====================================\n")
cat("TOP MODEL BY AUPRC\n")
cat("====================================\n")

print(by_auprc[1, required_columns])

cat("\nAll 32 configurations successfully collected.\n")
cat("Master table:", master_file, "\n")
cat("AUROC ranking:", file.path(out, "ML_ranked_by_AUROC.tsv"), "\n")
cat("AUPRC ranking:", file.path(out, "ML_ranked_by_AUPRC.tsv"), "\n")
