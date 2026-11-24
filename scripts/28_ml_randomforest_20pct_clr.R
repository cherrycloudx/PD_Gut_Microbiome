suppressPackageStartupMessages({
  library(SIAMCAT)
})

set.seed(56987)

input_path <- "results/machine_learning/PRJEB27564/input"

ridge_path <- "results/machine_learning/PRJEB27564/experiments/20pct_LOG_Ridge"

clr_path <- "results/machine_learning/PRJEB27564/experiments/20pct_CLR"

log_path <- "results/machine_learning/PRJEB27564/experiments/20pct_LOG_comparison"

out <- file.path(clr_path, "RandomForest")

dir.create(out, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ranger", quietly = TRUE)) {
  stop("The ranger package is not installed.")
}

cat("Loading SIAMCAT input...\n")

sc <- readRDS(
  file.path(input_path, "siamcat_genus_131.rds")
)

cat("\nApplying 20% prevalence filtering...\n")

sc <- filter.features(
  sc,
  filter.method = "prevalence",
  cutoff = 0.20,
  rm.unmapped = FALSE,
  feature.type = "original",
  verbose = 2
)

filtered <- get.filt_feat.matrix(sc)

cat("Samples:", ncol(filtered), "\n")
cat("Genera:", nrow(filtered), "\n")

if (nrow(filtered) != 133 || ncol(filtered) != 131) {
  stop("Unexpected filtered feature dimensions.")
}

positive <- filtered[filtered > 0]

if (length(positive) == 0) {
  stop("No positive feature abundances.")
}

log.n0 <- min(positive) / 100

cat("\nApplying CLR normalization...\n")

sc <- normalize.features(
  sc,
  norm.method = "log.clr",
  norm.param = list(log.n0 = log.n0),
  feature.type = "filtered",
  verbose = 2
)

norm_features <- get.norm_feat.matrix(sc)

if (any(!is.finite(norm_features))) {
  stop("Non-finite values in CLR matrix.")
}

saveRDS(
  norm_features,
  file.path(clr_path, "normalized_features_CLR.rds")
)

cat("\nLoading 20% Ridge cross-validation splits...\n")

ridge_ref <- readRDS(
  file.path(ridge_path, "ridge_final_SIAMCAT.rds")
)

shared_split <- data_split(ridge_ref, verbose = 0)

if (is.null(shared_split) ||
    shared_split$num.folds != 10 ||
    shared_split$num.resample != 10) {
  stop("Expected 10-fold x 10-repeat validation splits.")
}

reference_samples <- colnames(get.filt_feat.matrix(ridge_ref))

if (!identical(colnames(filtered), reference_samples)) {
  stop("Sample order differs from the 20% Ridge reference.")
}

data_split(sc) <- shared_split

cat("Using the same 10-fold x 10-repetition splits.\n")
cat("Training Random Forest...\n")

sc <- train.model(
  sc,
  method = "randomForest",
  measure = "classif.acc",
  feature.type = "normalized",
  verbose = 2
)

saveRDS(
  sc,
  file.path(out, "randomforest_trained.rds")
)

cat("\nGenerating predictions...\n")

sc <- make.predictions(sc, verbose = 2)

cat("\nEvaluating predictions...\n")

sc <- evaluate.predictions(sc, verbose = 2)

evaluation <- eval_data(sc)

auc <- as.numeric(evaluation$auroc.all)
auprc <- as.numeric(evaluation$auprc.all)

if (length(auc) != 10 || length(auprc) != 10) {
  stop("Unexpected number of cross-validation results.")
}

metrics <- data.frame(
  Repeat = seq_along(auc),
  AUROC = auc,
  AUPRC = auprc
)

write.table(
  metrics,
  file = file.path(out, "metrics_by_repeat.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

model_summary <- data.frame(
  Model = "RandomForest",
  Prevalence = "20%",
  Normalization = "log.clr",
  Samples = 131,
  Genera = 133,
  Folds = 10,
  Repetitions = 10,
  Mean_repeat_AUROC = mean(auc),
  SD_repeat_AUROC = sd(auc),
  Mean_curve_AUROC = as.numeric(evaluation$auroc),
  Mean_repeat_AUPRC = mean(auprc),
  SD_repeat_AUPRC = sd(auprc),
  Mean_curve_AUPRC = as.numeric(evaluation$auprc)
)

write.table(
  model_summary,
  file = file.path(out, "model_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  sc,
  file.path(out, "randomforest_final_SIAMCAT.rds")
)

cat("\nUpdating 20% CLR algorithm comparison...\n")

clr_comparison <- read.delim(
  file.path(clr_path, "CLR_algorithm_comparison.tsv"),
  check.names = FALSE
)

clr_comparison <- rbind(
  clr_comparison,
  model_summary
)

clr_comparison <- clr_comparison[
  order(-clr_comparison$Mean_curve_AUROC),
]

write.table(
  clr_comparison,
  file.path(clr_path, "CLR_algorithm_comparison.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nCombining all eight 20% model configurations...\n")

log_comparison <- read.delim(
  file.path(log_path, "algorithm_comparison.tsv"),
  check.names = FALSE
)

all_comparison <- rbind(
  log_comparison,
  clr_comparison
)

all_comparison <- all_comparison[
  order(-all_comparison$Mean_curve_AUROC),
]

comparison_file <- "results/machine_learning/PRJEB27564/experiments/20pct_LOG_vs_CLR_all_algorithms.tsv"

write.table(
  all_comparison,
  comparison_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n20% CLR Random Forest completed.\n")

cat("\nRandom Forest summary:\n")
print(model_summary)

cat("\nFINAL 20% LOG vs CLR COMPARISON:\n")

print(all_comparison[, c(
  "Model",
  "Normalization",
  "Mean_curve_AUROC",
  "Mean_curve_AUPRC"
)])

cat("\nOutputs saved in:", out, "\n")
cat("\nCombined comparison saved in:", comparison_file, "\n")
