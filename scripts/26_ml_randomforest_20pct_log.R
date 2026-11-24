suppressPackageStartupMessages({
  library(SIAMCAT)
})

set.seed(56987)

input_path <- "results/machine_learning/PRJEB27564/input"

ridge_path <- "results/machine_learning/PRJEB27564/experiments/20pct_LOG_Ridge"

comparison_path <- "results/machine_learning/PRJEB27564/experiments/20pct_LOG_comparison"

out <- file.path(comparison_path, "RandomForest")

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
  stop("Unexpected feature matrix dimensions.")
}

positive <- filtered[filtered > 0]

if (length(positive) == 0) {
  stop("No positive feature abundances.")
}

log.n0 <- min(positive) / 100

cat("\nApplying log-standardisation...\n")

sc <- normalize.features(
  sc,
  norm.method = "log.std",
  norm.param = list(
    log.n0 = log.n0,
    sd.min.q = 0.1
  ),
  feature.type = "filtered",
  verbose = 2
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

cat("Using the same Ridge 10-fold x 10-repeat splits.\n")
cat("Samples:", ncol(filtered), "\n")
cat("Genera:", nrow(filtered), "\n")

cat("\nTraining Random Forest...\n")

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
  Normalization = "log.std",
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

cat("\nUpdating 20% algorithm comparison...\n")

comparison <- read.delim(
  file.path(comparison_path, "algorithm_comparison.tsv"),
  check.names = FALSE
)

comparison <- rbind(comparison, model_summary)

comparison <- comparison[
  order(-comparison$Mean_curve_AUROC),
]

write.table(
  comparison,
  file = file.path(comparison_path, "algorithm_comparison.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n20% Random Forest experiment completed.\n")

cat("\nRandom Forest summary:\n")
print(model_summary)

cat("\nFINAL 20% LOG ALGORITHM COMPARISON:\n")

print(comparison[, c(
  "Model",
  "Mean_curve_AUROC",
  "Mean_curve_AUPRC"
)])

cat("\nOutputs saved in:", out, "\n")
