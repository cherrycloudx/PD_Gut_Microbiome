suppressPackageStartupMessages({
  library(SIAMCAT)
})

set.seed(56987)

input_file <- "results/machine_learning/PRJEB27564/input/siamcat_genus_131.rds"

reference_file <- "results/machine_learning/PRJEB27564/experiments/10pct_LOG_Ridge/ridge_final_SIAMCAT.rds"

out <- "results/machine_learning/PRJEB27564/experiments/20pct_LOG_Ridge"

dir.create(out, recursive = TRUE, showWarnings = FALSE)

required <- c("mlr3", "mlr3learners", "glmnet")

missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  stop("Missing required packages: ", paste(missing, collapse = ", "))
}

cat("Loading SIAMCAT input...\n")
sc <- readRDS(input_file)

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
cat("Genera retained:", nrow(filtered), "\n")

if (ncol(filtered) != 131) {
  stop("Unexpected sample count.")
}

if (nrow(filtered) != 133) {
  stop("Expected 133 genera. Inspect filtering before modelling.")
}

write.table(
  filtered,
  file = file.path(out, "filtered_genus_features.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

positive_values <- filtered[filtered > 0]

if (length(positive_values) == 0) {
  stop("No positive feature abundances.")
}

log.n0 <- min(positive_values) / 100

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

norm_features <- get.norm_feat.matrix(sc)

if (any(!is.finite(norm_features))) {
  stop("Non-finite values in normalized feature matrix.")
}

saveRDS(
  norm_features,
  file.path(out, "normalized_features_LOG.rds")
)

cat("\nLoading shared cross-validation splits...\n")

ridge_ref <- readRDS(reference_file)

shared_split <- data_split(ridge_ref, verbose = 0)

if (is.null(shared_split) ||
    shared_split$num.folds != 10 ||
    shared_split$num.resample != 10) {
  stop("Could not retrieve the expected 10 x 10 splits.")
}

reference_samples <- colnames(get.filt_feat.matrix(ridge_ref))

if (!identical(colnames(filtered), reference_samples)) {
  stop("Sample order differs from the 10% Ridge reference.")
}

data_split(sc) <- shared_split

cat("Using the same 10-fold x 10-repetition splits.\n")

saveRDS(
  sc,
  file.path(out, "preprocessed_SIAMCAT.rds")
)

cat("\nTraining Ridge regression...\n")

sc <- train.model(
  sc,
  method = "ridge",
  measure = "classif.acc",
  feature.type = "normalized",
  verbose = 2
)

saveRDS(
  sc,
  file.path(out, "ridge_trained_SIAMCAT.rds")
)

cat("\nGenerating predictions...\n")

sc <- make.predictions(sc, verbose = 2)

cat("\nEvaluating model...\n")

sc <- evaluate.predictions(sc, verbose = 2)

evaluation <- eval_data(sc)

if (is.null(evaluation$auroc.all)) {
  stop("Repeated cross-validation AUC results are missing.")
}

auc <- as.numeric(evaluation$auroc.all)
auprc <- as.numeric(evaluation$auprc.all)

if (length(auc) != 10 || length(auprc) != 10) {
  stop("Expected metrics from 10 repetitions.")
}

repeat_results <- data.frame(
  Repeat = seq_along(auc),
  AUROC = auc,
  AUPRC = auprc
)

write.table(
  repeat_results,
  file = file.path(out, "ridge_metrics_by_repeat.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

model_summary <- data.frame(
  Model = "Ridge",
  Prevalence = "20%",
  Normalization = "log.std",
  Samples = ncol(filtered),
  Genera = nrow(filtered),
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
  file = file.path(out, "ridge_model_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  sc,
  file.path(out, "ridge_final_SIAMCAT.rds")
)

cat("\n20% Ridge LOG experiment completed.\n")
cat("\nResults across repetitions:\n")
print(repeat_results)

cat("\nModel summary:\n")
print(model_summary)

cat("\nOutputs saved in:", out, "\n")
