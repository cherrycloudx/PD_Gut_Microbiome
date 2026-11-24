
suppressPackageStartupMessages({
  library(SIAMCAT)
})

set.seed(56987)

input_path <- "results/machine_learning/PRJEB27564/input"

ridge_path <- "results/machine_learning/PRJEB27564/experiments/5pct_LOG_Ridge"

comparison_path <- "results/machine_learning/PRJEB27564/experiments"

out <- file.path(comparison_path, "5pct_CLR/RandomForest")

dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Load original input
sc <- readRDS(
  file.path(input_path, "siamcat_genus_131.rds")
)

# Apply 5% prevalence filtering
sc <- filter.features(
  sc,
  filter.method = "prevalence",
  cutoff = 0.05,
  rm.unmapped = FALSE,
  feature.type = "original",
  verbose = 1
)

filtered <- get.filt_feat.matrix(sc)

# Determine pseudocount
positive <- filtered[filtered > 0]
log.n0 <- min(positive) / 100

# CLR normalization
sc <- normalize.features(
  sc,
  norm.method = "log.clr",
  norm.param = list(log.n0 = log.n0),
  feature.type = "filtered",
  verbose = 1
)

# Reuse Ridge cross-validation splits
ridge_ref <- readRDS(
  file.path(ridge_path, "ridge_final_SIAMCAT.rds")
)

shared_split <- data_split(ridge_ref, verbose = 0)

if (shared_split$num.folds != 10 ||
    shared_split$num.resample != 10) {
  stop("Unexpected cross-validation structure.")
}

data_split(sc) <- shared_split

if (nrow(filtered) != 198 || ncol(filtered) != 131) {
  stop("Unexpected feature matrix dimensions.")
}

cat("Samples:", ncol(filtered), "\n")
cat("Genera:", nrow(filtered), "\n")
cat("Normalization: log.clr\n")
cat("Cross-validation: 10 folds x 10 repetitions\n")

# Train Random Forest
cat("\nTraining Random Forest...\n")

sc <- train.model(
  sc,
  method = "randomForest",
  measure = "classif.acc",
  feature.type = "normalized",
  verbose = 2
)

saveRDS(sc, file.path(out, "randomforest_trained.rds"))

# Predictions and evaluation
cat("\nGenerating predictions...\n")

sc <- make.predictions(sc, verbose = 2)

cat("\nEvaluating model...\n")

sc <- evaluate.predictions(sc, verbose = 2)

evaluation <- eval_data(sc)

auc <- as.numeric(evaluation$auroc.all)
auprc <- as.numeric(evaluation$auprc.all)

if (length(auc) != 10 || length(auprc) != 10) {
  stop("Expected 10 repetition results.")
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

summary_result <- data.frame(
  Model = "RandomForest",
  Prevalence = "5%",
  Normalization = "log.clr",
  Samples = 131,
  Genera = 198,
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
  summary_result,
  file = file.path(out, "model_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  sc,
  file.path(out, "randomforest_final_SIAMCAT.rds")
)

# Complete CLR comparison
clr_path <- file.path(comparison_path, "5pct_CLR")

clr_comparison <- read.delim(
  file.path(clr_path, "CLR_algorithm_comparison.tsv")
)

clr_comparison <- rbind(clr_comparison, summary_result)

write.table(
  clr_comparison,
  file = file.path(clr_path, "CLR_algorithm_comparison.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Combine all eight model configurations
log_comparison <- read.delim(
  file.path(
    comparison_path,
    "5pct_LOG_comparison/algorithm_comparison.tsv"
  )
)

log_comparison$Normalization <- "log.std"

all_results <- rbind(log_comparison, clr_comparison)

all_results <- all_results[
  order(-all_results$Mean_curve_AUROC),
]

write.table(
  all_results,
  file = file.path(
    comparison_path,
    "5pct_LOG_vs_CLR_all_algorithms.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nRandom Forest CLR results:\n")
print(summary_result)

cat("\nAll eight model configurations:\n")
print(all_results[, c(
  "Model",
  "Normalization",
  "Mean_curve_AUROC",
  "Mean_curve_AUPRC"
)])

cat("\nOutputs saved in:", out, "\n")
