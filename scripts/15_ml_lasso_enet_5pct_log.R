
suppressPackageStartupMessages({
  library(SIAMCAT)
})

input_path <- "results/machine_learning/PRJEB27564/input"

ridge_path <- "results/machine_learning/PRJEB27564/experiments/5pct_LOG_Ridge"

out <- "results/machine_learning/PRJEB27564/experiments/5pct_LOG_comparison"

dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Load preprocessed data and Ridge reference
base_sc <- readRDS(file.path(input_path, "siamcat_genus_131.rds"))

# Recreate filtering and normalisation exactly as in the Ridge experiment
base_sc <- filter.features(
  base_sc,
  filter.method = "prevalence",
  cutoff = 0.05,
  rm.unmapped = FALSE,
  feature.type = "original",
  verbose = 1
)

filtered <- get.filt_feat.matrix(base_sc)
positive <- filtered[filtered > 0]
log.n0 <- min(positive) / 100

base_sc <- normalize.features(
  base_sc,
  norm.method = "log.std",
  norm.param = list(log.n0 = log.n0, sd.min.q = 0.1),
  feature.type = "filtered",
  verbose = 1
)

ridge_ref <- readRDS(
  file.path(ridge_path, "ridge_final_SIAMCAT.rds")
)

shared_split <- data_split(ridge_ref, verbose = 0)

if (is.null(shared_split) ||
    shared_split$num.folds != 10 ||
    shared_split$num.resample != 10) {
  stop("Ridge reference does not contain the expected 10 x 10 splits.")
}

if (nrow(filtered) != 198 || ncol(filtered) != 131) {
  stop("Unexpected filtered feature dimensions.")
}

cat("Using Ridge's existing 10-fold x 10-repeat splits.\n")
cat("Features:", nrow(filtered), "\n")
cat("Samples:", ncol(filtered), "\n")

run_model <- function(model_name, method, param.set = NULL) {

  cat("\n====================================\n")
  cat("Starting:", model_name, "\n")
  cat("====================================\n")

  sc <- base_sc

  # Reuse the exact Ridge train/test splits
  data_split(sc) <- shared_split

  set.seed(56987)

  sc <- train.model(
    sc,
    method = method,
    measure = "classif.acc",
    param.set = param.set,
    feature.type = "normalized",
    verbose = 2
  )

  sc <- make.predictions(sc, verbose = 2)
  sc <- evaluate.predictions(sc, verbose = 2)

  evaluation <- eval_data(sc)

  auc <- as.numeric(evaluation$auroc.all)
  pr <- as.numeric(evaluation$auprc.all)

  if (length(auc) != 10 || length(pr) != 10) {
    stop(paste("Unexpected repetition results for", model_name))
  }

  model_dir <- file.path(out, model_name)
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(sc, file.path(model_dir, "model_SIAMCAT.rds"))

  metrics <- data.frame(
    Repeat = seq_along(auc),
    AUROC = auc,
    AUPRC = pr
  )

  write.table(
    metrics,
    file.path(model_dir, "metrics_by_repeat.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  summary <- data.frame(
    Model = model_name,
    Prevalence = "5%",
    Normalization = "log.std",
    Samples = 131,
    Genera = 198,
    Folds = 10,
    Repetitions = 10,
    Mean_repeat_AUROC = mean(auc),
    SD_repeat_AUROC = sd(auc),
    Mean_curve_AUROC = as.numeric(evaluation$auroc),
    Mean_repeat_AUPRC = mean(pr),
    SD_repeat_AUPRC = sd(pr),
    Mean_curve_AUPRC = as.numeric(evaluation$auprc)
  )

  write.table(
    summary,
    file.path(model_dir, "model_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("\n", model_name, "results:\n")
  print(summary)

  summary
}

lasso_result <- run_model(
  "LASSO",
  "lasso"
)

enet_result <- run_model(
  "ElasticNet",
  "enet",
  param.set = list(alpha = 0.5)
)

ridge_result <- read.delim(
  file.path(ridge_path, "ridge_model_summary.tsv")
)

comparison <- rbind(
  ridge_result,
  lasso_result,
  enet_result
)

write.table(
  comparison,
  file.path(out, "algorithm_comparison.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nAlgorithm comparison completed:\n")
print(comparison[, c(
  "Model",
  "Mean_curve_AUROC",
  "Mean_curve_AUPRC"
)])
