
suppressPackageStartupMessages({
  library(SIAMCAT)
})

input_path <- "results/machine_learning/PRJEB27564/input"

ridge_path <- "results/machine_learning/PRJEB27564/experiments/5pct_LOG_Ridge"

out <- "results/machine_learning/PRJEB27564/experiments/5pct_CLR"

dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Load original SIAMCAT object
base_sc <- readRDS(
  file.path(input_path, "siamcat_genus_131.rds")
)

# Apply the same 5% prevalence filter
base_sc <- filter.features(
  base_sc,
  filter.method = "prevalence",
  cutoff = 0.05,
  rm.unmapped = FALSE,
  feature.type = "original",
  verbose = 1
)

filtered <- get.filt_feat.matrix(base_sc)

# Calculate the same pseudocount convention
positive <- filtered[filtered > 0]

if (length(positive) == 0) {
  stop("No positive values available.")
}

log.n0 <- min(positive) / 100

# CLR normalization
base_sc <- normalize.features(
  base_sc,
  norm.method = "log.clr",
  norm.param = list(log.n0 = log.n0),
  feature.type = "filtered",
  verbose = 1
)

norm_features <- get.norm_feat.matrix(base_sc)

if (any(!is.finite(norm_features))) {
  stop("Non-finite values in CLR matrix.")
}

saveRDS(
  norm_features,
  file.path(out, "normalized_features_CLR.rds")
)

# Reuse Ridge's original splits
ridge_ref <- readRDS(
  file.path(ridge_path, "ridge_final_SIAMCAT.rds")
)

shared_split <- data_split(ridge_ref, verbose = 0)

if (shared_split$num.folds != 10 ||
    shared_split$num.resample != 10) {
  stop("Unexpected cross-validation design.")
}

cat("Samples:", ncol(filtered), "\n")
cat("Genera:", nrow(filtered), "\n")
cat("Normalization: log.clr\n")
cat("Cross-validation: 10 x 10\n")

run_model <- function(model_name, method, param.set = NULL) {

  cat("\nStarting:", model_name, "\n")

  sc <- base_sc
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
    stop("Expected results from 10 validation repetitions.")
  }

  model_dir <- file.path(out, model_name)
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(
    sc,
    file.path(model_dir, "model_SIAMCAT.rds")
  )

  repeat_metrics <- data.frame(
    Repeat = seq_along(auc),
    AUROC = auc,
    AUPRC = pr
  )

  write.table(
    repeat_metrics,
    file.path(model_dir, "metrics_by_repeat.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  summary <- data.frame(
    Model = model_name,
    Prevalence = "5%",
    Normalization = "log.clr",
    Samples = 131,
    Genera = nrow(filtered),
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

  saveRDS(
    sc,
    file.path(model_dir, "final_SIAMCAT.rds")
  )

  print(summary)

  summary
}

results <- list()

results[[1]] <- run_model("Ridge", "ridge")
results[[2]] <- run_model("LASSO", "lasso")
results[[3]] <- run_model(
  "ElasticNet", "enet",
  param.set = list(alpha = 0.5)
)

comparison <- do.call(rbind, results)

write.table(
  comparison,
  file = file.path(out, "CLR_algorithm_comparison.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nCLR linear-model comparison completed.\n")
print(comparison[, c(
  "Model",
  "Mean_curve_AUROC",
  "Mean_curve_AUPRC"
)])
