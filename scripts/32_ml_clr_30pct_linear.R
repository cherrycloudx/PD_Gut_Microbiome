suppressPackageStartupMessages({
  library(SIAMCAT)
})

input_path <- "results/machine_learning/PRJEB27564/input"

ridge_path <- "results/machine_learning/PRJEB27564/experiments/30pct_LOG_Ridge"

out <- "results/machine_learning/PRJEB27564/experiments/30pct_CLR"

dir.create(out, recursive = TRUE, showWarnings = FALSE)

cat("Loading SIAMCAT input...\n")

base_sc <- readRDS(
  file.path(input_path, "siamcat_genus_131.rds")
)

cat("\nApplying 30% prevalence filtering...\n")

base_sc <- filter.features(
  base_sc,
  filter.method = "prevalence",
  cutoff = 0.30,
  rm.unmapped = FALSE,
  feature.type = "original",
  verbose = 2
)

filtered <- get.filt_feat.matrix(base_sc)

cat("Samples:", ncol(filtered), "\n")
cat("Genera:", nrow(filtered), "\n")

if (nrow(filtered) != 115 || ncol(filtered) != 131) {
  stop("Unexpected filtered feature dimensions.")
}

positive <- filtered[filtered > 0]

if (length(positive) == 0) {
  stop("No positive values available.")
}

log.n0 <- min(positive) / 100

cat("\nApplying CLR normalization...\n")

base_sc <- normalize.features(
  base_sc,
  norm.method = "log.clr",
  norm.param = list(log.n0 = log.n0),
  feature.type = "filtered",
  verbose = 2
)

norm_features <- get.norm_feat.matrix(base_sc)

if (any(!is.finite(norm_features))) {
  stop("Non-finite values in CLR matrix.")
}

saveRDS(
  norm_features,
  file.path(out, "normalized_features_CLR.rds")
)

cat("\nLoading 30% Ridge reference splits...\n")

ridge_ref <- readRDS(
  file.path(ridge_path, "ridge_final_SIAMCAT.rds")
)

shared_split <- data_split(ridge_ref, verbose = 0)

if (is.null(shared_split) ||
    shared_split$num.folds != 10 ||
    shared_split$num.resample != 10) {
  stop("Unexpected cross-validation design.")
}

reference_samples <- colnames(get.filt_feat.matrix(ridge_ref))

if (!identical(colnames(filtered), reference_samples)) {
  stop("Sample order differs from the 30% Ridge reference.")
}

cat("Using the same 10-fold x 10-repetition splits.\n")

run_model <- function(model_name, method, param.set = NULL) {

  cat("\n====================================\n")
  cat("Starting:", model_name, "\n")
  cat("====================================\n")

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
    stop(paste("Expected results from 10 repetitions for", model_name))
  }

  model_dir <- file.path(out, model_name)

  dir.create(
    model_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  saveRDS(
    sc,
    file.path(model_dir, "model_SIAMCAT.rds")
  )

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
    Prevalence = "30%",
    Normalization = "log.clr",
    Samples = 131,
    Genera = 115,
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

  cat("\n", model_name, "results:\n")
  print(summary)

  return(summary)
}

results <- list()

results[[1]] <- run_model("Ridge", "ridge")

results[[2]] <- run_model("LASSO", "lasso")

results[[3]] <- run_model(
  "ElasticNet",
  "enet",
  param.set = list(alpha = 0.5)
)

comparison <- do.call(rbind, results)

comparison <- comparison[
  order(-comparison$Mean_curve_AUROC),
]

write.table(
  comparison,
  file = file.path(out, "CLR_algorithm_comparison.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n30% CLR linear-model comparison completed.\n")

print(comparison[, c(
  "Model",
  "Mean_curve_AUROC",
  "Mean_curve_AUPRC"
)])

cat("\nOutputs saved in:", out, "\n")
