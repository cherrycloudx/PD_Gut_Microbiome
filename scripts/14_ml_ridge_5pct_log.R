
suppressPackageStartupMessages({
  library(SIAMCAT)
})

set.seed(56987)

input_file <- "results/machine_learning/PRJEB27564/input/siamcat_genus_131.rds"

out <- "results/machine_learning/PRJEB27564/experiments/5pct_LOG_Ridge"

dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Check essential modelling packages
required <- c("mlr3", "mlr3learners", "glmnet")

missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing, collapse = ", ")
  )
}

cat("Loading SIAMCAT input...\n")

sc <- readRDS(input_file)

# Step 1: 5% prevalence filtering
cat("\nApplying 5% prevalence filtering...\n")

sc <- filter.features(
  sc,
  filter.method = "prevalence",
  cutoff = 0.05,
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

if (nrow(filtered) != 198) {
  stop("Filtered genus count differs from expected 198. Inspect before modelling.")
}

write.table(
  filtered,
  file = file.path(out, "filtered_genus_features.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# Step 2: Calculate pseudocount from filtered data
positive_values <- filtered[filtered > 0]

if (length(positive_values) == 0) {
  stop("No positive feature abundances.")
}

log.n0 <- min(positive_values) / 100

cat("Log pseudocount:", format(log.n0, scientific = TRUE), "\n")

# Step 3: Log-standardisation
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

# Save preprocessing object
saveRDS(
  sc,
  file.path(out, "preprocessed_SIAMCAT.rds")
)

# Step 4: Repeated stratified 10-fold cross-validation
cat("\nPreparing repeated 10-fold cross-validation...\n")

sc <- create.data.split(
  sc,
  num.folds = 10,
  num.resample = 10,
  stratify = TRUE,
  verbose = 2
)

# Step 5: Train Ridge regression
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

# Step 6: Generate predictions
cat("\nGenerating cross-validation predictions...\n")

sc <- make.predictions(sc, verbose = 2)

# Step 7: Evaluate predictions
cat("\nEvaluating model performance...\n")

sc <- evaluate.predictions(sc, verbose = 2)

evaluation <- eval_data(sc)

if (is.null(evaluation$auroc.all)) {
  stop("Repeated cross-validation AUC results are missing.")
}

auc <- as.numeric(evaluation$auroc.all)
auprc <- as.numeric(evaluation$auprc.all)

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
  Prevalence = "5%",
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

cat("\nRidge experiment completed.\n")

cat("\nResults across repetitions:\n")
print(repeat_results)

cat("\nModel summary:\n")
print(model_summary)

cat("\nOutputs saved in:", out, "\n")
