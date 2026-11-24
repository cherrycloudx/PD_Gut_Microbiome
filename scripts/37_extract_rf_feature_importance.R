suppressPackageStartupMessages({
  library(SIAMCAT)
})

model_file <- "results/machine_learning/PRJEB27564/experiments/10pct_CLR/RandomForest/randomforest_final_SIAMCAT.rds"

out <- "results/machine_learning/PRJEB27564/experiments/10pct_CLR/RandomForest/feature_importance"

dir.create(out, recursive = TRUE, showWarnings = FALSE)

cat("Loading saved Random Forest model...\n")

sc <- readRDS(model_file)

cat("\nChecking model metadata...\n")

type <- model_type(sc, verbose = 0)
feat_type <- feature_type(sc, verbose = 0)

cat("Model type:", type, "\n")
cat("Feature type:", feat_type, "\n")

if (is.null(type) || type != "randomForest") {
  stop("The saved model is not a Random Forest model.")
}

cat("\nExtracting trained models...\n")

fit_list <- models(sc, verbose = 0)

if (is.null(fit_list) || length(fit_list) != 100) {
  stop("Expected 100 trained cross-validation models.")
}

cat("Number of trained models:", length(fit_list), "\n")

importance_list <- lapply(fit_list, function(fit) {

  importance <- fit$features

  if (is.null(importance) ||
      !is.numeric(importance) ||
      is.null(names(importance))) {
    stop("Invalid or missing feature importance scores.")
  }

  if (any(!is.finite(importance))) {
    stop("Non-finite importance scores detected.")
  }

  importance
})

feature_names <- names(importance_list[[1]])

if (length(feature_names) != 158) {
  stop("Expected 158 genera in the importance scores.")
}

consistent <- vapply(
  importance_list,
  function(x) setequal(names(x), feature_names),
  logical(1)
)

if (!all(consistent)) {
  stop("Feature names differ between trained models.")
}

importance_matrix <- do.call(
  rbind,
  lapply(importance_list, function(x) {
    x[feature_names]
  })
)

rownames(importance_matrix) <- names(fit_list)
colnames(importance_matrix) <- feature_names

cat("\nImportance matrix dimensions:\n")
print(dim(importance_matrix))

cat("\nCalculating importance ranks...\n")

rank_matrix <- t(
  apply(
    importance_matrix,
    1,
    function(x) rank(-x, ties.method = "average")
  )
)

summary <- data.frame(
  Genus = colnames(importance_matrix),
  Mean_permutation_importance = colMeans(importance_matrix),
  Median_permutation_importance = apply(
    importance_matrix, 2, median
  ),
  SD_permutation_importance = apply(
    importance_matrix, 2, sd
  ),
  Mean_importance_rank = colMeans(rank_matrix),
  Top10_count = colSums(rank_matrix <= 10),
  Top10_frequency_percent = colMeans(rank_matrix <= 10) * 100,
  Positive_importance_frequency_percent =
    colMeans(importance_matrix > 0) * 100,
  stringsAsFactors = FALSE
)

summary <- summary[
  order(
    -summary$Top10_frequency_percent,
    summary$Mean_importance_rank,
    -summary$Mean_permutation_importance
  ),
]

write.table(
  summary,
  file = file.path(out, "RF_permutation_importance_all_genera.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  importance_matrix,
  file = file.path(out, "RF_permutation_importance_by_model.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  rank_matrix,
  file = file.path(out, "RF_importance_ranks_by_model.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

top20 <- head(summary, 20)

write.table(
  top20,
  file = file.path(out, "RF_top20_important_genera.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pdf(
  file.path(out, "RF_top20_importance_stability.pdf"),
  width = 10,
  height = 8
)

par(mar = c(5, 14, 4, 2))

barplot(
  rev(top20$Top10_count),
  names.arg = rev(top20$Genus),
  horiz = TRUE,
  las = 1,
  xlim = c(0, 100),
  xlab = "Number of CV models ranking genus in top 10",
  main = "Random Forest: Top 20 genus importance stability"
)

dev.off()

cat("\n====================================\n")
cat("TOP 20 IMPORTANT GENERA\n")
cat("====================================\n")

print(top20, row.names = FALSE)

cat("\nFeature importance extraction completed.\n")
cat("Models analysed:", nrow(importance_matrix), "\n")
cat("Genera analysed:", ncol(importance_matrix), "\n")
cat("Output directory:", out, "\n")
