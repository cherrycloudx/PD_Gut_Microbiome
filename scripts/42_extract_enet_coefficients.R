suppressPackageStartupMessages({
  library(SIAMCAT)
})

model_file <- "results/machine_learning/PRJEB27564/experiments/30pct_LOG_comparison/ElasticNet/model_SIAMCAT.rds"

out <- "results/machine_learning/PRJEB27564/experiments/30pct_LOG_comparison/ElasticNet/feature_coefficients"

dir.create(out, recursive = TRUE, showWarnings = FALSE)

cat("Loading ElasticNet model...\n")

sc <- readRDS(model_file)

if (model_type(sc) != "enet") {
  stop("Expected an ElasticNet model.")
}

label_info <- methods::slot(sc, "label")$info

if (unname(label_info["no"]) != -1 ||
    unname(label_info["yes"]) != 1) {
  stop("Unexpected original label encoding.")
}

cat("Confirmed: Control = -1; PD = +1.\n")

fit_list <- models(sc, verbose = 0)

if (length(fit_list) != 100) {
  stop("Expected 100 trained models.")
}

cat("Trained models:", length(fit_list), "\n")

beta_list <- lapply(fit_list, function(fit) {

  stored <- fit$features

  if (is.null(stored) ||
      !is.numeric(stored) ||
      is.null(names(stored))) {
    stop("Invalid stored feature weights.")
  }

  class_names <- fit$model$model$glmnet.fit$classnames

  if (!identical(as.character(class_names), c("-1", "1"))) {
    stop("Unexpected GLMNET class ordering.")
  }

  beta <- -as.numeric(stored)
  names(beta) <- names(stored)

  beta
})

feature_names <- names(beta_list[[1]])

if (length(feature_names) != 115) {
  stop("Expected 115 genera.")
}

consistent <- vapply(
  beta_list,
  function(x) setequal(names(x), feature_names),
  logical(1)
)

if (!all(consistent)) {
  stop("Feature names differ across models.")
}

beta_matrix <- do.call(
  rbind,
  lapply(beta_list, function(x) x[feature_names])
)

rownames(beta_matrix) <- names(fit_list)

if (any(!is.finite(beta_matrix))) {
  stop("Non-finite coefficients detected.")
}

cat("\nCoefficient matrix dimensions:\n")
print(dim(beta_matrix))

cat("\nCalculating feature stability...\n")

selected <- abs(beta_matrix) > 1e-10

positive <- beta_matrix > 1e-10
negative <- beta_matrix < -1e-10

selected_count <- colSums(selected)

positive_count <- colSums(positive)
negative_count <- colSums(negative)

sign_consistency <- ifelse(
  selected_count > 0,
  pmax(positive_count, negative_count) / selected_count * 100,
  NA_real_
)

rank_matrix <- t(
  apply(
    abs(beta_matrix),
    1,
    function(x) rank(-x, ties.method = "average")
  )
)

total_abs <- rowSums(abs(beta_matrix))

relative_beta <- beta_matrix

for (i in seq_len(nrow(relative_beta))) {
  if (total_abs[i] > 0) {
    relative_beta[i, ] <- beta_matrix[i, ] / total_abs[i]
  } else {
    relative_beta[i, ] <- 0
  }
}

summary <- data.frame(
  Genus = feature_names,
  Mean_coefficient = colMeans(beta_matrix),
  Median_coefficient = apply(beta_matrix, 2, median),
  SD_coefficient = apply(beta_matrix, 2, sd),
  Mean_absolute_coefficient = colMeans(abs(beta_matrix)),
  Selection_count = selected_count,
  Selection_frequency_percent = selected_count / 100 * 100,
  Positive_PD_coefficient_count = positive_count,
  Negative_Control_coefficient_count = negative_count,
  Positive_PD_frequency_percent = positive_count,
  Negative_Control_frequency_percent = negative_count,
  Sign_consistency_among_selected_percent = sign_consistency,
  Mean_relative_signed_coefficient = colMeans(relative_beta),
  Mean_absolute_coefficient_rank = colMeans(rank_matrix),
  stringsAsFactors = FALSE
)

summary$Model_direction <- ifelse(
  positive_count > negative_count,
  "PD-associated model contribution",
  ifelse(
    negative_count > positive_count,
    "Control-associated model contribution",
    "Mixed or unselected"
  )
)

summary <- summary[
  order(
    -summary$Selection_frequency_percent,
    -summary$Mean_absolute_coefficient
  ),
]

write.table(
  summary,
  file = file.path(out, "ElasticNet_coefficient_stability_all.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  beta_matrix,
  file = file.path(out, "ElasticNet_signed_coefficients_by_model.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

top20 <- head(summary, 20)

write.table(
  top20,
  file = file.path(out, "ElasticNet_top20_stable_genera.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n====================================\n")
cat("TOP 20 ELASTICNET GENERA\n")
cat("====================================\n")

print(
  top20[, c(
    "Genus",
    "Mean_coefficient",
    "Selection_frequency_percent",
    "Positive_PD_frequency_percent",
    "Negative_Control_frequency_percent",
    "Sign_consistency_among_selected_percent",
    "Model_direction"
  )],
  row.names = FALSE
)

cat("\n====================================\n")
cat("EXTRACTION SUMMARY\n")
cat("====================================\n")

cat("Models analysed:", nrow(beta_matrix), "\n")
cat("Genera analysed:", ncol(beta_matrix), "\n")
cat("Output directory:", out, "\n")
