suppressPackageStartupMessages({
  library(SIAMCAT)
})

model_file <- "results/machine_learning/PRJEB27564/experiments/10pct_CLR/ElasticNet/final_SIAMCAT.rds"

out <- "results/machine_learning/PRJEB27564/experiments/10pct_CLR/ElasticNet/feature_coefficients"

dir.create(out, recursive = TRUE, showWarnings = FALSE)

cat("Loading 10% CLR ElasticNet model...\n")

sc <- readRDS(model_file)

if (model_type(sc) != "enet") {
  stop("Expected ElasticNet model.")
}

label_info <- methods::slot(sc, "label")$info

if (unname(label_info["no"]) != -1 ||
    unname(label_info["yes"]) != 1) {
  stop("Unexpected disease-status encoding.")
}

fit_list <- models(sc, verbose = 0)

if (length(fit_list) != 100) {
  stop("Expected 100 trained models.")
}

cat("Models:", length(fit_list), "\n")

beta_list <- lapply(fit_list, function(fit) {

  glmnet_model <- fit$model$model$glmnet.fit

  call_alpha <- fit$model$model$call$alpha

  if (is.null(call_alpha) || call_alpha != 0.5) {
    stop("A model does not have the expected alpha = 0.5.")
  }

  class_names <- glmnet_model$classnames

  if (!identical(as.character(class_names), c("-1", "1"))) {
    stop("Unexpected GLMNET class ordering.")
  }

  stored <- fit$features

  if (is.null(stored) ||
      is.null(names(stored)) ||
      !is.numeric(stored)) {
    stop("Invalid feature scores.")
  }

  beta <- -as.numeric(stored)
  names(beta) <- names(stored)

  beta
})

feature_names <- names(beta_list[[1]])

if (length(feature_names) != 158) {
  stop("Expected 158 genera.")
}

if (!all(vapply(
  beta_list,
  function(x) setequal(names(x), feature_names),
  logical(1)
))) {
  stop("Feature names differ between models.")
}

beta_matrix <- do.call(
  rbind,
  lapply(beta_list, function(x) x[feature_names])
)

rownames(beta_matrix) <- names(fit_list)

if (any(!is.finite(beta_matrix))) {
  stop("Non-finite coefficients detected.")
}

selected <- abs(beta_matrix) > 1e-10

positive <- beta_matrix > 1e-10
negative <- beta_matrix < -1e-10

selected_count <- colSums(selected)
positive_count <- colSums(positive)
negative_count <- colSums(negative)

rank_matrix <- t(
  apply(
    abs(beta_matrix),
    1,
    function(x) rank(-x, ties.method = "average")
  )
)

total_abs <- rowSums(abs(beta_matrix))

relative_beta <- beta_matrix

for (i in seq_len(nrow(beta_matrix))) {
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
  Selection_frequency_percent = selected_count,
  Positive_PD_count = positive_count,
  Negative_Control_count = negative_count,
  Positive_PD_frequency_percent = positive_count,
  Negative_Control_frequency_percent = negative_count,
  Sign_consistency_percent = ifelse(
    selected_count > 0,
    pmax(positive_count, negative_count) /
      selected_count * 100,
    NA_real_
  ),
  Mean_relative_signed_coefficient = colMeans(relative_beta),
  Mean_absolute_coefficient_rank = colMeans(rank_matrix)
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
  file = file.path(out, "ElasticNet_10pct_CLR_all_genera.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  beta_matrix,
  file = file.path(out, "ElasticNet_10pct_CLR_coefficients_by_model.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

top20 <- head(summary, 20)

write.table(
  top20,
  file = file.path(out, "ElasticNet_10pct_CLR_top20.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n====================================\n")
cat("TOP 20: 10% CLR ELASTICNET\n")
cat("====================================\n")

print(
  top20[, c(
    "Genus",
    "Mean_coefficient",
    "Selection_frequency_percent",
    "Positive_PD_frequency_percent",
    "Negative_Control_frequency_percent",
    "Sign_consistency_percent",
    "Model_direction"
  )],
  row.names = FALSE
)

cat("\nModels analysed:", nrow(beta_matrix), "\n")
cat("Genera analysed:", ncol(beta_matrix), "\n")
cat("Results saved in:", out, "\n")
