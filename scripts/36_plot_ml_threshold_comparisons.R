input_file <- "results/machine_learning/PRJEB27564/experiments/final_comparison/master_ML_model_comparison_32.tsv"

out <- "results/machine_learning/PRJEB27564/experiments/final_comparison"

dir.create(out, recursive = TRUE, showWarnings = FALSE)

results <- read.delim(
  input_file,
  stringsAsFactors = FALSE
)

thresholds <- c("5%", "10%", "20%", "30%")

models <- c(
  "Ridge",
  "LASSO",
  "ElasticNet",
  "RandomForest"
)

if (nrow(results) != 32) {
  stop("Expected 32 model configurations.")
}

plot_metric <- function(
    mean_col,
    sd_col,
    metric_label,
    output_file
) {

  values <- c(
    results[[mean_col]] - results[[sd_col]],
    results[[mean_col]] + results[[sd_col]]
  )

  y_range <- range(values, na.rm = TRUE)
  padding <- diff(y_range) * 0.10

  if (padding == 0) {
    padding <- 0.02
  }

  y_range <- y_range + c(-padding, padding)

  pdf(
    file.path(out, output_file),
    width = 11,
    height = 8
  )

  par(
    mfrow = c(2, 2),
    mar = c(4.5, 4.5, 3, 1)
  )

  for (model in models) {

    log_data <- results[
      results$Model == model &
      results$Normalization == "log.std",
    ]

    clr_data <- results[
      results$Model == model &
      results$Normalization == "log.clr",
    ]

    log_data <- log_data[
      match(thresholds, log_data$Prevalence),
    ]

    clr_data <- clr_data[
      match(thresholds, clr_data$Prevalence),
    ]

    x <- seq_along(thresholds)

    y_log <- log_data[[mean_col]]
    sd_log <- log_data[[sd_col]]

    y_clr <- clr_data[[mean_col]]
    sd_clr <- clr_data[[sd_col]]

    plot(
      x,
      y_log,
      type = "n",
      xaxt = "n",
      xlim = c(0.8, 4.2),
      ylim = y_range,
      xlab = "Prevalence threshold",
      ylab = metric_label,
      main = model
    )

    axis(
      1,
      at = x,
      labels = thresholds
    )

    abline(
      h = pretty(y_range),
      col = "grey85",
      lty = 3
    )

    lines(
      x,
      y_log,
      type = "b",
      col = "blue",
      lwd = 2,
      pch = 16
    )

    lines(
      x,
      y_clr,
      type = "b",
      col = "red",
      lwd = 2,
      pch = 17,
      lty = 2
    )

    arrows(
      x,
      y_log - sd_log,
      x,
      y_log + sd_log,
      angle = 90,
      code = 3,
      length = 0.04,
      col = "blue"
    )

    arrows(
      x,
      y_clr - sd_clr,
      x,
      y_clr + sd_clr,
      angle = 90,
      code = 3,
      length = 0.04,
      col = "red"
    )

    legend(
      "bottomright",
      legend = c("LOG standardisation", "CLR"),
      col = c("blue", "red"),
      lty = c(1, 2),
      pch = c(16, 17),
      lwd = 2,
      cex = 0.8,
      bty = "n"
    )
  }

  dev.off()
}

plot_metric(
  "Mean_repeat_AUROC",
  "SD_repeat_AUROC",
  "Mean repeated-CV AUROC",
  "ML_AUROC_prevalence_comparison.pdf"
)

plot_metric(
  "Mean_repeat_AUPRC",
  "SD_repeat_AUPRC",
  "Mean repeated-CV AUPRC",
  "ML_AUPRC_prevalence_comparison.pdf"
)

cat("\nBoth comparison figures created successfully.\n")

cat(
  "\nSaved in:",
  out,
  "\n"
)
