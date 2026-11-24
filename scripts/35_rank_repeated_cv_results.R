input <- "results/machine_learning/PRJEB27564/experiments/final_comparison/master_ML_model_comparison_32.tsv"

out <- "results/machine_learning/PRJEB27564/experiments/final_comparison"

x <- read.delim(input, check.names = FALSE)

metrics <- c(
  "Mean_repeat_AUROC",
  "Mean_repeat_AUPRC"
)

for (metric in metrics) {

  ranked <- x[order(-x[[metric]]), ]

  output_file <- file.path(
    out,
    paste0("ranked_by_", metric, ".tsv")
  )

  write.table(
    ranked,
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("\n====================================\n")
  cat("RANKING:", metric, "\n")
  cat("====================================\n")

  print(head(ranked[, c(
    "Model",
    "Prevalence",
    "Normalization",
    metric,
    if (metric == "Mean_repeat_AUROC") "SD_repeat_AUROC" else "SD_repeat_AUPRC"
  )], 10), row.names = FALSE)
}

cat("\nRepeated-CV rankings saved successfully.\n")
