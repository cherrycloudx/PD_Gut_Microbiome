suppressPackageStartupMessages({
  library(ANCOMBC)
})

result_file <- "results/statistics/PRJEB27564/ANCOMBC2_genus/ANCOMBC2_full_output.rds"

rf_file <- "results/machine_learning/PRJEB27564/experiments/10pct_CLR/RandomForest/feature_importance/RF_top20_important_genera.tsv"

out <- "results/statistics/PRJEB27564/ANCOMBC2_genus/final_summary"

dir.create(out, recursive = TRUE, showWarnings = FALSE)

cat("Loading ANCOM-BC2 results...\n")

fit <- readRDS(result_file)

res <- as.data.frame(fit$res)

lfc_col <- "lfc_parkinsonyes"
p_col <- "p_parkinsonyes"
q_col <- "q_parkinsonyes"
diff_col <- "diff_parkinsonyes"
sens_col <- "passed_ss_parkinsonyes"
robust_col <- "diff_robust_parkinsonyes"

required <- c(
  "taxon", lfc_col, p_col, q_col,
  diff_col, sens_col, robust_col
)

if (!all(required %in% names(res))) {
  stop("Required ANCOM-BC2 result columns are missing.")
}

cat("Taxa in ANCOM-BC2 results:", nrow(res), "\n")

res$Estimated_PD_to_Control_fold_change <- exp(res[[lfc_col]])

res$Direction <- ifelse(
  res[[lfc_col]] > 0,
  "Higher in PD",
  ifelse(
    res[[lfc_col]] < 0,
    "Higher in Control",
    "Equal"
  )
)

res$BH_significant <- res[[diff_col]]
res$Passed_sensitivity <- res[[sens_col]]
res$Robust_differential_abundance <- res[[robust_col]]

res <- res[
  order(res[[q_col]], na.last = TRUE),
]

write.table(
  res,
  file = file.path(out, "ANCOMBC2_all_results_interpreted.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

significant <- res[
  res[[diff_col]] %in% TRUE,
]

robust <- res[
  res[[robust_col]] %in% TRUE,
]

significant_not_robust <- res[
  res[[diff_col]] %in% TRUE &
    !(res[[robust_col]] %in% TRUE),
]

write.table(
  significant,
  file = file.path(out, "ANCOMBC2_BH_significant.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  robust,
  file = file.path(out, "ANCOMBC2_robust_significant.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  significant_not_robust,
  file = file.path(out, "ANCOMBC2_significant_not_robust.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n====================================\n")
cat("ANCOM-BC2 SUMMARY\n")
cat("====================================\n")

cat("Total tested taxa:", nrow(res), "\n")
cat("BH-significant taxa:", nrow(significant), "\n")
cat("Significant and robust taxa:", nrow(robust), "\n")
cat("Significant but not robust:", nrow(significant_not_robust), "\n")

cat("\nTOP RESULTS BY ADJUSTED P-VALUE:\n")

print(
  head(
    res[, c(
      "taxon",
      lfc_col,
      q_col,
      diff_col,
      sens_col,
      robust_col,
      "Direction",
      "Estimated_PD_to_Control_fold_change"
    )],
    20
  ),
  row.names = FALSE
)

cat("\n====================================\n")
cat("COMPARISON WITH RANDOM FOREST TOP 20\n")
cat("====================================\n")

rf <- read.delim(
  rf_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rf$Feature_key <- make.names(rf$Genus)
res$Feature_key <- make.names(res$taxon)

if (anyDuplicated(res$Feature_key)) {
  stop("Duplicate feature keys in ANCOM-BC2 results.")
}

lookup <- res[, c(
  "Feature_key",
  "taxon",
  lfc_col,
  q_col,
  diff_col,
  sens_col,
  robust_col,
  "Direction",
  "Estimated_PD_to_Control_fold_change"
)]

idx <- match(rf$Feature_key, lookup$Feature_key)

comparison <- cbind(
  rf,
  lookup[idx, setdiff(names(lookup), "Feature_key"), drop = FALSE]
)

comparison$ANCOMBC2_status <- ifelse(
  is.na(comparison$taxon),
  "Not retained or not tested",
  ifelse(
    comparison[[robust_col]] %in% TRUE,
    "Robust differential abundance",
    ifelse(
      comparison[[diff_col]] %in% TRUE,
      "Significant but not robust",
      "Not significant after BH and sensitivity analysis"
    )
  )
)

write.table(
  comparison,
  file = file.path(out, "RF_top20_vs_ANCOMBC2.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(
  comparison[, c(
    "Genus",
    "Top10_frequency_percent",
    lfc_col,
    q_col,
    robust_col,
    "ANCOMBC2_status"
  )],
  row.names = FALSE
)

cat("\n====================================\n")
cat("ANCOM-BC2 SUMMARY COMPLETED\n")
cat("====================================\n")

cat("Output directory:", out, "\n")
