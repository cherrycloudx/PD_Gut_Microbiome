ab_file <- "results/taxonomy/PRJEB27564/abundance_matrices_131/genus_relative_abundance_assigned_reads.tsv"

md_file <- "results/taxonomy/PRJEB27564/abundance_matrices_131/matched_sample_metadata.tsv"

importance_file <- "results/machine_learning/PRJEB27564/experiments/10pct_CLR/RandomForest/feature_importance/RF_top20_important_genera.tsv"

out <- "results/machine_learning/PRJEB27564/experiments/10pct_CLR/RandomForest/feature_importance/PD_vs_Control"

dir.create(out, recursive = TRUE, showWarnings = FALSE)

cat("Loading abundance data and metadata...\n")

ab <- read.delim(ab_file, check.names = FALSE)
md <- read.delim(md_file, check.names = FALSE)
importance <- read.delim(importance_file, check.names = FALSE)

if (!"SampleID" %in% names(ab)) {
  stop("SampleID column missing.")
}

if (!all(c("Run", "parkinson") %in% names(md))) {
  stop("Required metadata columns missing.")
}

if (anyDuplicated(ab$SampleID) || anyDuplicated(md$Run)) {
  stop("Duplicate sample identifiers found.")
}

if (!setequal(ab$SampleID, md$Run)) {
  stop("Abundance and metadata sample IDs do not match.")
}

md <- md[match(ab$SampleID, md$Run), ]

if (!identical(ab$SampleID, md$Run)) {
  stop("Sample alignment failed.")
}

group <- tolower(trimws(as.character(md$parkinson)))

if (!all(group %in% c("yes", "no"))) {
  stop("Unexpected disease-status values.")
}

pd <- group == "yes"
control <- group == "no"

cat("\nGroup sizes:\n")
print(table(group))

if (sum(pd) != 67 || sum(control) != 64) {
  stop("Unexpected group sizes.")
}

genus_cols <- setdiff(names(ab), "SampleID")

ab_mat <- as.matrix(ab[, genus_cols, drop = FALSE])
storage.mode(ab_mat) <- "numeric"

if (any(!is.finite(ab_mat)) || any(ab_mat < 0)) {
  stop("Invalid abundance values detected.")
}

row_totals <- rowSums(ab_mat)

cat("\nChecking relative abundance percentages...\n")
print(summary(row_totals))

if (any(abs(row_totals - 100) > 1e-3)) {
  stop("Abundance percentages do not sum to approximately 100.")
}

cat("All samples sum to approximately 100%.\n")

cat("\nTesting all 315 genera...\n")

stats <- lapply(seq_along(genus_cols), function(i) {

  x_pd <- ab_mat[pd, i]
  x_control <- ab_mat[control, i]

  test <- wilcox.test(
    x_pd,
    x_control,
    exact = FALSE
  )

  mean_pd <- mean(x_pd)
  mean_control <- mean(x_control)

  median_pd <- median(x_pd)
  median_control <- median(x_control)

  data.frame(
    Genus = genus_cols[i],
    Mean_PD_percent = mean_pd,
    Mean_Control_percent = mean_control,
    Median_PD_percent = median_pd,
    Median_Control_percent = median_control,
    Mean_difference_percentage_points =
      mean_pd - mean_control,
    Median_difference_percentage_points =
      median_pd - median_control,
    Prevalence_PD_percent = mean(x_pd > 0) * 100,
    Prevalence_Control_percent = mean(x_control > 0) * 100,
    P_value = test$p.value
  )
})

stats <- do.call(rbind, stats)

stats$BH_adjusted_p <- p.adjust(
  stats$P_value,
  method = "BH"
)

stats$Direction_by_median <- ifelse(
  stats$Median_difference_percentage_points > 0,
  "Higher in PD",
  ifelse(
    stats$Median_difference_percentage_points < 0,
    "Higher in Control",
    "Equal medians"
  )
)

stats <- stats[
  order(stats$BH_adjusted_p, stats$P_value),
]

write.table(
  stats,
  file = file.path(out, "all_315_genera_PD_vs_Control.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nMatching Random Forest top-20 genera...\n")

model_names <- make.names(genus_cols)

if (anyDuplicated(model_names)) {
  stop("Duplicate names after R feature-name conversion.")
}

importance_index <- match(
  importance$Genus,
  model_names
)

if (any(is.na(importance_index))) {
  cat("Unmatched genera:\n")
  print(importance$Genus[is.na(importance_index)])
  stop("Some model genera could not be matched.")
}

importance$Abundance_genus <- genus_cols[importance_index]

top20 <- stats[
  match(importance$Abundance_genus, stats$Genus),
]

if (any(is.na(top20$Genus))) {
  stop("Failed to match top-20 genera.")
}

top20 <- cbind(
  importance,
  top20[, setdiff(names(top20), "Genus"), drop = FALSE]
)

write.table(
  top20,
  file = file.path(out, "RF_top20_genera_PD_vs_Control.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n====================================\n")
cat("TOP 20 RANDOM FOREST GENERA\n")
cat("====================================\n")

print(top20[, c(
  "Genus",
  "Top10_frequency_percent",
  "Mean_PD_percent",
  "Mean_Control_percent",
  "Median_difference_percentage_points",
  "Direction_by_median",
  "P_value",
  "BH_adjusted_p"
)], row.names = FALSE)

cat("\n====================================\n")
cat("ALL-GENERA SIGNIFICANCE SUMMARY\n")
cat("====================================\n")

cat(
  "Genera with BH-adjusted p < 0.05:",
  sum(stats$BH_adjusted_p < 0.05),
  "\n"
)

cat("\nAnalysis completed.\n")
cat("Output directory:", out, "\n")
