rf_file <- "results/machine_learning/PRJEB27564/experiments/10pct_CLR/RandomForest/feature_importance/RF_top20_important_genera.tsv"

enet_file <- "results/machine_learning/PRJEB27564/experiments/30pct_LOG_comparison/ElasticNet/feature_coefficients/ElasticNet_top20_stable_genera.tsv"

ancom_file <- "results/statistics/PRJEB27564/ANCOMBC2_genus/final_summary/ANCOMBC2_all_results_interpreted.tsv"

out <- "results/machine_learning/PRJEB27564/experiments/final_comparison/cross_algorithm_consensus"

dir.create(out, recursive = TRUE, showWarnings = FALSE)

rf <- read.delim(rf_file, check.names = FALSE)
enet <- read.delim(enet_file, check.names = FALSE)
ancom <- read.delim(ancom_file, check.names = FALSE)

rf2 <- data.frame(
  Feature_key = rf$Genus,
  Genus = rf$Genus,
  RF_top10_frequency_percent = rf$Top10_frequency_percent,
  RF_mean_permutation_importance = rf$Mean_permutation_importance
)

enet2 <- data.frame(
  Feature_key = enet$Genus,
  ElasticNet_mean_coefficient = enet$Mean_coefficient,
  ElasticNet_selection_frequency_percent = enet$Selection_frequency_percent,
  ElasticNet_direction = enet$Model_direction
)

ancom$Feature_key <- make.names(ancom$taxon)

ancom2 <- ancom[, c(
  "Feature_key",
  "lfc_parkinsonyes",
  "q_parkinsonyes",
  "diff_parkinsonyes",
  "passed_ss_parkinsonyes",
  "diff_robust_parkinsonyes",
  "Direction"
)]

if (anyDuplicated(rf2$Feature_key) ||
    anyDuplicated(enet2$Feature_key) ||
    anyDuplicated(ancom2$Feature_key)) {
  stop("Duplicate feature keys detected.")
}

common <- merge(
  rf2,
  enet2,
  by = "Feature_key"
)

cat("Shared top-20 genera:", nrow(common), "\n")

common <- merge(
  common,
  ancom2,
  by = "Feature_key",
  all.x = TRUE
)

common$Direction_alignment <- ifelse(
  is.na(common$lfc_parkinsonyes),
  "Not tested or not retained",
  ifelse(
    common$ElasticNet_mean_coefficient > 0 &
      common$lfc_parkinsonyes > 0,
    "Aligned: PD direction",
    ifelse(
      common$ElasticNet_mean_coefficient < 0 &
        common$lfc_parkinsonyes < 0,
      "Aligned: control direction",
      "Discordant"
    )
  )
)

common$ANCOMBC2_status <- ifelse(
  is.na(common$diff_robust_parkinsonyes),
  "Not tested",
  ifelse(
    common$diff_robust_parkinsonyes,
    "Robust differential abundance",
    ifelse(
      common$diff_parkinsonyes,
      "Significant but not robust",
      "Not robustly significant"
    )
  )
)

common <- common[
  order(-common$RF_top10_frequency_percent,
        -common$ElasticNet_selection_frequency_percent),
]

write.table(
  common,
  file = file.path(out, "RF_ElasticNet_ANCOMBC2_consensus.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nCROSS-ALGORITHM CONSENSUS RESULTS\n")

print(
  common[, c(
    "Genus",
    "RF_top10_frequency_percent",
    "ElasticNet_mean_coefficient",
    "ElasticNet_selection_frequency_percent",
    "ElasticNet_direction",
    "lfc_parkinsonyes",
    "q_parkinsonyes",
    "ANCOMBC2_status",
    "Direction_alignment"
  )],
  row.names = FALSE
)

cat("\nDirection agreement summary:\n")
print(table(common$Direction_alignment))

cat("\nANCOM-BC2 status summary:\n")
print(table(common$ANCOMBC2_status))

cat("\nSaved to:", out, "\n")
