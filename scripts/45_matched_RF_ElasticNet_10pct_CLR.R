rf_file <- "results/machine_learning/PRJEB27564/experiments/10pct_CLR/RandomForest/feature_importance/RF_permutation_importance_all_genera.tsv"

enet_file <- "results/machine_learning/PRJEB27564/experiments/10pct_CLR/ElasticNet/feature_coefficients/ElasticNet_10pct_CLR_all_genera.tsv"

ancom_file <- "results/statistics/PRJEB27564/ANCOMBC2_genus/final_summary/ANCOMBC2_all_results_interpreted.tsv"

out <- "results/machine_learning/PRJEB27564/experiments/final_comparison/matched_10pct_CLR"

dir.create(out, recursive = TRUE, showWarnings = FALSE)

cat("Loading model results...\n")

rf <- read.delim(rf_file, check.names = FALSE)
enet <- read.delim(enet_file, check.names = FALSE)
ancom <- read.delim(ancom_file, check.names = FALSE)

if (!setequal(rf$Genus, enet$Genus)) {
  stop("Random Forest and ElasticNet feature sets differ.")
}

if (nrow(rf) != 158 || nrow(enet) != 158) {
  stop("Expected 158 genera in each model.")
}

rf$RF_rank <- rank(
  -rf$Mean_permutation_importance,
  ties.method = "min"
)

enet$ElasticNet_rank <- rank(
  -enet$Mean_absolute_coefficient,
  ties.method = "min"
)

rf$RF_top20 <- rf$RF_rank <= 20
enet$ElasticNet_top20 <- enet$ElasticNet_rank <= 20

cat("Merging matched model results...\n")

combined <- merge(
  rf,
  enet,
  by = "Genus",
  all = TRUE,
  sort = FALSE
)

if (nrow(combined) != 158) {
  stop("Unexpected number of merged genera.")
}

combined$Top20_overlap <- ifelse(
  combined$RF_top20 & combined$ElasticNet_top20,
  "Shared top 20",
  ifelse(
    combined$RF_top20,
    "Random Forest only",
    ifelse(
      combined$ElasticNet_top20,
      "ElasticNet only",
      "Neither top 20"
    )
  )
)

cat("\nMatching ANCOM-BC2 results...\n")

ancom$Feature_key <- make.names(ancom$taxon)
combined$Feature_key <- make.names(combined$Genus)

if (anyDuplicated(ancom$Feature_key)) {
  stop("Duplicate ANCOM-BC2 feature names.")
}

idx <- match(
  combined$Feature_key,
  ancom$Feature_key
)

combined$ANCOM_taxon <- ancom$taxon[idx]
combined$ANCOM_lfc <- ancom$lfc_parkinsonyes[idx]
combined$ANCOM_q <- ancom$q_parkinsonyes[idx]
combined$ANCOM_primary_significant <-
  ancom$diff_parkinsonyes[idx]
combined$ANCOM_passed_sensitivity <-
  ancom$passed_ss_parkinsonyes[idx]
combined$ANCOM_robust_significant <-
  ancom$diff_robust_parkinsonyes[idx]

combined$ANCOM_status <- ifelse(
  is.na(combined$ANCOM_taxon),
  "Not retained or tested",
  ifelse(
    combined$ANCOM_robust_significant %in% TRUE,
    "Robust differential abundance",
    ifelse(
      combined$ANCOM_primary_significant %in% TRUE,
      "Significant but not robust",
      "Not robustly significant"
    )
  )
)

combined <- combined[
  order(combined$RF_rank),
]

write.table(
  combined,
  file = file.path(out, "matched_10pct_CLR_all_158_genera.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

shared <- combined[
  combined$Top20_overlap == "Shared top 20",
]

write.table(
  shared,
  file = file.path(out, "matched_10pct_CLR_shared_top20.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n====================================\n")
cat("MATCHED MODEL COMPARISON\n")
cat("====================================\n")

cat("Total genera:", nrow(combined), "\n")
cat("Shared top-20 genera:", nrow(shared), "\n")

cat("\nTop shared genera:\n")

print(
  shared[, c(
    "Genus",
    "Top10_frequency_percent",
    "Mean_permutation_importance",
    "Mean_coefficient",
    "Selection_frequency_percent",
    "Model_direction",
    "ANCOM_q",
    "ANCOM_status"
  )],
  row.names = FALSE
)

cat("\nOverlap summary:\n")
print(table(combined$Top20_overlap))

cat("\nANCOM-BC2 status among shared genera:\n")
print(table(shared$ANCOM_status))

cat("\nMatched comparison completed.\n")
cat("Results saved in:", out, "\n")
