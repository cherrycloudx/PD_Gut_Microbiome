suppressPackageStartupMessages(library(vegan))

ab_file <- "results/taxonomy/PRJEB27564/abundance_matrices_131/genus_relative_abundance_assigned_reads.tsv"
md_file <- "results/taxonomy/PRJEB27564/abundance_matrices_131/matched_sample_metadata.tsv"

out_dir <- "results/statistics/PRJEB27564/beta_diversity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading abundance matrix and metadata...\n")

ab <- read.delim(ab_file, check.names = FALSE, stringsAsFactors = FALSE)
md <- read.delim(md_file, check.names = FALSE, stringsAsFactors = FALSE)

sample_ids <- as.character(ab$SampleID)

if (anyDuplicated(sample_ids)) {
  stop("Duplicate sample IDs in abundance matrix.")
}

if (anyDuplicated(as.character(md$Run))) {
  stop("Duplicate Run IDs in metadata.")
}

if (!setequal(sample_ids, as.character(md$Run))) {
  stop("Sample IDs do not match between abundance matrix and metadata.")
}

md <- md[match(sample_ids, as.character(md$Run)), , drop = FALSE]

if (!identical(sample_ids, as.character(md$Run))) {
  stop("Metadata sample order could not be matched.")
}

if (!all(as.character(md$parkinson) %in% c("no", "yes"))) {
  stop("Unexpected or missing Parkinson group labels.")
}

md$parkinson <- factor(
  as.character(md$parkinson),
  levels = c("no", "yes"),
  labels = c("Control", "PD")
)

taxa <- ab[, setdiff(names(ab), "SampleID"), drop = FALSE]

taxa[] <- lapply(taxa, function(x) as.numeric(as.character(x)))

if (anyNA(taxa) || any(!is.finite(as.matrix(taxa)))) {
  stop("Missing or non-finite abundance values detected.")
}

if (any(as.matrix(taxa) < 0)) {
  stop("Negative abundance values detected.")
}

rownames(taxa) <- sample_ids

if (any(abs(rowSums(taxa) - 100) > 0.001)) {
  stop("Some sample abundance totals are not approximately 100%.")
}

keep <- colSums(taxa) > 0
taxa <- taxa[, keep, drop = FALSE]

cat("Samples:", nrow(taxa), "\n")
cat("Genus features retained:", ncol(taxa), "\n")
print(table(md$parkinson))

cat("\nCalculating Bray-Curtis dissimilarity...\n")

bray <- vegan::vegdist(taxa, method = "bray")

dist_matrix <- as.matrix(bray)

write.table(
  data.frame(SampleID = rownames(dist_matrix),
             dist_matrix, check.names = FALSE),
  file = file.path(out_dir, "Bray_Curtis_distance_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Bray-Curtis distance matrix saved.\n")

cat("\nPerforming PCoA...\n")

pcoa <- cmdscale(bray, k = 2, eig = TRUE, add = TRUE)

scores <- as.data.frame(pcoa$points)
colnames(scores) <- c("PCoA1", "PCoA2")

eig <- pcoa$eig
positive_eig <- eig[eig > 0]

axis1_pct <- eig[1] / sum(positive_eig) * 100
axis2_pct <- eig[2] / sum(positive_eig) * 100

coords <- data.frame(
  SampleID = rownames(scores),
  PCoA1 = scores$PCoA1,
  PCoA2 = scores$PCoA2,
  parkinson = md$parkinson
)

write.table(
  coords,
  file = file.path(out_dir, "PCoA_coordinates.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pdf(
  file.path(out_dir, "Bray_Curtis_PCoA.pdf"),
  width = 8,
  height = 6
)

group_colours <- ifelse(
  md$parkinson == "Control",
  "#2878B5",
  "#D95F02"
)

plot(
  scores$PCoA1,
  scores$PCoA2,
  type = "n",
  xlab = paste0("PCoA1 (", round(axis1_pct, 2), "%)"),
  ylab = paste0("PCoA2 (", round(axis2_pct, 2), "%)"),
  main = "Bray-Curtis PCoA at Genus Level"
)

points(
  scores$PCoA1,
  scores$PCoA2,
  col = group_colours,
  pch = 19,
  cex = 1.1
)

legend(
  "topright",
  legend = c("Control", "PD"),
  col = c("#2878B5", "#D95F02"),
  pch = 19,
  bty = "n"
)

dev.off()

cat("PCoA plot saved.\n")

cat("\nRunning PERMANOVA...\n")

set.seed(20261122)

permanova <- vegan::adonis2(
  bray ~ parkinson,
  data = md,
  permutations = 9999
)

permanova_df <- as.data.frame(permanova)
permanova_df$Term <- rownames(permanova_df)

write.table(
  permanova_df,
  file = file.path(out_dir, "PERMANOVA_Bray_Curtis.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(permanova)

cat("\nRunning PERMDISP...\n")

dispersion <- vegan::betadisper(
  bray,
  group = md$parkinson,
  type = "median",
  bias.adjust = TRUE
)

set.seed(20261122)

dispersion_test <- vegan::permutest(
  dispersion,
  permutations = 9999
)

dispersion_df <- as.data.frame(dispersion_test$tab)
dispersion_df$Term <- rownames(dispersion_df)

write.table(
  dispersion_df,
  file = file.path(out_dir, "PERMDISP_Bray_Curtis.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(dispersion_test)

dispersion_samples <- data.frame(
  SampleID = sample_ids,
  Group = md$parkinson,
  Distance_to_group_median = dispersion$distances
)

write.table(
  dispersion_samples,
  file = file.path(out_dir, "PERMDISP_sample_distances.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_file <- file.path(out_dir, "beta_diversity_summary.txt")

sink(summary_file)

cat("BETA DIVERSITY ANALYSIS\n")
cat("=======================\n\n")
cat("Cohort: PRJEB27564 baseline\n")
cat("Samples:", nrow(taxa), "\n")
cat("Genus features:", ncol(taxa), "\n")
cat("Controls:", sum(md$parkinson == "Control"), "\n")
cat("PD:", sum(md$parkinson == "PD"), "\n")
cat("Distance: Bray-Curtis\n")
cat("Ordination: PCoA with additive correction\n")
cat("PCoA axis 1:", round(axis1_pct, 3), "%\n")
cat("PCoA axis 2:", round(axis2_pct, 3), "%\n")
cat("\nPERMANOVA:\n")
print(permanova)
cat("\nPERMDISP:\n")
print(dispersion_test)

sink()

cat("\n====================================\n")
cat("BETA DIVERSITY ANALYSIS COMPLETED\n")
cat("====================================\n")
cat("Results directory:", out_dir, "\n")
cat("PCoA axes:", round(axis1_pct, 2), "% and",
    round(axis2_pct, 2), "%\n")
cat("All output files saved successfully.\n")
