library(vegan)

ab <- read.delim(
  "results/taxonomy/PRJEB27564/abundance_matrices_131/genus_relative_abundance_assigned_reads.tsv",
  check.names = FALSE
)

md <- read.delim(
  "results/taxonomy/PRJEB27564/abundance_matrices_131/matched_sample_metadata.tsv",
  check.names = FALSE
)

md <- md[match(ab$SampleID, md$Run), ]

md$date <- as.Date(
  ifelse(md$Collection_Date == "0000-00-00",
         NA, md$Collection_Date)
)

keep <- !is.na(md$date)
md <- md[keep, ]

md$Group <- factor(
  ifelse(md$parkinson == "yes", "PD", "Control"),
  levels = c("Control", "PD")
)

md$Month <- factor(
  format(md$date, "%B"),
  levels = month.name
)

taxa <- ab[keep, -1, drop = FALSE]
taxa[] <- lapply(taxa, as.numeric)
rownames(taxa) <- md$Run

bray <- vegdist(taxa, method = "bray")
pcoa <- cmdscale(bray, k = 2, eig = TRUE, add = TRUE)

scores <- as.data.frame(pcoa$points)
axis1 <- pcoa$eig[1] / sum(pcoa$eig[pcoa$eig > 0]) * 100
axis2 <- pcoa$eig[2] / sum(pcoa$eig[pcoa$eig > 0]) * 100

out <- "results/statistics/PRJEB27564/beta_diversity"

pdf(
  file.path(out, "Bray_Curtis_PCoA_128_date_complete.pdf"),
  width = 9,
  height = 7
)

par(mar = c(5, 5, 4, 2))

group_colors <- c("Control" = "#2878B5", "PD" = "#D95F02")

plot(
  scores$V1, scores$V2,
  col = group_colors[md$Group],
  pch = 19,
  cex = 1.1,
  xlab = paste0("PCoA1 (", round(axis1, 2), "%)"),
  ylab = paste0("PCoA2 (", round(axis2, 2), "%)"),
  main = "Bray-Curtis PCoA: PD vs Control\n128 samples with valid collection dates"
)

legend(
  "topright",
  legend = names(group_colors),
  col = group_colors,
  pch = 19,
  bty = "n"
)

month_colors <- setNames(
  hcl.colors(12, palette = "Dark 3"),
  month.name
)

plot(
  scores$V1, scores$V2,
  col = month_colors[as.character(md$Month)],
  pch = 19,
  cex = 1.1,
  xlab = paste0("PCoA1 (", round(axis1, 2), "%)"),
  ylab = paste0("PCoA2 (", round(axis2, 2), "%)"),
  main = "Bray-Curtis PCoA: Collection Month"
)

legend(
  "right",
  legend = levels(md$Month),
  col = month_colors,
  pch = 19,
  cex = 0.75,
  bty = "n"
)

dev.off()

cat("PCoA figures generated successfully.\n")
cat("Samples:", nrow(md), "\n")
cat("Axis 1:", round(axis1, 2), "%\n")
cat("Axis 2:", round(axis2, 2), "%\n")
cat("Output:", file.path(out, "Bray_Curtis_PCoA_128_date_complete.pdf"), "\n")
