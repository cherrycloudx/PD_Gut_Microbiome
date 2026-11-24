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

md$parkinson <- factor(
  ifelse(md$parkinson == "yes", "PD", "Control"),
  levels = c("Control", "PD")
)

md$month <- factor(format(md$date, "%m"))

taxa <- ab[keep, -1, drop = FALSE]
taxa[] <- lapply(taxa, as.numeric)
rownames(taxa) <- md$Run

bray <- vegdist(taxa, method = "bray")

cat("\nSAMPLE COUNTS\n")
cat("Total:", nrow(md), "\n")
print(table(md$month, md$parkinson))

cat("\nRESTRICTED PERMANOVA\n")
cat("Permutations restricted within collection month\n")

set.seed(20261122)

result <- adonis2(
  bray ~ parkinson,
  data = md,
  permutations = 9999,
  method = "bray",
  strata = md$month
)

print(result)

out <- "results/statistics/PRJEB27564/collection_date_sensitivity"

write.table(
  as.data.frame(result),
  file.path(out, "PERMANOVA_month_stratified.tsv"),
  sep = "\t",
  quote = FALSE
)

cat("\nRestricted PERMANOVA completed.\n")
cat("Results saved in:", out, "\n")
