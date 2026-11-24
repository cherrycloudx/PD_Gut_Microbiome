library(vegan)

ab_file <- "results/taxonomy/PRJEB27564/abundance_matrices_131/genus_relative_abundance_assigned_reads.tsv"
md_file <- "results/taxonomy/PRJEB27564/abundance_matrices_131/matched_sample_metadata.tsv"

out <- "results/statistics/PRJEB27564/collection_date_sensitivity"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

ab <- read.delim(ab_file, check.names = FALSE)
md <- read.delim(md_file, check.names = FALSE)

md <- md[match(ab$SampleID, md$Run), ]

if (anyNA(md$Run)) stop("Sample matching failed.")

md$date <- as.Date(
  ifelse(md$Collection_Date == "0000-00-00",
         NA_character_, md$Collection_Date)
)

md$parkinson <- factor(
  ifelse(md$parkinson == "yes", "PD", "Control"),
  levels = c("Control", "PD")
)

md$month <- factor(
  format(md$date, "%m"),
  levels = sprintf("%02d", 1:12)
)

month_num <- as.integer(format(md$date, "%m"))

md$quarter <- factor(
  paste0("Q", ((month_num - 1) %/% 3) + 1),
  levels = paste0("Q", 1:4)
)

cat("\nDATE VALIDATION\n")
cat("Total samples:", nrow(md), "\n")
cat("Valid dates:", sum(!is.na(md$date)), "\n")
cat("Missing dates:", sum(is.na(md$date)), "\n")

cat("\nMONTH BY DISEASE GROUP\n")
month_tab <- table(md$month, md$parkinson)
print(month_tab)

cat("\nTesting month-group association...\n")
set.seed(20261122)

month_test <- fisher.test(
  month_tab,
  simulate.p.value = TRUE,
  B = 99999
)

print(month_test)

cat("\nQUARTER BY DISEASE GROUP\n")
quarter_tab <- table(md$quarter, md$parkinson)
print(quarter_tab)

set.seed(20261122)

quarter_test <- fisher.test(
  quarter_tab,
  simulate.p.value = TRUE,
  B = 99999
)

print(quarter_test)

write.table(
  as.data.frame.matrix(month_tab),
  file.path(out, "month_group_counts.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

write.table(
  as.data.frame.matrix(quarter_tab),
  file.path(out, "quarter_group_counts.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

date_tests <- data.frame(
  Test = c("Month versus disease group",
           "Quarter versus disease group"),
  P_value = c(month_test$p.value, quarter_test$p.value)
)

write.table(
  date_tests,
  file.path(out, "collection_date_group_tests.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\nPreparing complete-date abundance matrix...\n")

taxa <- ab[, -1, drop = FALSE]
taxa[] <- lapply(taxa, as.numeric)
rownames(taxa) <- ab$SampleID

md_cc <- md[!is.na(md$date), ]
taxa_cc <- taxa[md_cc$Run, , drop = FALSE]

if (anyNA(taxa_cc)) stop("Missing abundance values.")

bray_cc <- vegdist(taxa_cc, method = "bray")

cat("\nPERMANOVA: COMPLETE-DATE SAMPLES, UNADJUSTED\n")

set.seed(20261122)

permanova_unadj <- adonis2(
  bray_cc ~ parkinson,
  data = md_cc,
  permutations = 9999,
  by = "margin"
)

print(permanova_unadj)

cat("\nPERMANOVA: ADJUSTED FOR COLLECTION QUARTER\n")

set.seed(20261122)

permanova_adj <- adonis2(
  bray_cc ~ quarter + parkinson,
  data = md_cc,
  permutations = 9999,
  by = "margin"
)

print(permanova_adj)

write.table(
  as.data.frame(permanova_unadj),
  file.path(out, "PERMANOVA_complete_date_unadjusted.tsv"),
  sep = "\t", quote = FALSE
)

write.table(
  as.data.frame(permanova_adj),
  file.path(out, "PERMANOVA_adjusted_quarter.tsv"),
  sep = "\t", quote = FALSE
)

cat("\nPERMDISP: COMPLETE-DATE SAMPLES\n")

disp <- betadisper(
  bray_cc,
  md_cc$parkinson,
  type = "median",
  bias.adjust = TRUE
)

set.seed(20261122)

disp_test <- permutest(disp, permutations = 9999)

print(disp_test)

write.table(
  as.data.frame(disp_test$tab),
  file.path(out, "PERMDISP_complete_date.tsv"),
  sep = "\t", quote = FALSE
)

write.table(
  md_cc[, c("Run", "host_subject_id", "parkinson",
            "Collection_Date", "month", "quarter")],
  file.path(out, "complete_date_sample_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\n====================================\n")
cat("COLLECTION DATE SENSITIVITY COMPLETED\n")
cat("====================================\n")
cat("Valid-date samples:", nrow(md_cc), "\n")
cat("Results saved in:", out, "\n")
