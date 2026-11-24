
suppressPackageStartupMessages({
  library(dada2)
})

seqtab <- readRDS(
  "data/processed/PRJEB27564/dada2_chimera_131/seqtab_nochim.rds"
)

tax <- read.delim(
  "results/taxonomy/PRJEB27564/ASV_taxonomy_SILVA138.2.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

meta <- read.delim(
  "data/metadata/PRJEB27564_Romano_matched_131.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

out_path <- "results/taxonomy/PRJEB27564/abundance_matrices_131"
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

# Validate ASV correspondence
idx <- match(colnames(seqtab), tax$Sequence)

if (anyNA(idx) || length(idx) != nrow(tax)) {
  stop("ASV sequences do not match the taxonomy table.")
}

tax <- tax[idx, ]

if (nrow(seqtab) != 131 || ncol(seqtab) != 7501) {
  stop("Unexpected sequence table dimensions.")
}

# Validate metadata
if (!all(c("Run", "parkinson") %in% names(meta))) {
  stop("Required metadata columns Run and parkinson are missing.")
}

meta$Run <- as.character(meta$Run)

if (anyDuplicated(meta$Run)) {
  stop("Duplicate Run IDs in metadata.")
}

meta <- meta[match(rownames(seqtab), meta$Run), ]

if (anyNA(meta$Run)) {
  stop("Some ASV table samples are missing from metadata.")
}

write.table(
  meta,
  file = file.path(out_path, "matched_sample_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

asv_by_sample <- t(seqtab)
genus <- trimws(as.character(tax$Genus))
family <- trimws(as.character(tax$Family))

genus_assigned <- !is.na(genus) & nzchar(genus)
family_assigned <- !is.na(family) & nzchar(family)

# 1. Genus-level count matrix
genus_counts <- t(rowsum(
  asv_by_sample[genus_assigned, , drop = FALSE],
  group = genus[genus_assigned],
  reorder = TRUE
))

write.table(
  data.frame(SampleID = rownames(genus_counts),
             as.data.frame(genus_counts, check.names = FALSE),
             check.names = FALSE),
  file = file.path(out_path, "genus_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 2. Genus relative abundance using all non-chimeric reads
total_reads <- rowSums(seqtab)

genus_relative_total <- sweep(
  genus_counts, 1, total_reads[rownames(genus_counts)], "/"
) * 100

write.table(
  data.frame(SampleID = rownames(genus_relative_total),
             as.data.frame(genus_relative_total, check.names = FALSE),
             check.names = FALSE),
  file = file.path(out_path, "genus_relative_abundance_total_reads.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 3. Genus relative abundance among genus-assigned reads only
assigned_reads <- rowSums(genus_counts)

genus_relative_assigned <- sweep(
  genus_counts, 1, assigned_reads, "/"
) * 100

write.table(
  data.frame(SampleID = rownames(genus_relative_assigned),
             as.data.frame(genus_relative_assigned, check.names = FALSE),
             check.names = FALSE),
  file = file.path(out_path, "genus_relative_abundance_assigned_reads.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 4. Preserve ASVs without genus assignments
unassigned <- tax[!genus_assigned, ]

write.table(
  unassigned,
  file = file.path(out_path, "ASVs_unassigned_at_genus.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 5. Build a feature matrix using each ASV's deepest assigned rank
assigned_rank <- ifelse(
  genus_assigned, "Genus",
  ifelse(family_assigned, "Family_unresolved", "Unclassified_ASV")
)

taxon_label <- ifelse(
  genus_assigned, genus,
  ifelse(family_assigned, family, tax$ASV_ID)
)

feature_id <- ifelse(
  genus_assigned, paste0("Genus::", taxon_label),
  ifelse(
    family_assigned,
    paste0("Family_unclassified::", taxon_label),
    paste0("Unclassified_ASV::", taxon_label)
  )
)

feature_counts <- t(rowsum(
  asv_by_sample,
  group = feature_id,
  reorder = TRUE
))

write.table(
  data.frame(SampleID = rownames(feature_counts),
             as.data.frame(feature_counts, check.names = FALSE),
             check.names = FALSE),
  file = file.path(out_path, "deepest_rank_taxonomic_feature_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

feature_map_raw <- data.frame(
  feature_id = feature_id,
  assigned_rank = assigned_rank,
  taxon_label = taxon_label,
  ASV_ID = tax$ASV_ID,
  total_reads = as.numeric(colSums(seqtab))
)

feature_map <- aggregate(
  cbind(ASV_count = rep(1, nrow(feature_map_raw)),
        total_reads) ~ feature_id + assigned_rank + taxon_label,
  data = feature_map_raw,
  FUN = sum
)

write.table(
  feature_map,
  file = file.path(out_path, "taxonomic_feature_map.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 6. Overall genus abundance summary
genus_total <- sort(colSums(genus_counts), decreasing = TRUE)

genus_summary <- data.frame(
  Genus = names(genus_total),
  Total_reads = as.numeric(genus_total),
  Percentage_of_all_reads = round(
    100 * as.numeric(genus_total) / sum(seqtab), 3
  )
)

write.table(
  genus_summary,
  file = file.path(out_path, "genus_overall_abundance.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Abundance matrices completed.\n")
cat("Samples:", nrow(seqtab), "\n")
cat("Genus features:", ncol(genus_counts), "\n")
cat("Deepest-rank features:", ncol(feature_counts), "\n")
cat("Genus-assigned reads:", sum(genus_counts), "\n")
cat("Total reads:", sum(seqtab), "\n")
cat("Genus-assigned read percentage:",
    round(100 * sum(genus_counts) / sum(seqtab), 2), "%\n")

cat("\nTop 15 genera:\n")
print(head(genus_summary, 15))

cat("\nOutputs saved in:", out_path, "\n")
