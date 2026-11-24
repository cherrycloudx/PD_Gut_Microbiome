
suppressPackageStartupMessages({
  library(dada2)
})

seqtab_path <- "data/processed/PRJEB27564/dada2_chimera_131/seqtab_nochim.rds"

ref_path <- "references/databases/SILVA_138.2/silva_nr99_v138.2_toGenus_trainset.fa.gz"

out_path <- "results/taxonomy/PRJEB27564"

dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(seqtab_path)) {
  stop("Non-chimeric ASV table not found.")
}

if (!file.exists(ref_path)) {
  stop("SILVA reference database not found.")
}

seqtab <- readRDS(seqtab_path)

if (nrow(seqtab) != 131 || ncol(seqtab) != 7501) {
  stop("Unexpected sequence table dimensions. Check input.")
}

asv_seqs <- colnames(seqtab)

if (is.null(asv_seqs) || anyDuplicated(asv_seqs)) {
  stop("ASV sequences are missing or duplicated.")
}

asv_ids <- paste0("ASV_", seq_along(asv_seqs))

cat("Samples:", nrow(seqtab), "\n")
cat("ASVs:", length(asv_seqs), "\n")
cat("Reference: SILVA 138.2\n")
cat("\nStarting taxonomic assignment...\n")

tax.out <- assignTaxonomy(
  asv_seqs,
  ref_path,
  minBoot = 50,
  tryRC = TRUE,
  outputBootstraps = TRUE,
  taxLevels = c(
    "Kingdom", "Phylum", "Class",
    "Order", "Family", "Genus"
  ),
  multithread = FALSE,
  verbose = TRUE
)

saveRDS(
  tax.out,
  file.path(out_path, "taxonomy_output_with_bootstrap_SILVA138.2.rds")
)

tax <- tax.out$tax
boot <- tax.out$boot

if (nrow(tax) != length(asv_seqs)) {
  stop("Taxonomy result count does not match ASV count.")
}

rownames(tax) <- NULL
rownames(boot) <- NULL

taxonomy_table <- data.frame(
  ASV_ID = asv_ids,
  Sequence = asv_seqs,
  Total_reads = as.integer(colSums(seqtab)),
  as.data.frame(tax),
  check.names = FALSE
)

bootstrap_table <- data.frame(
  ASV_ID = asv_ids,
  as.data.frame(boot),
  check.names = FALSE
)

write.table(
  taxonomy_table,
  file = file.path(out_path, "ASV_taxonomy_SILVA138.2.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  bootstrap_table,
  file = file.path(out_path, "ASV_taxonomy_bootstrap_SILVA138.2.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

assigned <- vapply(
  seq_len(ncol(tax)),
  function(j) sum(!is.na(tax[, j]) & nzchar(tax[, j])),
  numeric(1)
)

rank_summary <- data.frame(
  Rank = colnames(tax),
  ASVs_assigned = assigned,
  ASVs_unassigned = nrow(tax) - assigned,
  Assigned_percentage = round(100 * assigned / nrow(tax), 2)
)

write.table(
  rank_summary,
  file = file.path(out_path, "taxonomy_rank_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  list(taxonomy = tax, bootstrap = boot),
  file = file.path(out_path, "taxonomy_SILVA138.2.rds")
)

cat("\nTaxonomic assignment completed.\n")
cat("ASVs classified:", nrow(tax), "\n\n")
print(rank_summary)
