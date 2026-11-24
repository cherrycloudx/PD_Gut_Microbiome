
suppressPackageStartupMessages({
  library(dada2)
})

in_path <- "data/processed/PRJEB27564/dada2_merge_131"
out_path <- "data/processed/PRJEB27564/dada2_chimera_131"
qc_path <- "results/qc/PRJEB27564/dada2_chimera_131"

dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_path, recursive = TRUE, showWarnings = FALSE)

seqtab <- readRDS(file.path(in_path, "seqtab_unfiltered.rds"))

if (nrow(seqtab) != 131) {
  stop("Expected 131 samples in the sequence table.")
}

if (ncol(seqtab) == 0 || is.null(rownames(seqtab))) {
  stop("Sequence table is empty or missing sample names.")
}

cat("Input sequence table:", nrow(seqtab), "samples x",
    ncol(seqtab), "ASVs\n")

reads_before <- rowSums(seqtab)

cat("Total merged reads:", sum(reads_before), "\n")
cat("\nStarting consensus chimera removal...\n")

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = FALSE,
  verbose = TRUE
)

reads_after <- rowSums(seqtab.nochim)

track <- data.frame(
  sample = rownames(seqtab),
  merged_reads = reads_before,
  nonchim_reads = reads_after,
  reads_removed = reads_before - reads_after,
  retained_pct = round(100 * reads_after / reads_before, 2)
)

write.table(
  track,
  file = file.path(qc_path, "chimera_tracking.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_table <- data.frame(
  stage = c("Before chimera removal", "After chimera removal"),
  samples = c(nrow(seqtab), nrow(seqtab.nochim)),
  ASVs = c(ncol(seqtab), ncol(seqtab.nochim)),
  total_reads = c(sum(reads_before), sum(reads_after))
)

write.table(
  summary_table,
  file = file.path(qc_path, "chimera_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  seqtab.nochim,
  file.path(out_path, "seqtab_nochim.rds")
)

cat("\nChimera removal completed.\n")
cat("Samples retained:", nrow(seqtab.nochim), "\n")
cat("ASVs before:", ncol(seqtab), "\n")
cat("ASVs after:", ncol(seqtab.nochim), "\n")
cat("Merged reads before:", sum(reads_before), "\n")
cat("Non-chimeric reads retained:", sum(reads_after), "\n")
cat("Overall read retention:",
    round(100 * sum(reads_after) / sum(reads_before), 2), "%\n")
