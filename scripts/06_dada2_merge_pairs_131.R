
suppressPackageStartupMessages({
  library(dada2)
})

base_path <- "data/processed/PRJEB27564"

in_path <- file.path(base_path, "dada2_inference_131")
out_path <- file.path(base_path, "dada2_merge_131")

qc_path <- "results/qc/PRJEB27564/dada2_merge_131"

dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_path, recursive = TRUE, showWarnings = FALSE)

sample.names <- read.delim(
  "results/qc/PRJEB27564/dada2_error_models_131/sample_manifest.tsv"
)$sample

derepFs <- readRDS(file.path(in_path, "derep_R1.rds"))
derepRs <- readRDS(file.path(in_path, "derep_R2.rds"))
dadaFs <- readRDS(file.path(in_path, "dada_R1.rds"))
dadaRs <- readRDS(file.path(in_path, "dada_R2.rds"))

objects <- list(derepFs, derepRs, dadaFs, dadaRs)

for (i in seq_along(objects)) {
  if (is.null(names(objects[[i]])) ||
      !all(sample.names %in% names(objects[[i]]))) {
    stop("Sample names do not match in saved DADA2 objects.")
  }
}

derepFs <- derepFs[sample.names]
derepRs <- derepRs[sample.names]
dadaFs <- dadaFs[sample.names]
dadaRs <- dadaRs[sample.names]

cat("Merging paired reads for 131 samples...\n")

mergers <- mergePairs(
  dadaFs,
  derepFs,
  dadaRs,
  derepRs,
  verbose = TRUE
)

names(mergers) <- sample.names

saveRDS(
  mergers,
  file.path(out_path, "merged_pairs_131.rds")
)

filtered_pairs <- vapply(
  derepFs,
  function(x) sum(x$uniques),
  numeric(1)
)

denoised_R1 <- vapply(
  dadaFs,
  function(x) sum(getUniques(x)),
  numeric(1)
)

denoised_R2 <- vapply(
  dadaRs,
  function(x) sum(getUniques(x)),
  numeric(1)
)

merged_pairs <- vapply(
  mergers,
  function(x) {
    if (is.null(x) || nrow(x) == 0) {
      return(0)
    }
    sum(x$abundance)
  },
  numeric(1)
)

track <- data.frame(
  sample = sample.names,
  filtered_pairs = filtered_pairs,
  denoised_R1 = denoised_R1,
  denoised_R2 = denoised_R2,
  merged_pairs = merged_pairs,
  merge_pct_filtered = round(
    100 * merged_pairs / filtered_pairs, 2
  )
)

write.table(
  track,
  file = file.path(qc_path, "merge_tracking.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nConstructing ASV sequence table...\n")

seqtab <- makeSequenceTable(mergers)

saveRDS(
  seqtab,
  file.path(out_path, "seqtab_unfiltered.rds")
)

cat("Samples with merged reads:",
    sum(merged_pairs > 0), "/131\n")

cat("Total merged read pairs:",
    sum(merged_pairs), "\n")

cat("ASV table dimensions:",
    nrow(seqtab), "samples x",
    ncol(seqtab), "ASVs\n")

cat("\nPaired-read merging completed.\n")
