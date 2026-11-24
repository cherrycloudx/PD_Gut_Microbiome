
suppressPackageStartupMessages({
  library(dada2)
})

filt_path <- "data/processed/PRJEB27564/dada2_filter_131"

error_path <- "results/qc/PRJEB27564/dada2_error_models_131"

out_path <- "data/processed/PRJEB27564/dada2_inference_131"

qc_path <- "results/qc/PRJEB27564/dada2_inference_131"

dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_path, recursive = TRUE, showWarnings = FALSE)

filtFs <- sort(list.files(
  filt_path,
  pattern = "_R1_filtered\\.fastq\\.gz$",
  full.names = TRUE
))

filtRs <- sort(list.files(
  filt_path,
  pattern = "_R2_filtered\\.fastq\\.gz$",
  full.names = TRUE
))

sample.names <- sub(
  "_R1_filtered\\.fastq\\.gz$",
  "",
  basename(filtFs)
)

sample.namesR <- sub(
  "_R2_filtered\\.fastq\\.gz$",
  "",
  basename(filtRs)
)

if (length(filtFs) != 131 || length(filtRs) != 131) {
  stop("Expected 131 forward and 131 reverse files.")
}

if (!identical(sample.names, sample.namesR)) {
  stop("Forward and reverse sample names do not match.")
}

errF <- readRDS(file.path(error_path, "error_model_R1.rds"))
errR <- readRDS(file.path(error_path, "error_model_R2.rds"))

cat("Starting forward dereplication...\n")

derepFs <- derepFastq(filtFs, verbose = TRUE)
names(derepFs) <- sample.names

saveRDS(derepFs, file.path(out_path, "derep_R1.rds"))

cat("Starting reverse dereplication...\n")

derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepRs) <- sample.names

saveRDS(derepRs, file.path(out_path, "derep_R2.rds"))

cat("Starting forward ASV inference...\n")

dadaFs <- dada(
  derepFs,
  err = errF,
  multithread = FALSE,
  verbose = TRUE
)

saveRDS(dadaFs, file.path(out_path, "dada_R1.rds"))

cat("Starting reverse ASV inference...\n")

dadaRs <- dada(
  derepRs,
  err = errR,
  multithread = FALSE,
  verbose = TRUE
)

saveRDS(dadaRs, file.path(out_path, "dada_R2.rds"))

getN <- function(x) sum(getUniques(x))

track <- data.frame(
  sample = sample.names,
  filtered_R1 = sapply(derepFs, function(x) sum(x$uniques)),
  filtered_R2 = sapply(derepRs, function(x) sum(x$uniques)),
  denoised_R1 = sapply(dadaFs, getN),
  denoised_R2 = sapply(dadaRs, getN)
)

write.table(
  track,
  file = file.path(qc_path, "asv_inference_tracking.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nASV inference completed for all 131 samples.\n")
