
suppressPackageStartupMessages({
  library(dada2)
})

filt_path <- "data/processed/PRJEB27564/dada2_filter_131"

qc_path <- "results/qc/PRJEB27564/dada2_error_models_131"

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
  stop("Expected 131 forward and 131 reverse filtered files.")
}

if (!identical(sample.names, sample.namesR)) {
  stop("Forward and reverse sample names do not match.")
}

write.table(
  data.frame(sample = sample.names),
  file = file.path(qc_path, "sample_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Samples:", length(sample.names), "\n")
cat("DADA2 version:", as.character(packageVersion("dada2")), "\n")

set.seed(1001)

cat("\nLearning forward-read error rates...\n")

errF <- learnErrors(
  filtFs,
  nbases = 1e8,
  randomize = TRUE,
  multithread = FALSE
)

saveRDS(errF, file.path(qc_path, "error_model_R1.rds"))

pdf(file.path(qc_path, "error_model_R1.pdf"))
plotErrors(errF, nominalQ = TRUE)
dev.off()

cat("\nLearning reverse-read error rates...\n")

errR <- learnErrors(
  filtRs,
  nbases = 1e8,
  randomize = TRUE,
  multithread = FALSE
)

saveRDS(errR, file.path(qc_path, "error_model_R2.rds"))

pdf(file.path(qc_path, "error_model_R2.pdf"))
plotErrors(errR, nominalQ = TRUE)
dev.off()

cat("\nBoth error models saved successfully.\n")
