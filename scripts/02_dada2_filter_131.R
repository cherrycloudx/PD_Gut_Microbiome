suppressPackageStartupMessages(library(dada2))

project_dir <- path.expand("~/PD_Gut_Microbiome")
setwd(project_dir)

run_file <- "data/metadata/PRJEB27564_Romano_matched_131_runs.txt"

input_dir <- "data/processed/PRJEB27564/cutadapt_trimmed_131"
output_dir <- "data/processed/PRJEB27564/dada2_filter_131"

stats_file <- "results/qc/PRJEB27564/dada2_filter_stats_131.tsv"

runs <- trimws(readLines(run_file))
runs <- runs[nzchar(runs)]

if (length(runs) != 131 || anyDuplicated(runs)) {
    stop("Sample list must contain 131 unique IDs.")
}

fnFs <- file.path(input_dir, paste0(runs, "_R1_trimmed.fastq.gz"))
fnRs <- file.path(input_dir, paste0(runs, "_R2_trimmed.fastq.gz"))

if (any(!file.exists(fnFs)) || any(!file.exists(fnRs))) {
    stop("One or more trimmed FASTQ files are missing.")
}

if (dir.exists(output_dir)) {
    stop("Output directory already exists. Stopping to prevent overwriting.")
}

if (file.exists(stats_file)) {
    stop("Statistics file already exists. Stopping to prevent overwriting.")
}

dir.create(output_dir, recursive = TRUE)

filtFs <- file.path(output_dir, paste0(runs, "_R1_filtered.fastq.gz"))
filtRs <- file.path(output_dir, paste0(runs, "_R2_filtered.fastq.gz"))

cat("Starting DADA2 filtering for", length(runs), "samples.\n")

filter_stats <- dada2::filterAndTrim(
    fwd = fnFs,
    filt = filtFs,
    rev = fnRs,
    filt.rev = filtRs,
    truncLen = c(280, 220),
    maxN = 0,
    maxEE = c(2, 4),
    truncQ = 2,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = FALSE,
    verbose = TRUE
)

stats <- as.data.frame(filter_stats)
stats$sample <- runs
stats$percent_retained <- round(
    100 * stats$reads.out / stats$reads.in, 2
)

stats <- stats[, c(
    "sample", "reads.in", "reads.out", "percent_retained"
)]

write.table(
    stats,
    stats_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("\nFiltering completed.\n")
cat("Samples:", nrow(stats), "\n")
cat("Input read pairs:", sum(stats$reads.in), "\n")
cat("Retained read pairs:", sum(stats$reads.out), "\n")
cat("Overall retention:",
    round(100 * sum(stats$reads.out) / sum(stats$reads.in), 2),
    "%\n")

cat("\nPer-sample retention summary:\n")
print(summary(stats$percent_retained))

cat("\nStatistics saved to:", stats_file, "\n")
