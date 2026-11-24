
seqtab <- readRDS(
  "data/processed/PRJEB27564/dada2_chimera_131/seqtab_nochim.rds"
)

meta <- read.delim(
  "results/taxonomy/PRJEB27564/abundance_matrices_131/matched_sample_metadata.tsv",
  stringsAsFactors = FALSE
)

out <- "results/statistics/PRJEB27564/alpha_diversity"

dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Analysis settings
target_depth <- 20000
set.seed(20261121)

# Match metadata
meta$Run <- as.character(meta$Run)
meta <- meta[match(rownames(seqtab), meta$Run), ]

if (anyNA(meta$Run)) {
  stop("Sample metadata matching failed.")
}

group <- ifelse(
  meta$parkinson == "yes", "PD",
  ifelse(meta$parkinson == "no", "Control", NA)
)

if (anyNA(group)) {
  stop("Unexpected group labels.")
}

depth <- rowSums(seqtab)
keep <- depth >= target_depth

cat("Original samples:", nrow(seqtab), "\n")
cat("Selected rarefaction depth:", target_depth, "\n")
cat("Samples retained:", sum(keep), "\n")
cat("Samples excluded:", sum(!keep), "\n")

cat("\nGroup distribution:\n")
print(table(group[keep]))

# Record excluded samples
excluded <- data.frame(
  SampleID = rownames(seqtab)[!keep],
  Group = group[!keep],
  Original_depth = depth[!keep]
)

write.table(
  excluded,
  file = file.path(out, "alpha_excluded_low_depth.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Rarefy without replacement
eligible <- seqtab[keep, , drop = FALSE]

rarefied <- matrix(
  0L,
  nrow = nrow(eligible),
  ncol = ncol(eligible),
  dimnames = dimnames(eligible)
)

for (i in seq_len(nrow(eligible))) {

  counts <- as.integer(eligible[i, ])

  # Expand ASV IDs into individual read labels
  read_labels <- rep(seq_along(counts), times = counts)

  # Randomly select exactly 20,000 reads
  sampled <- sample(
    read_labels,
    size = target_depth,
    replace = FALSE
  )

  rarefied[i, ] <- tabulate(
    sampled,
    nbins = ncol(eligible)
  )

  if (i %% 20 == 0) {
    cat("Rarefied samples:", i, "\n")
  }
}

if (any(rowSums(rarefied) != target_depth)) {
  stop("Rarefied sample depth validation failed.")
}

# Calculate diversity indices
calculate_diversity <- function(counts) {

  counts <- counts[counts > 0]
  p <- counts / sum(counts)

  observed <- length(counts)
  shannon <- -sum(p * log(p))
  simpson <- 1 - sum(p^2)

  pielou <- if (observed > 1) {
    shannon / log(observed)
  } else {
    NA_real_
  }

  c(
    Observed_ASVs = observed,
    Shannon = shannon,
    Simpson_1_minus_D = simpson,
    Pielou_evenness = pielou
  )
}

metrics <- t(apply(rarefied, 1, calculate_diversity))

alpha <- data.frame(
  SampleID = rownames(rarefied),
  Group = group[keep],
  Original_depth = depth[keep],
  Rarefaction_depth = rowSums(rarefied),
  metrics,
  row.names = NULL
)

# Save sample-wise diversity table
write.table(
  alpha,
  file = file.path(out, "alpha_diversity_20k.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Save rarefied count matrix
write.table(
  data.frame(
    SampleID = rownames(rarefied),
    as.data.frame(rarefied, check.names = FALSE),
    check.names = FALSE
  ),
  file = file.path(out, "rarefied_ASV_counts_20k.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

saveRDS(
  rarefied,
  file = file.path(out, "rarefied_ASV_counts_20k.rds")
)

# Group-wise descriptive statistics
indices <- c(
  "Observed_ASVs",
  "Shannon",
  "Simpson_1_minus_D",
  "Pielou_evenness"
)

group_summary <- do.call(rbind, lapply(indices, function(index) {

  do.call(rbind, lapply(c("Control", "PD"), function(g) {

    values <- alpha[alpha$Group == g, index]

    data.frame(
      Index = index,
      Group = g,
      Samples = length(values),
      Mean = mean(values),
      SD = sd(values),
      Median = median(values),
      Q1 = unname(quantile(values, 0.25)),
      Q3 = unname(quantile(values, 0.75))
    )
  }))
}))

write.table(
  group_summary,
  file = file.path(out, "alpha_group_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\nAlpha diversity completed.\n")
cat("Samples:", nrow(alpha), "\n")
cat("Rarefaction depth:", target_depth, "\n")

cat("\nGroup-wise diversity summary:\n")
print(group_summary, row.names = FALSE)

cat("\nOutputs saved in:", out, "\n")
