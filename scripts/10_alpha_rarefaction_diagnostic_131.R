
seqtab <- readRDS(
  "data/processed/PRJEB27564/dada2_chimera_131/seqtab_nochim.rds"
)

meta <- read.delim(
  "results/taxonomy/PRJEB27564/abundance_matrices_131/matched_sample_metadata.tsv",
  stringsAsFactors = FALSE
)

out <- "results/statistics/PRJEB27564/alpha_diversity"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

if (nrow(seqtab) != 131 || is.null(rownames(seqtab))) {
  stop("Unexpected sequence table.")
}

meta$Run <- as.character(meta$Run)
meta <- meta[match(rownames(seqtab), meta$Run), ]

if (anyNA(meta$Run)) stop("Metadata matching failed.")

group <- ifelse(meta$parkinson == "yes", "PD",
                ifelse(meta$parkinson == "no", "Control", NA))

if (anyNA(group)) stop("Unexpected group labels.")

depth <- rowSums(seqtab)

# Expected number of ASVs in a random subsample of n reads.
expected_richness <- function(counts, n) {
  counts <- counts[counts > 0]
  N <- sum(counts)

  if (n >= N) return(length(counts))
  if (n <= 0) return(0)

  possible <- (N - counts) >= n
  p_absent <- numeric(length(counts))

  p_absent[possible] <- exp(
    lchoose(N - counts[possible], n) - lchoose(N, n)
  )

  sum(1 - p_absent)
}

# Candidate common depths
fixed_depths <- c(
  2500, 5000, 10000, 15000, 20000,
  25000, 30000, 40000, 50000
)

curve_list <- vector("list", nrow(seqtab))
plateau_list <- vector("list", nrow(seqtab))

for (i in seq_len(nrow(seqtab))) {

  N <- depth[i]
  counts <- seqtab[i, ]

  depths <- unique(sort(c(
    round(seq(1000, N, length.out = 30)),
    fixed_depths[fixed_depths <= N],
    N
  )))

  richness <- vapply(
    depths,
    function(n) expected_richness(counts, n),
    numeric(1)
  )

  curve_list[[i]] <- data.frame(
    sample = rownames(seqtab)[i],
    group = group[i],
    depth = depths,
    expected_ASVs = richness
  )

  r80 <- expected_richness(counts, floor(0.8 * N))
  r100 <- sum(counts > 0)

  plateau_list[[i]] <- data.frame(
    sample = rownames(seqtab)[i],
    group = group[i],
    sequencing_depth = N,
    observed_ASVs = r100,
    expected_ASVs_at_80pct = r80,
    richness_gain_last_20pct = round(
      100 * (r100 - r80) / r100, 2
    )
  )
}

curves <- do.call(rbind, curve_list)
plateau <- do.call(rbind, plateau_list)

write.table(
  curves,
  file = file.path(out, "ASV_rarefaction_curves.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  plateau,
  file = file.path(out, "sample_rarefaction_plateau.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Candidate-depth sample retention
thresholds <- do.call(rbind, lapply(fixed_depths, function(d) {
  keep <- depth >= d

  data.frame(
    depth = d,
    samples = sum(keep),
    PD = sum(keep & group == "PD"),
    controls = sum(keep & group == "Control")
  )
}))

write.table(
  thresholds,
  file = file.path(out, "rarefaction_depth_retention.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Plot individual curves and group median curves
pdf(file.path(out, "ASV_rarefaction_curves.pdf"),
    width = 10, height = 7)

plot(
  NA,
  xlim = c(0, max(depth)),
  ylim = c(0, max(plateau$observed_ASVs) * 1.05),
  xlab = "Subsampled sequencing depth",
  ylab = "Expected observed ASV richness",
  main = "ASV rarefaction curves: 131 samples"
)

for (i in seq_len(nrow(seqtab))) {
  z <- curves[curves$sample == rownames(seqtab)[i], ]

  colour <- if (group[i] == "PD") {
    adjustcolor("firebrick", alpha.f = 0.16)
  } else {
    adjustcolor("steelblue", alpha.f = 0.16)
  }

  lines(z$depth, z$expected_ASVs, col = colour, lwd = 0.7)
}

medians <- list()

for (g in c("PD", "Control")) {
  for (d in fixed_depths) {
    z <- curves[curves$group == g & curves$depth == d, ]

    if (nrow(z) > 0) {
      medians[[length(medians) + 1]] <- data.frame(
        group = g,
        depth = d,
        median_richness = median(z$expected_ASVs),
        samples = nrow(z)
      )
    }
  }
}

medians <- do.call(rbind, medians)

for (g in c("PD", "Control")) {
  z <- medians[medians$group == g, ]
  z <- z[order(z$depth), ]

  colour <- if (g == "PD") "firebrick" else "steelblue"

  lines(z$depth, z$median_richness,
        col = colour, lwd = 3)

  points(z$depth, z$median_richness,
         col = colour, pch = 19)
}

abline(v = 20000, lty = 2, lwd = 2)

legend(
  "bottomright",
  legend = c("PD", "Control", "20,000-read candidate"),
  col = c("firebrick", "steelblue", "black"),
  lty = c(1, 1, 2),
  lwd = c(2, 2, 2),
  bty = "n"
)

dev.off()

write.table(
  medians,
  file = file.path(out, "group_rarefaction_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Rarefaction diagnostics completed.\n")
cat("Samples:", nrow(seqtab), "\n")
cat("\nCandidate-depth retention:\n")
print(thresholds, row.names = FALSE)

cat("\nRichness gain during final 20% of sequencing:\n")
print(summary(plateau$richness_gain_last_20pct))

cat("\nOutputs saved in:", out, "\n")
