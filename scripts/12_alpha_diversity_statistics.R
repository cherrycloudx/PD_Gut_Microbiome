
alpha <- read.delim(
  "results/statistics/PRJEB27564/alpha_diversity/alpha_diversity_20k.tsv"
)

indices <- c(
  "Observed_ASVs",
  "Shannon",
  "Simpson_1_minus_D",
  "Pielou_evenness"
)

results <- lapply(indices, function(index) {

  pd <- alpha[alpha$Group == "PD", index]
  control <- alpha[alpha$Group == "Control", index]

  test <- wilcox.test(
    pd,
    control,
    exact = FALSE,
    conf.int = TRUE
  )

  # Rank-biserial effect size: positive means higher values in PD
  ranks <- rank(c(pd, control))
  n_pd <- length(pd)
  n_control <- length(control)

  U <- sum(ranks[seq_along(pd)]) -
    n_pd * (n_pd + 1) / 2

  rank_biserial <- 2 * U / (n_pd * n_control) - 1

  data.frame(
    Index = index,
    PD_median = median(pd),
    Control_median = median(control),
    Median_difference_PD_minus_Control =
      median(pd) - median(control),
    PD_mean = mean(pd),
    Control_mean = mean(control),
    Rank_biserial = rank_biserial,
    Wilcoxon_W = unname(test$statistic),
    P_value = test$p.value
  )
})

results <- do.call(rbind, results)

results$Adjusted_P_value <- p.adjust(
  results$P_value,
  method = "BH"
)

results$Interpretation <- ifelse(
  results$Adjusted_P_value < 0.05,
  "Statistically significant after BH correction",
  "Not statistically significant after BH correction"
)

out <- "results/statistics/PRJEB27564/alpha_diversity"

write.table(
  results,
  file = file.path(out, "alpha_diversity_group_tests.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Alpha diversity statistical comparison\n")
cat("PD samples:", sum(alpha$Group == "PD"), "\n")
cat("Control samples:", sum(alpha$Group == "Control"), "\n\n")

print(results, row.names = FALSE)

# Generate four boxplots
pdf(
  file.path(out, "alpha_diversity_group_comparison.pdf"),
  width = 10,
  height = 8
)

par(mfrow = c(2, 2))

for (index in indices) {

  boxplot(
    alpha[[index]] ~ alpha$Group,
    col = c("lightblue", "salmon"),
    main = index,
    xlab = "Group",
    ylab = index
  )

}

dev.off()

cat("\nResults and plots saved in:", out, "\n")
