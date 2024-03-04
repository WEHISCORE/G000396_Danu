# NOTE: This is to replace the 'Time course analysis' section in DE_analysis.R

# Time course DE analysis (polynomial fit) -------------------------------------

Timepoint <- as.integer(y$samples$timepoint)
# NOTE: The coefficients X1 and X2 represent the linear and quadratic time
#       effect, respectively.
X <- poly(Timepoint, degree = 2)
Group <- relevel(y$samples$cell_line, "WT")
design_tc <- model.matrix(~Group * X)
colnames(design_tc) <- make.names(sub("Group", "", colnames(design_tc)))

fit_tc <- voomLmFit(
  y,
  design_tc,
  block = y$samples$cell_line_rep,
  sample.weights = TRUE,
  plot = TRUE)
fit_tc <- eBayes(fit_tc)

# NOTE: This is an ANOVA-like test to identify different (linear) time trends
#       in any one of the KOs vs. WT.
tt_linear <- topTable(
  fit_tc,
  coef = c("GID1KO.X1", "GID2KO.X1", "GID7KO.X1", "GID8KO.X1", "GID9KO.X1"),
  number = Inf,
  p.value = 0.05)
dim(tt_linear)

cm_tc <- makeContrasts(
  # Different linear trends for KOs vs. WTs
  GID1KO_vs_WT.linear = GID1KO.X1,
  GID2KO_vs_WT.linear = GID2KO.X1,
  GID7KO_vs_WT.linear = GID7KO.X1,
  GID8KO_vs_WT.linear = GID8KO.X1,
  GID9KO_vs_WT.linear = GID9KO.X1,
  KOs_vs_WT.linear =
    (GID1KO.X1 + GID2KO.X1 + GID7KO.X1 + GID8KO.X1 + GID9KO.X1) / 5,

  # Different quadratic trends for KOs vs. WTs
  GID1KO_vs_WT.quadratic = GID1KO.X2,
  GID2KO_vs_WT.quadratic = GID2KO.X2,
  GID7KO_vs_WT.quadratic = GID7KO.X2,
  GID8KO_vs_WT.quadratic = GID8KO.X2,
  GID9KO_vs_WT.quadratic = GID9KO.X2,
  KOs_vs_WT.quadratic =
    (GID1KO.X2 + GID2KO.X2 + GID7KO.X2 + GID8KO.X2 + GID9KO.X2) / 5,

  levels = design_tc)
cfit_tc <- contrasts.fit(fit_tc, cm_tc)
cfit_tc <- eBayes(cfit_tc)
summary(decideTests(cfit_tc))

# Outputs of time course DE analysis -------------------------------------------

# DEG lists as CSVs
dir.create(here("output", "timecourse"))
l_deg_tc_summary_df <- lapply(colnames(cfit_tc), function(j) {
  message(j)
  tt <- topTable(cfit_tc, j, number = Inf, sort.by = "none")
  if (grepl("KOs", j)) {
    # Add pairwise logFCs to the 'average contrast' to enable post-hoc
    # filtering of genes based on 'consistent' logFCs in the pairwise
    # comparisons.
    pairwise_coefs <- sapply(
      setdiff(levels(y$samples$cell_line), "WT"),
      function(cl) gsub("KOs", cl, j))
    n_pairwise <- 5
    names(pairwise_coefs) <- paste0(pairwise_coefs, ".logFC")
    pairwise_lfcCons <- lapply(pairwise_coefs, function(coef) {
      topTable(cfit_tc, coef = coef, sort.by = "none", n = Inf)$logFC
    })
    tt <- cbind(tt, do.call(cbind, pairwise_lfcCons))

    # Summarise the number of DEGs in each comparison.
    # NOTE: Mimics summary(decideTests(cfit)) value and then repeats with added
    #       'consistency of pairwise logFCs' constraint (i.e. 'lfcCon'.)
    deg_summary_df <- data.frame(
      Down = c(
        sum(tt$adj.P.Val < 0.05 & tt$logFC < 0),
        sum(
          tt$adj.P.Val < 0.05 &
            tt$logFC < 0 &
            rowSums(do.call(cbind, pairwise_lfcCons) < 0) == n_pairwise)),
      NotSig = c(
        sum(tt$adj.P.Val >= 0.05),
        sum(
          tt$adj.P.Val > 0.05 |
            (tt$adj.P.Val < 0.05 &
               tt$logFC < 0 &
               rowSums(do.call(cbind, pairwise_lfcCons) < 0) < n_pairwise) |
            tt$adj.P.Val < 0.05 &
            tt$logFC > 0 &
            rowSums(do.call(cbind, pairwise_lfcCons) > 0) < n_pairwise)),
      Up = c(
        sum(tt$adj.P.Val < 0.05 & tt$logFC > 0),
        sum(
          tt$adj.P.Val < 0.05 &
            tt$logFC > 0 &
            rowSums(do.call(cbind, pairwise_lfcCons) > 0) == n_pairwise)),
      row.names = c(j, paste0(j, ".lfcCon")))
  } else {
    # Summarise the number of DEGs in each comparison
    # NOTE: Mimics summary(decideTests(cfit)) value.
    deg_summary_df <- data.frame(
      Down = sum(tt$adj.P.Val < 0.05 & tt$logFC < 0),
      NotSig = sum(tt$adj.P.Val >= 0.05 ),
      Up = sum(tt$adj.P.Val < 0.05 & tt$logFC > 0),
      row.names = j)
  }

  # Write topTable results to CSV
  write.csv(
    flattenDF(tt[order(tt$adj.P.Val), ]),
    here("output", "timecourse", paste0(j, ".timecourse_DEGs.csv")))

  return(deg_summary_df)
})
# Summarise the number of DEGs in each comparison
do.call(rbind, l_deg_tc_summary_df)

# Plots
lcpm <- cpm(y, log = TRUE)
lcpm_fit_tc <- coef(fit_tc) %*% t(design_tc)
dir.create(here("output", "timecourse"))
tt <- read.csv(
  here("output", "timecourse", "KOs_vs_WT.linear.timecourse_DEGs.csv"),
  row.names = 1)

# NOTE: Identify and plot DEGs where the WT has (linear) trend with the
#       opposite sign to the KOs.
# NOTE: This uses the average of the pairwise logFCs.
same_sign <- sign(coef(fit_tc)[rownames(tt), "X1"] + tt$logFC) ==
  sign(coef(fit_tc)[rownames(tt), "X1"])
opposite_sign <- !same_sign
# NOTE: This stricter definition checks all pairwise comparisons.
same_sign_strict <-
  sign(coef(fit_tc)[rownames(tt), "X1"] + tt$GID1KO_vs_WT.linear.logFC) ==
  sign(coef(fit_tc)[rownames(tt), "X1"]) &
  sign(coef(fit_tc)[rownames(tt), "X1"] + tt$GID2KO_vs_WT.linear.logFC) ==
  sign(coef(fit_tc)[rownames(tt), "X1"]) &
  sign(coef(fit_tc)[rownames(tt), "X1"] + tt$GID7KO_vs_WT.linear.logFC) ==
  sign(coef(fit_tc)[rownames(tt), "X1"]) &
  sign(coef(fit_tc)[rownames(tt), "X1"] + tt$GID8KO_vs_WT.linear.logFC) ==
  sign(coef(fit_tc)[rownames(tt), "X1"]) &
  sign(coef(fit_tc)[rownames(tt), "X1"] + tt$GID9KO_vs_WT.linear.logFC) ==
  sign(coef(fit_tc)[rownames(tt), "X1"])
opposite_sign_strict <-
  sign(coef(fit_tc)[rownames(tt), "X1"] + tt$GID1KO_vs_WT.linear.logFC) !=
  sign(coef(fit_tc)[rownames(tt), "X1"]) &
  sign(coef(fit_tc)[rownames(tt), "X1"] + tt$GID2KO_vs_WT.linear.logFC) !=
  sign(coef(fit_tc)[rownames(tt), "X1"]) &
  sign(coef(fit_tc)[rownames(tt), "X1"] + tt$GID7KO_vs_WT.linear.logFC) !=
  sign(coef(fit_tc)[rownames(tt), "X1"]) &
  sign(coef(fit_tc)[rownames(tt), "X1"] + tt$GID8KO_vs_WT.linear.logFC) !=
  sign(coef(fit_tc)[rownames(tt), "X1"]) &
  sign(coef(fit_tc)[rownames(tt), "X1"] + tt$GID9KO_vs_WT.linear.logFC) !=
  sign(coef(fit_tc)[rownames(tt), "X1"])


table(same_sign)
table(same_sign_strict)
table(same_sign, same_sign_strict, tt$P.Value < 0.05)
table(opposite_sign, opposite_sign_strict, tt$P.Value < 0.05)




head(tt[same_sign & tt$adj.P.Val < 0.05, ])
head(tt[!same_sign & tt$adj.P.Val < 0.05, ])

pdf(
  here("tmp", "same_sign.linear.pdf"),
  width = 9,
  height = 3)
# TODO: Remove subsetting.
for (g in rownames(tt[same_sign & tt$adj.P.Val < 0.05, ])[1:100]) {
  message(g)
  x1 <- cbind(data.frame(expression = lcpm[g, ]), y$samples)
  x2 <- cbind(
    data.frame(expression = lcpm_fit_tc[g, ]),
    y$samples[, c("cell_line", "timepoint")])
  p <- ggplot(
    data = x1,
    aes(
      x = timepoint,
      y = expression,
      colour = cell_line,
      group = cell_line_rep)) +
    geom_point() +
    geom_line(lty = 2, alpha = 0.5) +
    geom_smooth(
      aes(
        x = timepoint,
        y = expression,
        colour = cell_line,
        group = cell_line),
      data = x2,
      method = "lm",
      formula = y ~ poly(x, 1),
      se = FALSE,
      lwd = 2) +
    facet_grid(~cell_line) +
    scale_colour_manual(values = cell_line_colours) +
    scale_fill_manual(values = group_colours) +
    ggtitle(g) +
    guides(colour = "none", fill = "none") +
    theme_cowplot() +
    panel_border() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  print(p)
}
dev.off()

pdf(
  here("tmp", "opposite_sign.linear.pdf"),
  width = 9,
  height = 3)
# TODO: Remove subsetting.
for (g in rownames(tt[opposite_sign & tt$adj.P.Val < 0.05, ])[1:100]) {
  message(g)
  x1 <- cbind(data.frame(expression = lcpm[g, ]), y$samples)
  x2 <- cbind(
    data.frame(expression = lcpm_fit_tc[g, ]),
    y$samples[, c("cell_line", "timepoint")])
  p <- ggplot(
    data = x1,
    aes(
      x = timepoint,
      y = expression,
      colour = cell_line,
      group = cell_line_rep)) +
    geom_point() +
    geom_line(lty = 2, alpha = 0.5) +
    geom_smooth(
      aes(
        x = timepoint,
        y = expression,
        colour = cell_line,
        group = cell_line),
      data = x2,
      method = "lm",
      formula = y ~ poly(x, 1),
      se = FALSE,
      lwd = 2) +
    facet_grid(~cell_line) +
    scale_colour_manual(values = cell_line_colours) +
    scale_fill_manual(values = group_colours) +
    ggtitle(g) +
    guides(colour = "none", fill = "none") +
    theme_cowplot() +
    panel_border() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  print(p)
}
dev.off()

pdf(
  here("tmp", "same_sign_strict.linear.pdf"),
  width = 9,
  height = 3)
# TODO: Remove subsetting.
for (g in rownames(tt[same_sign_strict & tt$adj.P.Val < 0.05, ])[1:100]) {
  message(g)
  x1 <- cbind(data.frame(expression = lcpm[g, ]), y$samples)
  x2 <- cbind(
    data.frame(expression = lcpm_fit_tc[g, ]),
    y$samples[, c("cell_line", "timepoint")])
  p <- ggplot(
    data = x1,
    aes(
      x = timepoint,
      y = expression,
      colour = cell_line,
      group = cell_line_rep)) +
    geom_point() +
    geom_line(lty = 2, alpha = 0.5) +
    geom_smooth(
      aes(
        x = timepoint,
        y = expression,
        colour = cell_line,
        group = cell_line),
      data = x2,
      method = "lm",
      formula = y ~ poly(x, 1),
      se = FALSE,
      lwd = 2) +
    facet_grid(~cell_line) +
    scale_colour_manual(values = cell_line_colours) +
    scale_fill_manual(values = group_colours) +
    ggtitle(g) +
    guides(colour = "none", fill = "none") +
    theme_cowplot() +
    panel_border() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  print(p)
}
dev.off()

pdf(
  here("tmp", "opposite_sign_strict.linear.pdf"),
  width = 9,
  height = 3)
# TODO: Remove subsetting.
for (g in rownames(tt[opposite_sign_strict & tt$adj.P.Val < 0.05, ])[1:100]) {
  message(g)
  x1 <- cbind(data.frame(expression = lcpm[g, ]), y$samples)
  x2 <- cbind(
    data.frame(expression = lcpm_fit_tc[g, ]),
    y$samples[, c("cell_line", "timepoint")])
  p <- ggplot(
    data = x1,
    aes(
      x = timepoint,
      y = expression,
      colour = cell_line,
      group = cell_line_rep)) +
    geom_point() +
    geom_line(lty = 2, alpha = 0.5) +
    geom_smooth(
      aes(
        x = timepoint,
        y = expression,
        colour = cell_line,
        group = cell_line),
      data = x2,
      method = "lm",
      formula = y ~ poly(x, 1),
      se = FALSE,
      lwd = 2) +
    facet_grid(~cell_line) +
    scale_colour_manual(values = cell_line_colours) +
    scale_fill_manual(values = group_colours) +
    ggtitle(g) +
    guides(colour = "none", fill = "none") +
    theme_cowplot() +
    panel_border() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  print(p)
}
dev.off()

# TODO: Make Just 2-3 PDFs of time course data?
#       1. Containing 'linear' smooth of all genes (regardless of DE status)
#       2. Containing 'quadratic' smooth of all genes (regardless of DE status)
#       3. Containing 'linear' + 'quadratic' smooth of all genes (regardless of
#          DE status); is this even possible using ggplot? Or perhaps that's
#          what I'm already doing with poly(x, 2)?



# OLD CODE -----------



  p <- ggplot(
    data = cbind(data.frame(lcpm = lcpm[g, ]), y$samples),
    mapping = aes(x = timepoint, y = lcpm, group = cell_line_rep)) +
    geom_smooth(
      aes(x = timepoint, y = lcpm, group = cell_line, colour = cell_line),
      se = FALSE,
      lwd = 2) +
    scale_colour_manual(values = cell_line_colours) +
    geom_point(aes(fill = group, colour = cell_line), shape = 21, size = 2) +
    scale_fill_manual(values = group_colours) +
    geom_line(aes(colour = cell_line), alpha = 0.5, lty = 2) +
    facet_grid(~cell_line) +
    ggtitle(g) +
    guides(colour = "none", fill = "none") +
    theme_cowplot() +
    panel_border() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  print(p)
}
dev.off()




tt_gid1 <- topTable(
  fit_tc,
  coef = "GID1KO.X1",
  sort.by = "T")



logCPM.obs <- cpm(y, log = TRUE)
lcpm_fit <- coef(fit) %*% t(design)

tt <- topTable(
  fit,
  coef = "GID9KO:X",
  number = Inf,
  p.value = 0.05)
table(sign(tt$logFC))

g <- rownames(tt[tt$logFC > 0 , ])[1]
x1 <- cbind(
  data.frame(
    expression = logCPM.obs[g, ]),
  y$samples)
x2 <- cbind(
  data.frame(
    expression = lcpm_fit[g, ]),
  y$samples[, c("cell_line", "timepoint", "cell_line_rep")])
ggplot(
  data = x1,
  aes(x = timepoint, y = expression, colour = cell_line, group = cell_line_rep)) +
  geom_point() +
  geom_line(lty = 2, alpha = 0.5) +
  geom_smooth(
    aes(x = timepoint, y = expression, colour = cell_line, group = cell_line_rep),
    data = x2,
    method = "lm",
    formula = y ~ poly(x, 1),
    se = FALSE,
    lwd = 2) +
  facet_grid(~cell_line) +
  scale_colour_manual(values = cell_line_colours) +
  scale_fill_manual(values = group_colours) +
  ggtitle(g) +
  guides(colour = "none", fill = "none") +
  theme_cowplot() +
  panel_border() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# NOTE: This is an ANOVA-like test to identify different (linear) time trends
#       in any one of the KOs vs. WT.
tt <- topTable(
  fit,
  coef = c("GID1KO:X", "GID2KO:X", "GID7KO:X", "GID8KO:X", "GID9KO:X"),
  number = Inf,
  p.value = 0.05)
tt$aveLogFC <- rowMeans(tt[, c("GID1KO.X", "GID2KO.X", "GID7KO.X", "GID8KO.X", "GID9KO.X")])
# NOTE: A sort of pairwise consistency of logFC.
# TODO: A better filter for this (e.g., that takes into pairwise significance)?
tt$keep <- matrixStats::rowAlls(
  tt[, c("GID1KO.X", "GID2KO.X", "GID7KO.X", "GID8KO.X", "GID9KO.X")] > 0) |
  matrixStats::rowAlls(
    tt[, c("GID1KO.X", "GID2KO.X", "GID7KO.X", "GID8KO.X", "GID9KO.X")] < 0)
g <- rownames(tt[tt$aveLogFC > 0 , ])[1]
x1 <- cbind(
  data.frame(
    expression = logCPM.obs[g, ]),
  y$samples)
x2 <- cbind(
  data.frame(
    expression = lcpm_fit[g, ]),
  y$samples[, c("cell_line", "timepoint", "cell_line_rep")])
ggplot(
  data = x1,
  aes(x = timepoint, y = expression, colour = cell_line, group = cell_line_rep)) +
  geom_point() +
  geom_line(lty = 2, alpha = 0.5) +
  geom_smooth(
    aes(x = timepoint, y = expression, colour = cell_line, group = cell_line_rep),
    data = x2,
    method = "lm",
    formula = y ~ poly(x, 1),
    se = FALSE,
    lwd = 2) +
  facet_grid(~cell_line) +
  scale_colour_manual(values = cell_line_colours) +
  scale_fill_manual(values = group_colours) +
  ggtitle(g) +
  guides(colour = "none", fill = "none") +
  theme_cowplot() +
  panel_border() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# TODO: Find genes with opposite trends over time in WT and KOs amongst the DEGs.
# NOTE: This is a stronger requirement than the above ANOVA-like test because
#       in the above a gene could have a negative (resp. positive) trend over
#       time in all genotypes but it is more strongly negative (resp. positive)
#       in WT than any of the KOs.
diffs <- apply(
  coef(fit)[, c("GID1KO:X", "GID2KO:X", "GID7KO:X", "GID8KO:X", "GID9KO:X")],
  2,
  `+`,
  coef(fit)[, "X", drop = FALSE])
tt$opposite_trend <- rownames(tt) %in%
  names(
    which(
      sign(coef(fit)[, "X"]) != sign(diffs[, 1]) &
        (matrixStats::rowAlls(sign(diffs) > 0) |
           matrixStats::rowAlls(sign(diffs) < 0))))
# TODO: Is there a way to achieve the above with contrasts?
table(tt$opposite_trend)
g <- rownames(tt[tt$opposite_trend, ])[2]
x1 <- cbind(
  data.frame(
    expression = logCPM.obs[g, ]),
  y$samples)
x2 <- cbind(
  data.frame(
    expression = lcpm_fit[g, ]),
  y$samples[, c("cell_line", "timepoint", "cell_line_rep")])
ggplot(
  data = x1,
  aes(x = timepoint, y = expression, colour = cell_line, group = cell_line_rep)) +
  geom_point() +
  geom_line(lty = 2, alpha = 0.5) +
  geom_smooth(
    aes(x = timepoint, y = expression, colour = cell_line, group = cell_line_rep),
    data = x2,
    method = "lm",
    formula = y ~ poly(x, 1),
    se = FALSE,
    lwd = 2) +
  facet_grid(~cell_line) +
  scale_colour_manual(values = cell_line_colours) +
  scale_fill_manual(values = group_colours) +
  ggtitle(g) +
  guides(colour = "none", fill = "none") +
  theme_cowplot() +
  panel_border() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
