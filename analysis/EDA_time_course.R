# NOTE: Starting after filtering out low-quality replicates and summing
#       technical replicates

library(splines)
X <- ns(as.integer(y$samples$timepoint), df = 3)
Group <- relevel(y$samples$cell_line, "WT")

design <- model.matrix(~Group*X)
colnames(design) <- sub("Group", "", colnames(design))

fit0 <- voomLmFit(
  y,
  design,
  block = y$samples$cell_line_rep,
  sample.weights = TRUE,
  plot = TRUE)
fit0 <- eBayes(fit0)

tt <- topTable(
  fit0,
  coef = c("GID1KO:X1", "GID1KO:X2", "GID1KO:X3"),
  number = Inf,
  p.value = 0.05)
nrow(tt)

# This does a t-test rather than F-test, so no good
# fit <- contrasts.fit(
#   fit0,
#   coef = match(c("GID1KO:X1", "GID1KO:X2", "GID1KO:X3"), colnames(fit0)))
# fit <- eBayes(fit)
# summary(decideTests(fit))

# Plots ------------------------------------------------------------------------

lcpm <- cpm(y, log = TRUE)

g <- rownames(tt)[1]
ggplot(
  data = cbind(data.frame(lcpm = lcpm[g, ]), y$samples),
  mapping = aes(x = timepoint, y = lcpm, colour = cell_line, group = cell_line_rep)) +
  geom_point() +
  geom_line(alpha = 0.5, lty = 2) +
  # TODO: Is there a way to use the spline fit?
  geom_smooth(
    aes(x = timepoint, y = lcpm, group = cell_line),
    se = FALSE,
    lwd = 2) +
  facet_grid(~cell_line) +
  ggtitle(g) +
  guides(colour = "none") +
  theme_cowplot() +
  panel_border() +
  scale_colour_manual(values = cell_line_colours)
ggplot(
  data = cbind(data.frame(lcpm = lcpm[g, ]), y$samples),
  mapping = aes(x = timepoint, y = lcpm, group = cell_line_rep)) +
  # TODO: Is there a way to use the spline fit?
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
  panel_border()

# TODO: There are 9 possible patterns of the signs for a spline fit with 3
#       terms (-,-,-), (-,-,+), ..., (+,+,+); see if there's any pattern to the
#       genes in each category.

# NOTE: Corresponds to (+,+,+)
i1 <- which(rowMins(as.matrix(tt[, c("GID1KO.X1", "GID1KO.X2", "GID1KO.X3")])) > 0)
g <- row.names(tt)[i1[1]]
ggplot(
  data = cbind(data.frame(lcpm = lcpm[g, ]), y$samples),
  mapping = aes(x = timepoint, y = lcpm, colour = cell_line, group = cell_line_rep)) +
  geom_point() +
  geom_line(alpha = 0.5, lty = 2) +
  # TODO: Is there a way to use the spline fit?
  geom_smooth(
    aes(x = timepoint, y = lcpm, group = cell_line),
    se = FALSE,
    lwd = 2) +
  facet_grid(~cell_line) +
  ggtitle(g) +
  guides(colour = "none") +
  theme_cowplot() +
  panel_border() +
  scale_colour_manual(values = cell_line_colours)
ggplot(
  data = cbind(data.frame(lcpm = lcpm[g, ]), y$samples),
  mapping = aes(x = timepoint, y = lcpm, group = cell_line_rep)) +
  # TODO: Is there a way to use the spline fit?
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
  panel_border()

# NOTE: Corresponds to (-,-,-)
i2 <- which(rowMaxs(as.matrix(tt[, c("GID1KO.X1", "GID1KO.X2", "GID1KO.X3")])) < 0)
g <- row.names(tt)[i2[1]]
ggplot(
  data = cbind(data.frame(lcpm = lcpm[g, ]), y$samples),
  mapping = aes(x = timepoint, y = lcpm, colour = cell_line, group = cell_line_rep)) +
  geom_point() +
  geom_line(alpha = 0.5, lty = 2) +
  # TODO: Is there a way to use the spline fit?
  geom_smooth(
    aes(x = timepoint, y = lcpm, group = cell_line),
    se = FALSE,
    lwd = 2) +
  facet_grid(~cell_line) +
  ggtitle(g) +
  guides(colour = "none") +
  theme_cowplot() +
  panel_border() +
  scale_colour_manual(values = cell_line_colours)
ggplot(
  data = cbind(data.frame(lcpm = lcpm[g, ]), y$samples),
  mapping = aes(x = timepoint, y = lcpm, group = cell_line_rep)) +
  # TODO: Is there a way to use the spline fit?
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
  panel_border()

# A non-DE gene
g <- tail(
  rownames(
    topTable(
      fit0,
      coef = c("GID1KO:X1", "GID1KO:X2", "GID1KO:X3"),
      number = Inf,
      p.value = 1)),
  1)
ggplot(
  data = cbind(data.frame(lcpm = lcpm[g, ]), y$samples),
  mapping = aes(x = timepoint, y = lcpm, colour = cell_line, group = cell_line_rep)) +
  geom_point() +
  geom_line(alpha = 0.5, lty = 2) +
  # TODO: Is there a way to use the spline fit?
  geom_smooth(
    aes(x = timepoint, y = lcpm, group = cell_line),
    se = FALSE,
    lwd = 2) +
  facet_grid(~cell_line) +
  ggtitle(g) +
  guides(colour = "none") +
  theme_cowplot() +
  panel_border() +
  scale_colour_manual(values = cell_line_colours)
ggplot(
  data = cbind(data.frame(lcpm = lcpm[g, ]), y$samples),
  mapping = aes(x = timepoint, y = lcpm, group = cell_line_rep)) +
  # TODO: Is there a way to use the spline fit?
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
  panel_border()

# TODO: Could also extract genes that change over time within a given cell line
#       (rather than those that change where there change is different to the
#       change within WT).

# Heatmaps ---------------------------------------------------------------------

# TODO: The signs of the spline terms don't seem very useful/interpretable.
z <- sign(tt[, c("GID1KO.X1", "GID1KO.X2", "GID1KO.X3")])
zz <- apply(z, 1, function(x) paste0(x, collapse = "."))
pheatmap::pheatmap(
  lcpm[head(rownames(tt), 100), order(y$samples$group)],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = y8$samples[, c("timepoint", "cell_line")],
  # TODO: Main
  # main = gsub("_vs_", " vs. ", j),
  annotation_colors = list(
    cell_line = cell_line_colours,
    timepoint = timepoint_colour,
    group = group_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  annotation_row = data.frame(pattern = zz))

lcpmbg <- cpmByGroup(y, group = y$samples$group, log = TRUE)
pheatmap::pheatmap(
  lcpmbg[rownames(tt)[order(zz)], ],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    group = colnames(lcpmbg),
    row.names = colnames(lcpmbg)),
  # TODO: Main
  # main = gsub("_vs_", " vs. ", j),
  annotation_colors = list(group = group_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_rownames = FALSE,
  annotation_row = data.frame(pattern = zz))
pheatmap::pheatmap(
  lcpmbg[rownames(tt)[order(zz)], ],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    group = colnames(lcpmbg),
    row.names = colnames(lcpmbg)),
  # TODO: Main
  # main = gsub("_vs_", " vs. ", j),
  annotation_colors = list(group = group_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  show_rownames = FALSE,
  annotation_row = data.frame(pattern = zz))

# Just showing the groups involved in the contrast
pheatmap::pheatmap(
  lcpmbg[rownames(tt)[order(zz)], grep("WT|GID1KO", colnames(lcpmbg))],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    group = colnames(lcpmbg),
    row.names = colnames(lcpmbg))[
      grep("WT|GID1KO", colnames(lcpmbg), value = TRUE), , drop = FALSE],
  # TODO: Main
  # main = gsub("_vs_", " vs. ", j),
  annotation_colors = list(group = group_colours[
    grep("WT|GID1KO", colnames(lcpmbg), value = TRUE)]),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_rownames = FALSE,
  annotation_row = data.frame(pattern = zz))
pheatmap::pheatmap(
  lcpmbg[rownames(tt)[order(zz)], grep("WT|GID1KO", colnames(lcpmbg))],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    group = colnames(lcpmbg),
    row.names = colnames(lcpmbg))[
      grep("WT|GID1KO", colnames(lcpmbg), value = TRUE), , drop = FALSE],
  # TODO: Main
  # main = gsub("_vs_", " vs. ", j),
  annotation_colors = list(group = group_colours[
    grep("WT|GID1KO", colnames(lcpmbg), value = TRUE)]),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  show_rownames = FALSE,
  annotation_row = data.frame(pattern = zz))

# TODO: Could cluster genes with a significant time effect to group them into
#       common patterns.

# df = 3 -----------------------------------------------------------------------

X_3 <- ns(as.integer(y$samples$timepoint), df = 3)
Group <- relevel(y$samples$cell_line, "WT")
design_tc_3 <- model.matrix(~Group * X_3)
colnames(design_tc_3) <- sub("Group", "", colnames(design_tc_3))

fit_tc_3 <- voomLmFit(
  y,
  design_tc_3,
  block = y$samples$cell_line_rep,
  sample.weights = TRUE,
  plot = TRUE)
fit_tc_3 <- eBayes(fit_tc_3)

# Summarise the number of DEGs in each comparison
tt_gid1_3 <- topTable(
  fit_tc_3,
  coef = c("GID1KO:X_31", "GID1KO:X_32", "GID1KO:X_33"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
tt_gid2_3  <- topTable(
  fit_tc_3,
  coef = c("GID2KO:X_31", "GID2KO:X_32", "GID2KO:X_33"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
tt_gid7_3  <- topTable(
  fit_tc_3,
  coef = c("GID7KO:X_31", "GID7KO:X_32", "GID7KO:X_33"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
tt_gid8_3  <- topTable(
  fit_tc_3,
  coef = c("GID8KO:X_31", "GID8KO:X_32", "GID8KO:X_33"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
tt_gid9_3  <- topTable(
  fit_tc_3,
  coef = c("GID9KO:X_31", "GID9KO:X_32", "GID9KO:X_33"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
# NOTE: This is an ANOVA-like test.
tt_any_ko_3 <-  topTable(
  fit_tc_3,
  coef = c(
    "GID1KO:X_31", "GID1KO:X_32", "GID1KO:X_33",
    "GID2KO:X_31", "GID2KO:X_32", "GID2KO:X_33",
    "GID7KO:X_31", "GID7KO:X_32", "GID7KO:X_33",
    "GID8KO:X_31", "GID8KO:X_32", "GID8KO:X_33",
    "GID9KO:X_31", "GID9KO:X_32", "GID9KO:X_33"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
sdt_3 <- data.frame(
  Sig = c(
    sum(tt_gid1_3$adj.P.Val < 0.05),
    sum(tt_gid2_3$adj.P.Val < 0.05),
    sum(tt_gid7_3$adj.P.Val < 0.05),
    sum(tt_gid8_3$adj.P.Val < 0.05),
    sum(tt_gid9_3$adj.P.Val < 0.05),
    sum(tt_any_ko_3$adj.P.Val < 0.05)),
  NotSig = c(
    sum(tt_gid1_3$adj.P.Val >= 0.05),
    sum(tt_gid2_3$adj.P.Val >= 0.05),
    sum(tt_gid7_3$adj.P.Val >= 0.05),
    sum(tt_gid8_3$adj.P.Val >= 0.05),
    sum(tt_gid9_3$adj.P.Val >= 0.05),
    sum(tt_any_ko_3$adj.P.Val >= 0.05)),
  row.names = c(
    "GID1KO_vs_WT", "GID2KO_vs_WT", "GID7KO_vs_WT",
    "GID8KO_vs_WT", "GID9KO_vs_WT", "any_KO_vs_WT"))

g <- rownames(tt_any_ko_3)[1]
p3 <- ggplot(
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

# df = 2 -----------------------------------------------------------------------

X_2 <- ns(as.integer(y$samples$timepoint), df = 2)
Group <- relevel(y$samples$cell_line, "WT")
design_tc_2 <- model.matrix(~Group * X_2)
colnames(design_tc_2) <- sub("Group", "", colnames(design_tc_2))

fit_tc_2 <- voomLmFit(
  y,
  design_tc_2,
  block = y$samples$cell_line_rep,
  sample.weights = TRUE,
  plot = TRUE)
fit_tc_2 <- eBayes(fit_tc_2)

# Summarise the number of DEGs in each comparison
tt_gid1_2 <- topTable(
  fit_tc_2,
  coef = c("GID1KO:X_21", "GID1KO:X_22"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
tt_gid2_2  <- topTable(
  fit_tc_2,
  coef = c("GID2KO:X_21", "GID2KO:X_22"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
tt_gid7_2  <- topTable(
  fit_tc_2,
  coef = c("GID7KO:X_21", "GID7KO:X_22"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
tt_gid8_2  <- topTable(
  fit_tc_2,
  coef = c("GID8KO:X_21", "GID8KO:X_22"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
tt_gid9_2  <- topTable(
  fit_tc_2,
  coef = c("GID9KO:X_21", "GID9KO:X_22"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
# NOTE: This is an ANOVA-like test.
tt_any_ko_2 <-  topTable(
  fit_tc_2,
  coef = c(
    "GID1KO:X_21", "GID1KO:X_22",
    "GID2KO:X_21", "GID2KO:X_22",
    "GID7KO:X_21", "GID7KO:X_22",
    "GID8KO:X_21", "GID8KO:X_22",
    "GID9KO:X_21", "GID9KO:X_22"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
sdt_2 <- data.frame(
  Sig = c(
    sum(tt_gid1_2$adj.P.Val < 0.05),
    sum(tt_gid2_2$adj.P.Val < 0.05),
    sum(tt_gid7_2$adj.P.Val < 0.05),
    sum(tt_gid8_2$adj.P.Val < 0.05),
    sum(tt_gid9_2$adj.P.Val < 0.05),
    sum(tt_any_ko_2$adj.P.Val < 0.05)),
  NotSig = c(
    sum(tt_gid1_2$adj.P.Val >= 0.05),
    sum(tt_gid2_2$adj.P.Val >= 0.05),
    sum(tt_gid7_2$adj.P.Val >= 0.05),
    sum(tt_gid8_2$adj.P.Val >= 0.05),
    sum(tt_gid9_2$adj.P.Val >= 0.05),
    sum(tt_any_ko_2$adj.P.Val >= 0.05)),
  row.names = c(
    "GID1KO_vs_WT", "GID2KO_vs_WT", "GID7KO_vs_WT",
    "GID8KO_vs_WT", "GID9KO_vs_WT", "any_KO_vs_WT"))

p2 <- ggplot(
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

# df = 1 -----------------------------------------------------------------------

X_1 <- ns(as.integer(y$samples$timepoint), df = 1)
Group <- relevel(y$samples$cell_line, "WT")
design_tc_1 <- model.matrix(~Group * X_1)
colnames(design_tc_1) <- sub("Group", "", colnames(design_tc_1))

fit_tc_1 <- voomLmFit(
  y,
  design_tc_1,
  block = y$samples$cell_line_rep,
  sample.weights = TRUE,
  plot = TRUE)
fit_tc_1 <- eBayes(fit_tc_1)

# Summarise the number of DEGs in each comparison
tt_gid1_1 <- topTable(
  fit_tc_1,
  coef = c("GID1KO:X_1"),
  number = Inf,
  p.value = 1,
  sort.by = "t")
tt_gid2_1  <- topTable(
  fit_tc_1,
  coef = c("GID2KO:X_1"),
  number = Inf,
  p.value = 1,
  sort.by = "t")
tt_gid7_1  <- topTable(
  fit_tc_1,
  coef = c("GID7KO:X_1"),
  number = Inf,
  p.value = 1,
  sort.by = "t")
tt_gid8_1  <- topTable(
  fit_tc_1,
  coef = c("GID8KO:X_1"),
  number = Inf,
  p.value = 1,
  sort.by = "t")
tt_gid9_1  <- topTable(
  fit_tc_1,
  coef = c("GID9KO:X_1"),
  number = Inf,
  p.value = 1,
  sort.by = "t")
# NOTE: This is an ANOVA-like test.
tt_any_ko_1 <-  topTable(
  fit_tc_1,
  coef = c(
    "GID1KO:X_1",
    "GID2KO:X_1",
    "GID7KO:X_1",
    "GID8KO:X_1",
    "GID9KO:X_1"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
sdt_1 <- data.frame(
  Sig = c(
    sum(tt_gid1_1$adj.P.Val < 0.05),
    sum(tt_gid2_1$adj.P.Val < 0.05),
    sum(tt_gid7_1$adj.P.Val < 0.05),
    sum(tt_gid8_1$adj.P.Val < 0.05),
    sum(tt_gid9_1$adj.P.Val < 0.05),
    sum(tt_any_ko_1$adj.P.Val < 0.05)),
  NotSig = c(
    sum(tt_gid1_1$adj.P.Val >= 0.05),
    sum(tt_gid2_1$adj.P.Val >= 0.05),
    sum(tt_gid7_1$adj.P.Val >= 0.05),
    sum(tt_gid8_1$adj.P.Val >= 0.05),
    sum(tt_gid9_1$adj.P.Val >= 0.05),
    sum(tt_any_ko_1$adj.P.Val >= 0.05)),
  row.names = c(
    "GID1KO_vs_WT", "GID2KO_vs_WT", "GID7KO_vs_WT",
    "GID8KO_vs_WT", "GID9KO_vs_WT", "any_KO_vs_WT"))

p1 <- ggplot(
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

# Quadratic fit ----------------------------------------------------------------

XX_2 <- poly(as.integer(y$samples$timepoint), degree = 2)
Group <- relevel(y$samples$cell_line, "WT")
design_tc_22 <- model.matrix(~Group * XX_2)
colnames(design_tc_22) <- sub("Group", "", colnames(design_tc_22))

fit_tc_22 <- voomLmFit(
  y,
  design_tc_22,
  block = y$samples$cell_line_rep,
  sample.weights = TRUE,
  plot = TRUE)
fit_tc_22 <- eBayes(fit_tc_22)

# Summarise the number of DEGs in each comparison
tt_gid1_22 <- topTable(
  fit_tc_22,
  coef = c("GID1KO:XX_21", "GID1KO:XX_22"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
tt_gid2_22  <- topTable(
  fit_tc_22,
  coef = c("GID2KO:XX_21", "GID2KO:XX_22"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
tt_gid7_22  <- topTable(
  fit_tc_22,
  coef = c("GID7KO:XX_21", "GID7KO:XX_22"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
tt_gid8_22  <- topTable(
  fit_tc_22,
  coef = c("GID8KO:XX_21", "GID8KO:XX_22"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
tt_gid9_22  <- topTable(
  fit_tc_22,
  coef = c("GID9KO:XX_21", "GID9KO:XX_22"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
# NOTE: This is an ANOVA-like test.
tt_any_ko_22 <-  topTable(
  fit_tc_22,
  coef = c(
    "GID1KO:XX_21", "GID1KO:XX_22",
    "GID2KO:XX_21", "GID2KO:XX_22",
    "GID7KO:XX_21", "GID7KO:XX_22",
    "GID8KO:XX_21", "GID8KO:XX_22",
    "GID9KO:XX_21", "GID9KO:XX_22"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
sdt_22 <- data.frame(
  Sig = c(
    sum(tt_gid1_22$adj.P.Val < 0.05),
    sum(tt_gid2_22$adj.P.Val < 0.05),
    sum(tt_gid7_22$adj.P.Val < 0.05),
    sum(tt_gid8_22$adj.P.Val < 0.05),
    sum(tt_gid9_22$adj.P.Val < 0.05),
    sum(tt_any_ko_22$adj.P.Val < 0.05)),
  NotSig = c(
    sum(tt_gid1_22$adj.P.Val >= 0.05),
    sum(tt_gid2_22$adj.P.Val >= 0.05),
    sum(tt_gid7_22$adj.P.Val >= 0.05),
    sum(tt_gid8_22$adj.P.Val >= 0.05),
    sum(tt_gid9_22$adj.P.Val >= 0.05),
    sum(tt_any_ko_22$adj.P.Val >= 0.05)),
  row.names = c(
    "GID1KO_vs_WT", "GID2KO_vs_WT", "GID7KO_vs_WT",
    "GID8KO_vs_WT", "GID9KO_vs_WT", "any_KO_vs_WT"))

tt <- topTable(fit_tc_22, coef = c("GID9KO:XX_22"), p.value = 0.05, number = Inf)
g <- rownames(tt[tt$logFC < 0, ]) # ∩ shaped logFC
ggplot(
  data = cbind(data.frame(lcpm = lcpm[g[1], ]), y$samples),
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
  ggtitle(g[1]) +
  guides(colour = "none", fill = "none") +
  theme_cowplot() +
  panel_border() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

coef(fit_tc_22)[g[1], c("GID1KO")]

# 1. Samples ordered by `group`.
hm <- pheatmap(
  lcpm[
    g,
    y$samples$cell_line %in% c("GID1KO", "WT")][
      ,
      order(y$samples$group[y$samples$cell_line %in% c("GID1KO", "WT")])],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = y$samples[, "group", drop = FALSE],
  annotation_colors = list(group = group_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  width = 12,
  height = 12,
  cutree_rows = 3,
  show_rownames = TRUE,
  silent = FALSE)
gc <- cutree(hm$tree_row, k = 3)

p <- lapply(g[gc == 1][1:9], function(gg) {
  ggplot(
    data = cbind(data.frame(lcpm = lcpm[gg, ]), y$samples),
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
    ggtitle(gg) +
    guides(colour = "none", fill = "none") +
    theme_cowplot() +
    panel_border() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
})
wrap_plots(p, ncol = 3) + plot_annotation("gc = 1")

# Linear fit -------------------------------------------------------------------

XX_1 <- poly(as.integer(y$samples$timepoint), degree = 1)
Group <- relevel(y$samples$cell_line, "WT")
design_tc_12 <- model.matrix(~Group * XX_1)
colnames(design_tc_12) <- sub("Group", "", colnames(design_tc_12))

fit_tc_12 <- voomLmFit(
  y,
  design_tc_12,
  block = y$samples$cell_line_rep,
  sample.weights = TRUE,
  plot = TRUE)
fit_tc_12 <- eBayes(fit_tc_12)

# Summarise the number of DEGs in each comparison
tt_gid1_12 <- topTable(
  fit_tc_12,
  coef = c("GID1KO:XX_1"),
  number = Inf,
  p.value = 1,
  sort.by = "t")
tt_gid2_12  <- topTable(
  fit_tc_12,
  coef = c("GID2KO:XX_1"),
  number = Inf,
  p.value = 1,
  sort.by = "t")
tt_gid7_12  <- topTable(
  fit_tc_12,
  coef = c("GID7KO:XX_1"),
  number = Inf,
  p.value = 1,
  sort.by = "t")
tt_gid8_12  <- topTable(
  fit_tc_12,
  coef = c("GID8KO:XX_1"),
  number = Inf,
  p.value = 1,
  sort.by = "t")
tt_gid9_12  <- topTable(
  fit_tc_12,
  coef = c("GID9KO:XX_1"),
  number = Inf,
  p.value = 1,
  sort.by = "t")
# NOTE: This is an ANOVA-like test.
tt_any_ko_12 <-  topTable(
  fit_tc_12,
  coef = c(
    "GID1KO:XX_1",
    "GID2KO:XX_1",
    "GID7KO:XX_1",
    "GID8KO:XX_1",
    "GID9KO:XX_1"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
sdt_12 <- data.frame(
  Sig = c(
    sum(tt_gid1_12$adj.P.Val < 0.05),
    sum(tt_gid2_12$adj.P.Val < 0.05),
    sum(tt_gid7_12$adj.P.Val < 0.05),
    sum(tt_gid8_12$adj.P.Val < 0.05),
    sum(tt_gid9_12$adj.P.Val < 0.05),
    sum(tt_any_ko_12$adj.P.Val < 0.05)),
  NotSig = c(
    sum(tt_gid1_12$adj.P.Val >= 0.05),
    sum(tt_gid2_12$adj.P.Val >= 0.05),
    sum(tt_gid7_12$adj.P.Val >= 0.05),
    sum(tt_gid8_12$adj.P.Val >= 0.05),
    sum(tt_gid9_12$adj.P.Val >= 0.05),
    sum(tt_any_ko_12$adj.P.Val >= 0.05)),
  row.names = c(
    "GID1KO_vs_WT", "GID2KO_vs_WT", "GID7KO_vs_WT",
    "GID8KO_vs_WT", "GID9KO_vs_WT", "any_KO_vs_WT"))

tt <- topTable(fit_tc_12, coef = c("GID9KO:XX_1"), p.value = 0.05, number = Inf)
g <- rownames(tt[tt$logFC > 0, ]) # positive slope logFC
ggplot(
  data = cbind(data.frame(lcpm = lcpm[g[1], ]), y$samples),
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
  ggtitle(g[1]) +
  guides(colour = "none", fill = "none") +
  theme_cowplot() +
  panel_border() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# 1. Samples ordered by `group`.
hm <- pheatmap(
  lcpm[
    g,
    y$samples$cell_line %in% c("GID1KO", "WT")][
      ,
      order(y$samples$group[y$samples$cell_line %in% c("GID1KO", "WT")])],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = y$samples[, "group", drop = FALSE],
  annotation_colors = list(group = group_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  width = 12,
  height = 12,
  cutree_rows = 4,
  show_rownames = TRUE,
  silent = FALSE)
gc <- cutree(hm$tree_row, k = 4)

p <- lapply(g[gc == 1][1:9], function(gg) {
  ggplot(
    data = cbind(data.frame(lcpm = lcpm[gg, ]), y$samples),
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
    ggtitle(gg) +
    guides(colour = "none", fill = "none") +
    theme_cowplot() +
    panel_border() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
})
wrap_plots(p, ncol = 3) + plot_annotation("gc = 1")

p <- lapply(g[gc == 2][1:9], function(gg) {
  ggplot(
    data = cbind(data.frame(lcpm = lcpm[gg, ]), y$samples),
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
    ggtitle(gg) +
    guides(colour = "none", fill = "none") +
    theme_cowplot() +
    panel_border() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
})
wrap_plots(p, ncol = 3) + plot_annotation("gc = 2")

p <- lapply(g[gc == 3][1:9], function(gg) {
  ggplot(
    data = cbind(data.frame(lcpm = lcpm[gg, ]), y$samples),
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
    ggtitle(gg) +
    guides(colour = "none", fill = "none") +
    theme_cowplot() +
    panel_border() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
})
wrap_plots(p, ncol = 3) + plot_annotation("gc = 3")

p <- lapply(g[gc == 4][1:9], function(gg) {
  ggplot(
    data = cbind(data.frame(lcpm = lcpm[gg, ]), y$samples),
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
    ggtitle(gg) +
    guides(colour = "none", fill = "none") +
    theme_cowplot() +
    panel_border() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
})
wrap_plots(p, ncol = 3) + plot_annotation("gc = 4")


# TODO: UP TO HERE. Take a look at the above plots to understand better what
#       we're looking at.
# TODO: Can we do it using the ANOVA-style test rather than each KO vs. WT?

# 'Simple' interaction model ---------------------------------------------------

# TODO: Gordon suggested just looking at genes with a significant
#       genotype:timepoint interaction, which is simpler than the spline/poly
#       fits but as far as I can tell doesn't give a simple way to prioritise
#       genes by the pattern-over-time of their expression or logFC.
