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

# This do t-test rather than F-test, so no good
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
