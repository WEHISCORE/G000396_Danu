# EDA plotting logFCs rather than logCPM in heatmaps.
# Peter Hickey
# 2024-02-16

# Lasender sex genes (consecutive timepoints) ----------------------------------

# NOTE: Using logFCs from consecutive timepoints, i.e.
#       Day_6_vs_Day3, Day_9_vs_Day_6, Day_12_vs_Day_9

# TODO: A more robust way of pulling out the relevant contrasts.
j <- colnames(cfit)[1:18]
pheatmap(
  coef(cfit)[lasonder_genes, j],
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(
    -max(abs(coef(cfit)[lasonder_genes, j])),
    max(abs(coef(cfit)[lasonder_genes, j])),
    length.out = 101),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    comparison = factor(
      c(rep("6_vs_3", 6), rep("9_vs_6", 6), rep("12_vs_9", 6)),
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    cell_line = factor(
      sapply(strsplit(j, "\\."), "[[", 1),
      levels(y$samples$cell_line)),
    row.names = j),
  annotation_row = data.frame(
    sex = c(
      rep("Male", length(intersect(male_genes$`Male genes`, rownames(lcpm)))),
      rep(
        "Female",
        length(intersect(female_genes$`Female genes`, rownames(lcpm))))),
    row.names = lasonder_genes),
  main = "Lasonder sex genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    comparison = setNames(
      timepoint_colours[c("Day_6", "Day_9", "Day_12")],
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    sex = c(Male = "skyblue3", Female = "deeppink")),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE)

pheatmap(
  coef(cfit)[lasonder_genes, j],
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(
    -max(abs(coef(cfit)[lasonder_genes, j])),
    max(abs(coef(cfit)[lasonder_genes, j])),
    length.out = 101),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    comparison = factor(
      c(rep("6_vs_3", 6), rep("9_vs_6", 6), rep("12_vs_9", 6)),
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    cell_line = factor(
      sapply(strsplit(j, "\\."), "[[", 1),
      levels(y$samples$cell_line)),
    row.names = j),
  annotation_row = data.frame(
    sex = c(
      rep("Male", length(intersect(male_genes$`Male genes`, rownames(lcpm)))),
      rep(
        "Female",
        length(intersect(female_genes$`Female genes`, rownames(lcpm))))),
    row.names = lasonder_genes),
  main = "Lasonder sex genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    comparison = setNames(
      timepoint_colours[c("Day_6", "Day_9", "Day_12")],
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    sex = c(Male = "skyblue3", Female = "deeppink")),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = TRUE)

# Lasender sex genes (relative to Day_3) ---------------------------------------

# NOTE: Using logFCs from comparison to Day_3, i.e.
#       Day_6_vs_Day3, Day_9_vs_Day_3, Day_12_vs_Day_3

cm2 <- makeContrasts(
  # Day6 vs. Day3 (within cell line)
  GID1KO.Day_6_vs_GID1KO.Day_3 = GID1KO.Day_6 - GID1KO.Day_3,
  GID2KO.Day_6_vs_GID2KO.Day_3 = GID2KO.Day_6 - GID2KO.Day_3,
  GID7KO.Day_6_vs_GID7KO.Day_3 = GID7KO.Day_6 - GID7KO.Day_3,
  GID8KO.Day_6_vs_GID8KO.Day_3 = GID8KO.Day_6 - GID8KO.Day_3,
  GID9KO.Day_6_vs_GID9KO.Day_3 = GID9KO.Day_6 - GID9KO.Day_3,
  WT.Day_6_vs_WT.Day_3 = WT.Day_6 - WT.Day_3,

  # Day9 vs. Day3 (within cell line)
  GID1KO.Day_9_vs_GID1KO.Day_3 = GID1KO.Day_9 - GID1KO.Day_3,
  GID2KO.Day_9_vs_GID2KO.Day_3 = GID2KO.Day_9 - GID2KO.Day_3,
  GID7KO.Day_9_vs_GID7KO.Day_3 = GID7KO.Day_9 - GID7KO.Day_3,
  GID8KO.Day_9_vs_GID8KO.Day_3 = GID8KO.Day_9 - GID8KO.Day_3,
  GID9KO.Day_9_vs_GID9KO.Day_3 = GID9KO.Day_9 - GID9KO.Day_3,
  WT.Day_9_vs_WT.Day_3 = WT.Day_9 - WT.Day_3,

  # Day12 vs. Day3 (within cell line)
  GID1KO.Day_12_vs_GID1KO.Day_3 = GID1KO.Day_12 - GID1KO.Day_3,
  GID2KO.Day_12_vs_GID2KO.Day_3 = GID2KO.Day_12 - GID2KO.Day_3,
  GID7KO.Day_12_vs_GID7KO.Day_3 = GID7KO.Day_12 - GID7KO.Day_3,
  GID8KO.Day_12_vs_GID8KO.Day_3 = GID8KO.Day_12 - GID8KO.Day_3,
  GID9KO.Day_12_vs_GID9KO.Day_3 = GID9KO.Day_12 - GID9KO.Day_3,
  WT.Day_12_vs_WT.Day_3 = WT.Day_12 - WT.Day_3,

  levels = design)
cfit2 <- contrasts.fit(fit, cm2)

# TODO: A more robust way of pulling out the relevant contrasts.
j2 <- colnames(cfit2)[1:18]
pheatmap(
  coef(cfit2)[lasonder_genes, j2],
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(
    -max(abs(coef(cfit2)[lasonder_genes, j2])),
    max(abs(coef(cfit2)[lasonder_genes, j2])),
    length.out = 101),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    comparison = factor(
      c(rep("6_vs_3", 6), rep("9_vs_3", 6), rep("12_vs_3", 6)),
      c("6_vs_3", "9_vs_3", "12_vs_3")),
    cell_line = factor(
      sapply(strsplit(j2, "\\."), "[[", 1),
      levels(y$samples$cell_line)),
    row.names = j2),
  annotation_row = data.frame(
    sex = c(
      rep("Male", length(intersect(male_genes$`Male genes`, rownames(lcpm)))),
      rep(
        "Female",
        length(intersect(female_genes$`Female genes`, rownames(lcpm))))),
    row.names = lasonder_genes),
  main = "Lasonder sex genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    comparison = setNames(
      timepoint_colours[c("Day_6", "Day_9", "Day_12")],
      c("6_vs_3", "9_vs_3", "12_vs_3")),
    sex = c(Male = "skyblue3", Female = "deeppink")),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE)

pheatmap(
  coef(cfit2)[lasonder_genes, j2],
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(
    -max(abs(coef(cfit2)[lasonder_genes, j2])),
    max(abs(coef(cfit2)[lasonder_genes, j2])),
    length.out = 101),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    comparison = factor(
      c(rep("6_vs_3", 6), rep("9_vs_3", 6), rep("12_vs_3", 6)),
      c("6_vs_3", "9_vs_3", "12_vs_3")),
    cell_line = factor(
      sapply(strsplit(j2, "\\."), "[[", 1),
      levels(y$samples$cell_line)),
    row.names = j2),
  annotation_row = data.frame(
    sex = c(
      rep("Male", length(intersect(male_genes$`Male genes`, rownames(lcpm)))),
      rep(
        "Female",
        length(intersect(female_genes$`Female genes`, rownames(lcpm))))),
    row.names = lasonder_genes),
  main = "Lasonder sex genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    comparison = setNames(
      timepoint_colours[c("Day_6", "Day_9", "Day_12")],
      c("6_vs_3", "9_vs_3", "12_vs_3")),
    sex = c(Male = "skyblue3", Female = "deeppink")),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = TRUE)

# Timecourse genes ---------------------------------------------------------=---

tt <- topTable(
  fit_tc,
  coef = c(
    "GID1KO:X1", "GID1KO:X2", "GID1KO:X3",
    "GID2KO:X1", "GID2KO:X2", "GID2KO:X3",
    "GID7KO:X1", "GID7KO:X2", "GID7KO:X3",
    "GID8KO:X1", "GID8KO:X2", "GID8KO:X3",
    "GID9KO:X1", "GID9KO:X2", "GID9KO:X3"),
  p.value = 0.05,
  number = Inf)
timecourse_genes <- rownames(tt)

# Timecourse genes (consecitive timepoints) ------------------------------------

j <- colnames(cfit)[1:18]
hm <- pheatmap(
  coef(cfit)[timecourse_genes, j],
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(
    -max(abs(coef(cfit)[timecourse_genes, j])),
    max(abs(coef(cfit)[timecourse_genes, j])),
    length.out = 101),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    comparison = factor(
      c(rep("6_vs_3", 6), rep("9_vs_6", 6), rep("12_vs_9", 6)),
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    cell_line = factor(
      sapply(strsplit(j, "\\."), "[[", 1),
      levels(y$samples$cell_line)),
    row.names = j),
  main = "Timecourse genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    comparison = setNames(
      timepoint_colours[c("Day_6", "Day_9", "Day_12")],
      c("6_vs_3", "9_vs_6", "12_vs_9"))),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  show_rownames = FALSE)
# TODO: How many clusters?
cutree_rows <- 4
gc <- cutree(hm$tree_row, k = cutree_rows)
gc_colours <- setNames(
  palette.colors(cutree_rows, "Okabe-Ito"),
  seq_len(cutree_rows))
pheatmap(
  coef(cfit)[timecourse_genes, j],
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(
    -max(abs(coef(cfit)[timecourse_genes, j])),
    max(abs(coef(cfit)[timecourse_genes, j])),
    length.out = 101),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    comparison = factor(
      c(rep("6_vs_3", 6), rep("9_vs_6", 6), rep("12_vs_9", 6)),
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    cell_line = factor(
      sapply(strsplit(j, "\\."), "[[", 1),
      levels(y$samples$cell_line)),
    row.names = j),
  main = "Timecourse genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    comparison = setNames(
      timepoint_colours[c("Day_6", "Day_9", "Day_12")],
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    cluster = gc_colours,
    lasonder_gene = c("FALSE" = "white", "TRUE" = "red")),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  annotation_row = data.frame(
    cluster = factor(gc),
    lasonder_gene = factor(names(gc) %in% lasonder_genes),
    row.names = names(gc)))

# TODO: If it's useful to highlight Lasonder genes in timecourse heatmap then
#       need to include it in DE_analysis.R.
table(gc, names(gc) %in% lasonder_genes)

pheatmap(
  coef(cfit)[timecourse_genes, j],
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(
    -max(abs(coef(cfit)[timecourse_genes, j])),
    max(abs(coef(cfit)[timecourse_genes, j])),
    length.out = 101),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    comparison = factor(
      c(rep("6_vs_3", 6), rep("9_vs_6", 6), rep("12_vs_9", 6)),
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    cell_line = factor(
      sapply(strsplit(j, "\\."), "[[", 1),
      levels(y$samples$cell_line)),
    row.names = j),
  main = "Timecourse genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    comparison = setNames(
      timepoint_colours[c("Day_6", "Day_9", "Day_12")],
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    cluster = gc_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_row = data.frame(cluster = factor(gc), row.names = names(gc)))

# Timecourse genes (relative to Day_3) -----------------------------------------

# NOTE: Using logFCs from comparison to Day_3, i.e.
#       Day_6_vs_Day3, Day_9_vs_Day_3, Day_12_vs_Day_3

# TODO: A more robust way of pulling out the relevant contrasts.
j2 <- colnames(cfit2)[1:18]
hm2 <- pheatmap(
  coef(cfit2)[timecourse_genes, j2],
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(
    -max(abs(coef(cfit2)[timecourse_genes, j2])),
    max(abs(coef(cfit2)[timecourse_genes, j2])),
    length.out = 101),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    comparison = factor(
      c(rep("6_vs_3", 6), rep("9_vs_3", 6), rep("12_vs_3", 6)),
      c("6_vs_3", "9_vs_3", "12_vs_3")),
    cell_line = factor(
      sapply(strsplit(j2, "\\."), "[[", 1),
      levels(y$samples$cell_line)),
    row.names = j2),
  main = "Timecourse genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    comparison = setNames(
      timepoint_colours[c("Day_6", "Day_9", "Day_12")],
      c("6_vs_3", "9_vs_3", "12_vs_3"))),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  show_rownames = FALSE)
# TODO: How many clusters?
cutree_rows2 <- 4
gc2 <- cutree(hm2$tree_row, k = cutree_rows2)
gc2_colours <- setNames(
  palette.colors(cutree_rows2, "Dark2"),
  seq_len(cutree_rows2))
pheatmap(
  coef(cfit2)[timecourse_genes, j2],
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(
    -max(abs(coef(cfit2)[timecourse_genes, j2])),
    max(abs(coef(cfit2)[timecourse_genes, j2])),
    length.out = 101),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    comparison = factor(
      c(rep("6_vs_3", 6), rep("9_vs_3", 6), rep("12_vs_3", 6)),
      c("6_vs_3", "9_vs_3", "12_vs_3")),
    cell_line = factor(
      sapply(strsplit(j2, "\\."), "[[", 1),
      levels(y$samples$cell_line)),
    row.names = j2),
  main = "Timecourse genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    comparison = setNames(
      timepoint_colours[c("Day_6", "Day_9", "Day_12")],
      c("6_vs_3", "9_vs_3", "12_vs_3")),
    cluster = gc2_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  annotation_row = data.frame(cluster = factor(gc2), row.names = names(gc2)))

pheatmap(
  coef(cfit2)[timecourse_genes, j2],
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(
    -max(abs(coef(cfit2)[timecourse_genes, j2])),
    max(abs(coef(cfit2)[timecourse_genes, j2])),
    length.out = 101),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    comparison = factor(
      c(rep("6_vs_3", 6), rep("9_vs_3", 6), rep("12_vs_3", 6)),
      c("6_vs_3", "9_vs_3", "12_vs_3")),
    cell_line = factor(
      sapply(strsplit(j2, "\\."), "[[", 1),
      levels(y$samples$cell_line)),
    row.names = j2),
  main = "Timecourse sex genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    comparison = setNames(
      timepoint_colours[c("Day_6", "Day_9", "Day_12")],
      c("6_vs_3", "9_vs_3", "12_vs_3")),
    cluster = gc2_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_row = data.frame(cluster = factor(gc2), row.names = names(gc2)))

# Timecourse genes (excluding Day_6_vs_Day_6; consecutive timepoints) ----------

# NOTE: Dropping Day_6_vs_Day_3 because a lot of the logFCs are similar across
#       cell lines in this comparison in order  to focus on changes later in
#       the timecourse.

j3 <- grep("Day_3", colnames(cfit)[1:18], invert = TRUE, value = TRUE)
hm3 <- pheatmap(
  coef(cfit)[timecourse_genes, j3],
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(
    -max(abs(coef(cfit)[timecourse_genes, j3])),
    max(abs(coef(cfit)[timecourse_genes, j3])),
    length.out = 101),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    comparison = factor(
      c(rep("9_vs_6", 6), rep("12_vs_9", 6)),
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    cell_line = factor(
      sapply(strsplit(j3, "\\."), "[[", 1),
      levels(y$samples$cell_line)),
    row.names = j3),
  main = "Timecourse genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    comparison = setNames(
      timepoint_colours[c("Day_6", "Day_9", "Day_12")],
      c("6_vs_3", "9_vs_6", "12_vs_9"))),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  show_rownames = FALSE)
# TODO: How many clusters?
cutree_rows3 <- 4
gc3 <- cutree(hm3$tree_row, k = cutree_rows3)
gc3_colours <- setNames(
  palette.colors(cutree_rows3, "R4"),
  seq_len(cutree_rows3))
pheatmap(
  coef(cfit)[timecourse_genes, j3],
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(
    -max(abs(coef(cfit)[timecourse_genes, j3])),
    max(abs(coef(cfit)[timecourse_genes, j3])),
    length.out = 101),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    comparison = factor(
      c(rep("9_vs_6", 6), rep("12_vs_9", 6)),
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    cell_line = factor(
      sapply(strsplit(j3, "\\."), "[[", 1),
      levels(y$samples$cell_line)),
    row.names = j3),
  main = "Timecourse genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    comparison = setNames(
      timepoint_colours[c("Day_6", "Day_9", "Day_12")],
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    cluster = gc3_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  annotation_row = data.frame(cluster = factor(gc3), row.names = names(gc3)))

pheatmap(
  coef(cfit)[timecourse_genes, j3],
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(
    -max(abs(coef(cfit)[timecourse_genes, j3])),
    max(abs(coef(cfit)[timecourse_genes, j3])),
    length.out = 101),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    comparison = factor(
      c(rep("9_vs_6", 6), rep("12_vs_9", 6)),
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    cell_line = factor(
      sapply(strsplit(j3, "\\."), "[[", 1),
      levels(y$samples$cell_line)),
    row.names = j3),
  main = "Timecourse genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    comparison = setNames(
      timepoint_colours[c("Day_6", "Day_9", "Day_12")],
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    cluster = gc3_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_row = data.frame(cluster = factor(gc3), row.names = names(gc3)))

# All genes (consecutive timepoints) -------------------------------------------

all_genes <- rownames(y)
j <- colnames(cfit)[1:18]
hm <- pheatmap(
  coef(cfit)[all_genes, j],
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(
    -max(abs(coef(cfit)[all_genes, j])),
    max(abs(coef(cfit)[all_genes, j])),
    length.out = 101),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    comparison = factor(
      c(rep("6_vs_3", 6), rep("9_vs_6", 6), rep("12_vs_9", 6)),
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    cell_line = factor(
      sapply(strsplit(j, "\\."), "[[", 1),
      levels(y$samples$cell_line)),
    row.names = j),
  main = "All genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    comparison = setNames(
      timepoint_colours[c("Day_6", "Day_9", "Day_12")],
      c("6_vs_3", "9_vs_6", "12_vs_9"))),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  show_rownames = FALSE)
# TODO: How many clusters?
cutree_rows <- 9
gc <- cutree(hm$tree_row, k = cutree_rows)
gc_colours <- setNames(
  palette.colors(cutree_rows, "Okabe-Ito"),
  seq_len(cutree_rows))
pheatmap(
  coef(cfit)[all_genes, j],
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(
    -max(abs(coef(cfit)[all_genes, j])),
    max(abs(coef(cfit)[all_genes, j])),
    length.out = 101),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    comparison = factor(
      c(rep("6_vs_3", 6), rep("9_vs_6", 6), rep("12_vs_9", 6)),
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    cell_line = factor(
      sapply(strsplit(j, "\\."), "[[", 1),
      levels(y$samples$cell_line)),
    row.names = j),
  main = "All genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    comparison = setNames(
      timepoint_colours[c("Day_6", "Day_9", "Day_12")],
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    cluster = gc_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  annotation_row = data.frame(
    cluster = factor(gc),
    DE = factor(names(gc) %in% timecourse_genes),
    row.names = names(gc)))

pheatmap(
  coef(cfit)[all_genes, j],
  color = colorRampPalette(c("blue","white","red"))(100),
  breaks = seq(
    -max(abs(coef(cfit)[timecourse_genes, j])),
    max(abs(coef(cfit)[timecourse_genes, j])),
    length.out = 101),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = data.frame(
    comparison = factor(
      c(rep("6_vs_3", 6), rep("9_vs_6", 6), rep("12_vs_9", 6)),
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    cell_line = factor(
      sapply(strsplit(j, "\\."), "[[", 1),
      levels(y$samples$cell_line)),
    row.names = j),
  main = "Timecourse genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    comparison = setNames(
      timepoint_colours[c("Day_6", "Day_9", "Day_12")],
      c("6_vs_3", "9_vs_6", "12_vs_9")),
    cluster = gc_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_row = data.frame(
    cluster = factor(gc),
    DE = factor(names(gc) %in% timecourse_genes),
    row.names = names(gc)))
