# DE analysis of mini-bulk data for G000396_Danu
# Peter Hickey
# 2024-11-26

# Setup ------------------------------------------------------------------------

library(here)
library(SingleCellExperiment)
library(edgeR)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggrepel)
library(Glimma)
library(scuttle)
library(splines)
library(pheatmap)

# Take a DataFrame with AtomicList columns and return a data.frame where these
# columns have been flattened by paste-ing together the elements separated by
# `sep`.
flattenDF <- function(x, sep = "; ") {
  as.data.frame(
    DataFrame(
      endoapply(x, function(xx) {
        if (!is(xx, "list")) {
          return(xx)
        }
        unstrsplit(as(xx, "CharacterList"), sep = sep)
      }),
      row.names = rownames(x)))
}

# Load data --------------------------------------------------------------------

sce <- readRDS(here("data", "SCEs", "G000396_Danu.preprocessed.SCE.rds"))
# Define counts assay for use with edgeR::SE2DGEList()
# NOTE: Using UMI counts results in a weird voom plot, notably an up-tick on
#       the right hand side of the voom plot.
counts(sce) <- assay(sce, "read_counts")
# NOTE: Only retain relevant colData columns.
cd <- colData(sce)
colData(sce) <- cd[
  ,
  c("cell_line", "timepoint", "biological_replicate", "technical_replicate",
    "sample", "group", "cell_line_rep")]

# Some useful colours
# NOTE: First 5 colours are based on Set1
cell_line_colours <- c(
  "GID1KO" = "#e41a1c",
  "GID2KO" = "#377eb8",
  "GID7KO" = "#4daf4a",
  "GID8KO" = "#984ea3",
  "GID9KO" = "#ff7f00",
  "WT" = "black")
timepoint_colours <- setNames(
  palette.colors(nlevels(sce$timepoint), "Set2"),
  levels(sce$timepoint))
group_colours <- setNames(
  unlist(
    lapply(
      cell_line_colours[1:6],
      function(x) {
        unlist(lapply(x, colorspace::lighten, amount = seq(0, 0.75, 0.25)))
      }
    )
  ),
  levels(sce$group))

# Some useful shapes
cell_line_shapes <- setNames(c(21:25, 1), levels(sce$cell_line))

# Setup DGEList object, filter, and normalize ----------------------------------

y <- SE2DGEList(sce)

# Subset to relevant samples
# NOTE: Not using criteria from preprocessing.Rmd because it is stricter than
#       I want/need for DE analysis. Here, I really just want to remove those
#       libraries with crap library size.
libsize_drop <- isOutlier(colSums(y$counts), type = "lower", log = TRUE)
keep_rep <- !libsize_drop
summary(keep_rep)
y <- y[, keep_rep]
y$samples <- droplevels(y$samples)

# Sum technical replicates
y <- sumTechReps(y, y$samples$sample)

# Filter out lowly-expressed genes
keep <- filterByExpr(y, group = y$samples$group)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]

# NOTE: No need for TMMwsp because 96% of the counts in the unfiltered count
#       matrix are non-zero.
y <- normLibSizes(y, method = "TMM")
ggplot(
  y$samples,
  aes(
    x = norm.factors,
    y = lib.size,
    colour = group,
    label = ifelse(
      isOutlier(y$samples$lib.size, log = TRUE, type = "lower"),
      colnames(y),
      ""))) +
  geom_point() +
  scale_y_log10(labels = scales::comma) +
  scale_x_log10() +
  geom_label_repel(label.size = 0.1, max.overlaps = 100) +
  theme_cowplot() +
  scale_colour_manual(values = group_colours)

# MDS plots --------------------------------------------------------------------

dir.create(here("output", "MDS"))
dir.create(here("output", "Glimma"))

# Original data
# Static plots
dev.off()
pdf(here("output", "MDS", "MDS.pdf"), width = 8, height = 7)
par(mfrow = c(1, 1))
mds <- plotMDS(
  y,
  col = timepoint_colours[as.character(y$samples$timepoint)],
  main = "Overall\nColoured by timepoint")
mds_df <- cbind(
  data.frame(x = mds$x, y = mds$y),
  y$samples)
p <- ggplot(mds_df, aes(x = x, y = y)) +
  geom_point(
    aes(colour = timepoint, fill = timepoint, shape = cell_line),
    size = 4) +
  xlab(
    paste0(
      "Leading logFC dim 1 (",
      round(100 * mds$var.explained[1], 0),
      "%)")) +
  ylab(
    paste0(
      "Leading logFC dim 2 (",
      round(100 * mds$var.explained[2], 0),
      "%)")) +
  ggtitle("Overall") +
  scale_colour_manual(values = timepoint_colours) +
  scale_fill_manual(values = timepoint_colours) +
  scale_shape_manual(values = cell_line_shapes) +
  guides(shape = guide_legend(
    override.aes = list(shape = cell_line_shapes, fill = "darkgrey"))) +
  theme_bw()
print(p)

mds <- plotMDS(
  y[, y$samples$timepoint == "Day_3"],
  col = cell_line_colours[
    as.character(y$samples$cell_line[y$samples$timepoint == "Day_3"])],
  main = "Day 3\nColoured by cell line")
mds_df <- cbind(
  data.frame(x = mds$x, y = mds$y),
  y$samples[y$samples$timepoint == "Day_3", ])
p <- ggplot(mds_df, aes(x = x, y = y)) +
  geom_point(
    aes(colour = cell_line, fill = cell_line, shape = cell_line),
    size = 4) +
  xlab(
    paste0(
      "Leading logFC dim 1 (",
      round(100 * mds$var.explained[1], 0),
      "%)")) +
  ylab(
    paste0(
      "Leading logFC dim 2 (",
      round(100 * mds$var.explained[2], 0),
      "%)")) +
  ggtitle("Day 3") +
  scale_colour_manual(values = cell_line_colours) +
  scale_fill_manual(values = cell_line_colours) +
  scale_shape_manual(values = cell_line_shapes) +
  theme_bw()
print(p)

mds <- plotMDS(
  y[, y$samples$timepoint == "Day_6"],
  col = cell_line_colours[
    as.character(y$samples$cell_line[y$samples$timepoint == "Day_6"])],
  main = "Day 6\nColoured by cell line")
mds_df <- cbind(
  data.frame(x = mds$x, y = mds$y),
  y$samples[y$samples$timepoint == "Day_6", ])
p <- ggplot(mds_df, aes(x = x, y = y)) +
  geom_point(
    aes(colour = cell_line, fill = cell_line, shape = cell_line),
    size = 4) +
  xlab(
    paste0(
      "Leading logFC dim 1 (",
      round(100 * mds$var.explained[1], 0),
      "%)")) +
  ylab(
    paste0(
      "Leading logFC dim 2 (",
      round(100 * mds$var.explained[2], 0),
      "%)")) +
  ggtitle("Day 6") +
  scale_colour_manual(values = cell_line_colours) +
  scale_fill_manual(values = cell_line_colours) +
  scale_shape_manual(values = cell_line_shapes) +
  theme_bw()
print(p)

mds <- plotMDS(
  y[, y$samples$timepoint == "Day_9"],
  col = cell_line_colours[
    as.character(y$samples$cell_line[y$samples$timepoint == "Day_9"])],
  main = "Day 9\nColoured by cell line")
mds_df <- cbind(
  data.frame(x = mds$x, y = mds$y),
  y$samples[y$samples$timepoint == "Day_9", ])
p <- ggplot(mds_df, aes(x = x, y = y)) +
  geom_point(
    aes(colour = cell_line, fill = cell_line, shape = cell_line),
    size = 4) +
  xlab(
    paste0(
      "Leading logFC dim 1 (",
      round(100 * mds$var.explained[1], 0),
      "%)")) +
  ylab(
    paste0(
      "Leading logFC dim 2 (",
      round(100 * mds$var.explained[2], 0),
      "%)")) +
  ggtitle("Day 9") +
  scale_colour_manual(values = cell_line_colours) +
  scale_fill_manual(values = cell_line_colours) +
  scale_shape_manual(values = cell_line_shapes) +
  theme_bw()
print(p)

mds <- plotMDS(
  y[, y$samples$timepoint == "Day_12"],
  col = cell_line_colours[
    as.character(y$samples$cell_line[y$samples$timepoint == "Day_12"])],
  main = "Day 12\nColoured by cell line")
mds_df <- cbind(
  data.frame(x = mds$x, y = mds$y),
  y$samples[y$samples$timepoint == "Day_12", ])
p <- ggplot(mds_df, aes(x = x, y = y)) +
  geom_point(
    aes(colour = cell_line, fill = cell_line, shape = cell_line),
    size = 4) +
  xlab(
    paste0(
      "Leading logFC dim 1 (",
      round(100 * mds$var.explained[1], 0),
      "%)")) +
  ylab(
    paste0(
      "Leading logFC dim 2 (",
      round(100 * mds$var.explained[2], 0),
      "%)")) +
  ggtitle("Day 12") +
  scale_colour_manual(values = cell_line_colours) +
  scale_fill_manual(values = cell_line_colours) +
  scale_shape_manual(values = cell_line_shapes) +
  theme_bw()
print(p)
dev.off()

# Interactive
glimmaMDS(y, html = here("output", "Glimma", "overall.MDS.html"))
glimmaMDS(
  y[, y$samples$timepoint == "Day_3"],
  html = here("output", "Glimma", "Day_3.MDS.html"))
glimmaMDS(
  y[, y$samples$timepoint == "Day_6"],
  html = here("output", "Glimma", "Day_6.MDS.html"))
glimmaMDS(
  y[, y$samples$timepoint == "Day_9"],
  html = here("output", "Glimma", "Day_9.MDS.html"))
glimmaMDS(
  y[, y$samples$timepoint == "Day_12"],
  html = here("output", "Glimma", "Day_12.MDS.html"))
# NOTE: Don't actually require Glimma's '_files' directory, which contains
#       .js and .css files, because Glimma actually generates standalone HTML
#       files. Unfortunately, we have to do this clean up
#       manually (at least for now; see
#       https://github.com/hasaru-k/GlimmaV2/issues/84).
unlink(list.dirs(here("output", "Glimma"), recursive = FALSE), recursive = TRUE)

# Adjusting for timepoint
# NOTE: Some evidence in this plot for the claim that "GIKO cell lines seem to
#       arrest (compared to WT) somewhere between day 6 and 9" in that all
#       the samples cluster together except for the WT.Day_9 and WT.Day_12
#       samples.
y_adj <- removeBatchEffect(
  voom(y),
  batch = y$samples$timepoint,
  group = y$samples$cell_line)
# Static plots
pdf(
  here("output", "MDS", "MDS.adjusted_for_timepoint.pdf"),
  width = 8,
  height = 7)
par(mfrow = c(1, 1))
mds <- plotMDS(
  y_adj,
  col = timepoint_colours[as.character(y$samples$timepoint)],
  main = "Overall adjusted for timepoint\nColoured by timepoint")
mds_df <- cbind(
  data.frame(x = mds$x, y = mds$y),
  y$samples)
p <- ggplot(mds_df, aes(x = x, y = y)) +
  geom_point(
    aes(colour = timepoint, fill = timepoint, shape = cell_line),
    size = 4) +
  xlab(
    paste0(
      "Leading logFC dim 1 (",
      round(100 * mds$var.explained[1], 0),
      "%)")) +
  ylab(
    paste0(
      "Leading logFC dim 2 (",
      round(100 * mds$var.explained[2], 0),
      "%)")) +
  ggtitle("Overall adjusted for timepoint") +
  scale_colour_manual(values = timepoint_colours) +
  scale_fill_manual(values = timepoint_colours) +
  scale_shape_manual(values = cell_line_shapes) +
  guides(shape = guide_legend(
    override.aes = list(shape = cell_line_shapes, fill = "darkgrey"))) +
  theme_bw()
print(p)
mds <- plotMDS(
  y_adj,
  col = cell_line_colours[as.character(y$samples$cell_line)],
  main = "Overall adjusted for timepoint\nColoured by cell_line")
dev.off()

# Interactive plots
glimmaMDS(
  y_adj,
  groups = y$samples,
  labels = rownames(y$samples),
  html = here("output", "Glimma", "overall_adjusted_for_Day.MDS.html"))
# NOTE: Don't actually require Glimma's '_files' directory, which contains
#       .js and .css files, because Glimma actually generates standalone HTML
#       files. Unfortunately, we have to do this clean up
#       manually (at least for now; see
#       https://github.com/hasaru-k/GlimmaV2/issues/84).
unlink(list.dirs(here("output", "Glimma"), recursive = FALSE), recursive = TRUE)

# Multi-level DE analysis ------------------------------------------------------

# NOTE: Using voomLmFit with quality weights after filtering out low-quality
#       replicates

design <- model.matrix(~0 + group, y$samples)
colnames(design) <- sub("group", "", colnames(design))

fit <- voomLmFit(
  y,
  design,
  block = y$samples$cell_line_rep,
  sample.weights = TRUE,
  plot = TRUE,
  keep.EList = TRUE)
# NOTE: Very low correlation (0.02)
fit$correlation

cm <- makeContrasts(
  # Day6 vs. Day3 (within cell line)
  GID1KO.Day_6_vs_GID1KO.Day_3 = GID1KO.Day_6 - GID1KO.Day_3,
  GID2KO.Day_6_vs_GID2KO.Day_3 = GID2KO.Day_6 - GID2KO.Day_3,
  GID7KO.Day_6_vs_GID7KO.Day_3 = GID7KO.Day_6 - GID7KO.Day_3,
  GID8KO.Day_6_vs_GID8KO.Day_3 = GID8KO.Day_6 - GID8KO.Day_3,
  GID9KO.Day_6_vs_GID9KO.Day_3 = GID9KO.Day_6 - GID9KO.Day_3,
  WT.Day_6_vs_WT.Day_3 = WT.Day_6 - WT.Day_3,

  # Day9 vs. Day6 (within cell line)
  GID1KO.Day_9_vs_GID1KO.Day_6 = GID1KO.Day_9 - GID1KO.Day_6,
  GID2KO.Day_9_vs_GID2KO.Day_6 = GID2KO.Day_9 - GID2KO.Day_6,
  GID7KO.Day_9_vs_GID7KO.Day_6 = GID7KO.Day_9 - GID7KO.Day_6,
  GID8KO.Day_9_vs_GID8KO.Day_6 = GID8KO.Day_9 - GID8KO.Day_6,
  GID9KO.Day_9_vs_GID9KO.Day_6 = GID9KO.Day_9 - GID9KO.Day_6,
  WT.Day_9_vs_WT.Day_6 = WT.Day_9 - WT.Day_6,

  # Day12 vs. Day9 (within cell line)
  GID1KO.Day_12_vs_GID1KO.Day_9 = GID1KO.Day_12 - GID1KO.Day_9,
  GID2KO.Day_12_vs_GID2KO.Day_9 = GID2KO.Day_12 - GID2KO.Day_9,
  GID7KO.Day_12_vs_GID7KO.Day_9 = GID7KO.Day_12 - GID7KO.Day_9,
  GID8KO.Day_12_vs_GID8KO.Day_9 = GID8KO.Day_12 - GID8KO.Day_9,
  GID9KO.Day_12_vs_GID9KO.Day_9 = GID9KO.Day_12 - GID9KO.Day_9,
  WT.Day_12_vs_WT.Day_9 = WT.Day_12 - WT.Day_9,

  # Comparisons at Day 3 (cell lines vs. WT)
  GID1KO.Day_3_vs_WT.Day_3 = GID1KO.Day_3 - WT.Day_3,
  GID2KO.Day_3_vs_WT.Day_3 = GID2KO.Day_3 - WT.Day_3,
  GID7KO.Day_3_vs_WT.Day_3 = GID7KO.Day_3 - WT.Day_3,
  GID8KO.Day_3_vs_WT.Day_3 = GID8KO.Day_3 - WT.Day_3,
  GID9KO.Day_3_vs_WT.Day_3 = GID9KO.Day_3 - WT.Day_3,
  KOs.Day_3_vs_WT.Day_3 =
    (GID1KO.Day_3 + GID2KO.Day_3 + GID7KO.Day_3 + GID8KO.Day_3 +
       GID9KO.Day_3) / 5 - WT.Day_3,
  KOs_except_GID1.Day_3_vs_WT.Day_3 =
    (GID2KO.Day_3 + GID7KO.Day_3 + GID8KO.Day_3 + GID9KO.Day_3) / 4 - WT.Day_3,
  KOs_except_GID9.Day_3_vs_WT.Day_3 =
    (GID1KO.Day_3 + GID2KO.Day_3 + GID7KO.Day_3 + GID8KO.Day_3) / 4 - WT.Day_3,

  # Comparisons at Day 6 (cell lines vs. WT)
  GID1KO.Day_6_vs_WT.Day_6 = GID1KO.Day_6 - WT.Day_6,
  GID2KO.Day_6_vs_WT.Day_6 = GID2KO.Day_6 - WT.Day_6,
  GID7KO.Day_6_vs_WT.Day_6 = GID7KO.Day_6 - WT.Day_6,
  GID8KO.Day_6_vs_WT.Day_6 = GID8KO.Day_6 - WT.Day_6,
  GID9KO.Day_6_vs_WT.Day_6 = GID9KO.Day_6 - WT.Day_6,
  KOs.Day_6_vs_WT.Day_6 =
    (GID1KO.Day_6 + GID2KO.Day_6 + GID7KO.Day_6 + GID8KO.Day_6 +
       GID9KO.Day_6) / 5 - WT.Day_6,
  KOs_except_GID1.Day_6_vs_WT.Day_6 =
    (GID2KO.Day_6 + GID7KO.Day_6 + GID8KO.Day_6 + GID9KO.Day_6) / 4 - WT.Day_6,
  KOs_except_GID9.Day_6_vs_WT.Day_6 =
    (GID1KO.Day_6 + GID2KO.Day_6 + GID7KO.Day_6 + GID8KO.Day_6) / 4 - WT.Day_6,

  # Comparisons at Day 9 (cell lines vs. WT)
  GID1KO.Day_9_vs_WT.Day_9 = GID1KO.Day_9 - WT.Day_9,
  GID2KO.Day_9_vs_WT.Day_9 = GID2KO.Day_9 - WT.Day_9,
  GID7KO.Day_9_vs_WT.Day_9 = GID7KO.Day_9 - WT.Day_9,
  GID8KO.Day_9_vs_WT.Day_9 = GID8KO.Day_9 - WT.Day_9,
  GID9KO.Day_9_vs_WT.Day_9 = GID9KO.Day_9 - WT.Day_9,
  KOs.Day_9_vs_WT.Day_9 =
    (GID1KO.Day_9 + GID2KO.Day_9 + GID7KO.Day_9 + GID8KO.Day_9 +
       GID9KO.Day_9) / 5 - WT.Day_9,
  KOs_except_GID1.Day_9_vs_WT.Day_9 =
    (GID2KO.Day_9 + GID7KO.Day_9 + GID8KO.Day_9 + GID9KO.Day_9) / 4 - WT.Day_9,
  KOs_except_GID9.Day_9_vs_WT.Day_9 =
    (GID1KO.Day_9 + GID2KO.Day_9 + GID7KO.Day_9 + GID8KO.Day_9) / 4 - WT.Day_9,

  # Comparisons at Day 12 (cell lines vs. WT)
  GID1KO.Day_12_vs_WT.Day_12 = GID1KO.Day_12 - WT.Day_12,
  GID2KO.Day_12_vs_WT.Day_12 = GID2KO.Day_12 - WT.Day_12,
  GID7KO.Day_12_vs_WT.Day_12 = GID7KO.Day_12 - WT.Day_12,
  GID8KO.Day_12_vs_WT.Day_12 = GID8KO.Day_12 - WT.Day_12,
  GID9KO.Day_12_vs_WT.Day_12 = GID9KO.Day_12 - WT.Day_12,
  KOs.Day_12_vs_WT.Day_12 =
    (GID1KO.Day_12 + GID2KO.Day_12 + GID7KO.Day_12 + GID8KO.Day_12 +
       GID9KO.Day_12) / 5 - WT.Day_12,
  KOs_except_GID1.Day_12_vs_WT.Day_12 =
    (GID2KO.Day_12 + GID7KO.Day_12 + GID8KO.Day_12 + GID9KO.Day_12) / 4 -
    WT.Day_12,
  KOs_except_GID9.Day_12_vs_WT.Day_12 =
    (GID1KO.Day_12 + GID2KO.Day_12 + GID7KO.Day_12 + GID8KO.Day_12) / 4 -
    WT.Day_12,

  # Interactions (Day6 vs. Day3; cell lines vs. WT)
  `(GID1KO.Day_6_vs_GID1KO.Day_3)_vs_(WT.Day_6_vs_WT.Day_3)` =
    (GID1KO.Day_6 - GID1KO.Day_3) - (WT.Day_6 - WT.Day_3),
  `(GID2KO.Day_6_vs_GID2KO.Day_3)_vs_(WT.Day_6_vs_WT.Day_3)` =
    (GID2KO.Day_6 - GID2KO.Day_3) - (WT.Day_6 - WT.Day_3),
  `(GID7KO.Day_6_vs_GID7KO.Day_3)_vs_(WT.Day_6_vs_WT.Day_3)` =
    (GID7KO.Day_6 - GID7KO.Day_3) - (WT.Day_6 - WT.Day_3),
  `(GID8KO.Day_6_vs_GID8KO.Day_3)_vs_(WT.Day_6_vs_WT.Day_3)` =
    (GID8KO.Day_6 - GID8KO.Day_3) - (WT.Day_6 - WT.Day_3),
  `(GID9KO.Day_6_vs_GID9KO.Day_3)_vs_(WT.Day_6_vs_WT.Day_3)` =
    (GID9KO.Day_6 - GID9KO.Day_3) - (WT.Day_6 - WT.Day_3),
  `(KOs.Day_6_vs_KOs.Day_3)_vs_(WT.Day_6_vs_WT.Day_3)` =
    ((GID1KO.Day_6 - GID1KO.Day_3) + (GID2KO.Day_6 - GID2KO.Day_3) +
       (GID7KO.Day_6 - GID7KO.Day_3) + (GID8KO.Day_6 - GID8KO.Day_3) +
       (GID9KO.Day_6 - GID9KO.Day_3)) / 5 -
    (WT.Day_6 - WT.Day_3),
  `(KOs_except_GID1.Day_6_vs_KOs_except_GID1.Day_3)_vs_(WT.Day_6_vs_WT.Day_3)` =
    ((GID2KO.Day_6 - GID2KO.Day_3) + (GID7KO.Day_6 - GID7KO.Day_3) +
       (GID8KO.Day_6 - GID8KO.Day_3) + (GID9KO.Day_6 - GID9KO.Day_3)) / 4 -
    (WT.Day_6 - WT.Day_3),
  `(KOs_except_GID9.Day_6_vs_KOs_except_GID9.Day_3)_vs_(WT.Day_6_vs_WT.Day_3)` =
    ((GID1KO.Day_6 - GID1KO.Day_3) + (GID2KO.Day_6 - GID2KO.Day_3) +
       (GID7KO.Day_6 - GID7KO.Day_3) + (GID8KO.Day_6 - GID8KO.Day_3)) / 4 -
    (WT.Day_6 - WT.Day_3),

  # Interactions (Day9 vs. Day6; cell lines vs. WT)
  `(GID1KO.Day_9_vs_GID1KO.Day_6)_vs_(WT.Day_9_vs_WT.Day_6)` =
    (GID1KO.Day_9 - GID1KO.Day_6) - (WT.Day_9 - WT.Day_6),
  `(GID2KO.Day_9_vs_GID2KO.Day_6)_vs_(WT.Day_9_vs_WT.Day_6)` =
    (GID2KO.Day_9 - GID2KO.Day_6) - (WT.Day_9 - WT.Day_6),
  `(GID7KO.Day_9_vs_GID7KO.Day_6)_vs_(WT.Day_9_vs_WT.Day_6)` =
    (GID7KO.Day_9 - GID7KO.Day_6) - (WT.Day_9 - WT.Day_6),
  `(GID8KO.Day_9_vs_GID8KO.Day_6)_vs_(WT.Day_9_vs_WT.Day_6)` =
    (GID8KO.Day_9 - GID8KO.Day_6) - (WT.Day_9 - WT.Day_6),
  `(GID9KO.Day_9_vs_GID9KO.Day_6)_vs_(WT.Day_9_vs_WT.Day_6)` =
    (GID9KO.Day_9 - GID9KO.Day_6) - (WT.Day_9 - WT.Day_6),
  `(KOs.Day_9_vs_KOs.Day_6)_vs_(WT.Day_9_vs_WT.Day_6)` =
    ((GID1KO.Day_9 - GID1KO.Day_6) + (GID2KO.Day_9 - GID2KO.Day_6) +
       (GID7KO.Day_9 - GID7KO.Day_6) + (GID8KO.Day_9 - GID8KO.Day_6) +
       (GID9KO.Day_9 - GID9KO.Day_6)) / 5 -
    (WT.Day_9 - WT.Day_6),
  `(KOs_except_GID1.Day_9_vs_KOs_except_GID1.Day_6)_vs_(WT.Day_9_vs_WT.Day_6)` =
    ((GID2KO.Day_9 - GID2KO.Day_6) + (GID7KO.Day_9 - GID7KO.Day_6) +
       (GID8KO.Day_9 - GID8KO.Day_6) + (GID9KO.Day_9 - GID9KO.Day_6)) / 4 -
    (WT.Day_9 - WT.Day_6),
  `(KOs_except_GID9.Day_9_vs_KOs_except_GID9.Day_6)_vs_(WT.Day_9_vs_WT.Day_6)` =
    ((GID1KO.Day_9 - GID1KO.Day_6) + (GID2KO.Day_9 - GID2KO.Day_6) +
       (GID7KO.Day_9 - GID7KO.Day_6) + (GID8KO.Day_9 - GID8KO.Day_6) +
       (GID9KO.Day_9 - GID9KO.Day_6)) / 4 -
    (WT.Day_9 - WT.Day_6),

  # Interactions (Day12 vs. Day9; cell lines vs. WT)
  `(GID1KO.Day_12_vs_GID1KO.Day_9)_vs_(WT.Day_12_vs_WT.Day_9)` =
    (GID1KO.Day_12 - GID1KO.Day_9) - (WT.Day_12 - WT.Day_9),
  `(GID2KO.Day_12_vs_GID2KO.Day_9)_vs_(WT.Day_12_vs_WT.Day_9)` =
    (GID2KO.Day_12 - GID2KO.Day_9) - (WT.Day_12 - WT.Day_9),
  `(GID7KO.Day_12_vs_GID7KO.Day_9)_vs_(WT.Day_12_vs_WT.Day_9)` =
    (GID7KO.Day_12 - GID7KO.Day_9) - (WT.Day_12 - WT.Day_9),
  `(GID8KO.Day_12_vs_GID8KO.Day_9)_vs_(WT.Day_12_vs_WT.Day_9)` =
    (GID8KO.Day_12 - GID8KO.Day_9) - (WT.Day_12 - WT.Day_9),
  `(GID9KO.Day_12_vs_GID9KO.Day_9)_vs_(WT.Day_12_vs_WT.Day_9)` =
    (GID9KO.Day_12 - GID9KO.Day_9) - (WT.Day_12 - WT.Day_9),
  `(KOs.Day_12_vs_KOs.Day_9)_vs_(WT.Day_12_vs_WT.Day_9)` =
    ((GID1KO.Day_12 - GID1KO.Day_9) + (GID2KO.Day_12 - GID2KO.Day_9) +
       (GID7KO.Day_12 - GID7KO.Day_9) + (GID8KO.Day_12 - GID8KO.Day_9) +
       (GID9KO.Day_12 - GID9KO.Day_9)) / 5 -
    (WT.Day_12 - WT.Day_9),
  `(KOs_except_GID1.Day_12_vs_KOs_except_GID1.Day_9)_vs_(WT.Day_12_vs_WT.Day_9)` =
    ((GID2KO.Day_12 - GID2KO.Day_9) + (GID7KO.Day_12 - GID7KO.Day_9) +
       (GID8KO.Day_12 - GID8KO.Day_9) + (GID9KO.Day_12 - GID9KO.Day_9)) / 4 -
    (WT.Day_12 - WT.Day_9),
  `(KOs_except_GID9.Day_12_vs_KOs_except_GID9.Day_9)_vs_(WT.Day_12_vs_WT.Day_9)` =
    ((GID1KO.Day_12 - GID1KO.Day_9) + (GID2KO.Day_12 - GID2KO.Day_9) +
       (GID7KO.Day_12 - GID7KO.Day_9) + (GID8KO.Day_12 - GID8KO.Day_9) +
       (GID9KO.Day_12 - GID9KO.Day_9)) / 4 -
    (WT.Day_12 - WT.Day_9),

  levels = design)
cfit <- contrasts.fit(fit, cm)
cfit <- eBayes(cfit)

# Outputs of Multi-level DE analysis -------------------------------------------

# DEG lists as CSVs
dir.create(here("output", "DEGs"))
l_deg_summary_df <- lapply(colnames(cfit), function(j) {
  message(j)
  tt <- topTable(cfit, j, number = Inf, sort.by = "none")
  if (grepl("KOs", j)) {
    # Add pairwise logFCs to the 'average contrast' to enable post-hoc
    # filtering of genes based on 'consistent' logFCs in the pairwise
    # comparisons.
    if (grepl("except_GID1", j)) {
      pairwise_coefs <- sapply(
        setdiff(levels(y$samples$cell_line), c("WT", "GID1KO")),
        function(cl)  gsub("KOs\\_except\\_GID1", cl, j))
      n_pairwise <- 4
    } else if (grepl("except_GID9", j)) {
      pairwise_coefs <- sapply(
        setdiff(levels(y$samples$cell_line), c("WT", "GID9KO")),
        function(cl)  gsub("KOs\\_except\\_GID9", cl, j))
      n_pairwise <- 4
    } else {
      pairwise_coefs <- sapply(
        setdiff(levels(y$samples$cell_line), "WT"),
        function(cl)  gsub("KOs", cl, j))
      n_pairwise <- 5
    }
    names(pairwise_coefs) <- paste0(pairwise_coefs, ".logFC")
    pairwise_lfcCons <- lapply(pairwise_coefs, function(coef) {
      topTable(cfit, coef = coef, sort.by = "none", n = Inf)$logFC
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
    here("output", "DEGs", paste0(j, ".DEGs.csv")))

  return(deg_summary_df)
})
# Summarise the number of DEGs in each comparison
do.call(rbind, l_deg_summary_df)

# Glimma MA plots
dir.create(here("output", "Glimma"))
lapply(colnames(cfit), function(j) {
  message(j)
  glimmaMA(
    x = cfit,
    dge = y,
    coef = j,
    anno = y$genes[, c("GENEID", "Name", "description")],
    sample.cols = unname(group_colours[as.character(y$samples$group)]),
    html = here(
      "output",
      "Glimma",
      paste0(j, ".html")),
    main = j)
})
# NOTE: Don't actually require Glimma's '_files' directory, which contains
#       .js and .css files, because Glimma actually generates standalone HTML
#       files. Unfortunately, we have to do this clean up
#       manually (at least for now; see
#       https://github.com/hasaru-k/GlimmaV2/issues/84).
unlink(list.dirs(here("output", "Glimma"), recursive = FALSE), recursive = TRUE)

# Heatmaps
lcpm <- edgeR::cpm(y, log = TRUE)
dir.create(here("output", "heatmaps"))
lapply(colnames(cfit), function(j) {
  message(j)

  # 1. Samples ordered by `group`.
  pheatmap(
    lcpm[
      rownames(topTable(cfit, coef = j, n = 50, p.value = 0.05)),
      order(y$samples$group)],
    scale = "row",
    color = colorRampPalette(c("blue","white","red"))(100),
    fontsize_row = 6,
    fontsize_col = 5,
    fontsize = 6,
    annotation_col = y$samples[, "group", drop = FALSE],
    main = gsub("_vs_", " vs. ", j),
    annotation_colors = list(group = group_colours),
    angle_col = 45,
    treeheight_row = 30,
    treeheight_col = 30,
    cluster_cols = FALSE,
    filename = here("output", "heatmaps", paste0(j, ".ordered.pdf")),
    width = 12,
    height = 12)

  # 2. Samples clustered by expression pattern.
  pheatmap(
    lcpm[rownames(topTable(cfit, coef = j, n = 50, p.value = 0.05)), ],
    scale = "row",
    color = colorRampPalette(c("blue","white","red"))(100),
    fontsize_row = 6,
    fontsize_col = 5,
    fontsize = 6,
    annotation_col = y$samples[, c("timepoint", "cell_line")],
    main = gsub("_vs_", " vs. ", j),
    annotation_colors = list(
      cell_line = cell_line_colours,
      timepoint = timepoint_colours),
    angle_col = 45,
    treeheight_row = 30,
    treeheight_col = 30,
    cluster_cols = TRUE,
    filename = here("output", "heatmaps", paste0(j, ".clustered.pdf")),
    width = 12,
    height = 12)

  # 3. Samples subsetted to those involved in the comparison and then clustered
  #    by expression pattern.
  jj <- y$samples$group %in% names(cm[, j][cm[, j] != 0])
  pheatmap(
    lcpm[rownames(topTable(cfit, coef = j, n = 50, p.value = 0.05)), jj],
    scale = "row",
    color = colorRampPalette(c("blue","white","red"))(100),
    fontsize_row = 6,
    fontsize_col = 5,
    fontsize = 6,
    annotation_col = droplevels(y$samples[jj, c("timepoint", "cell_line")]),
    main = gsub("_vs_", " vs. ", j),
    annotation_colors = list(
      cell_line = cell_line_colours[
        levels(droplevels(y$samples$cell_line[jj]))],
      timepoint = timepoint_colours[
        levels(droplevels(y$samples$timepoint[jj]))]),
    angle_col = 45,
    treeheight_row = 30,
    treeheight_col = 30,
    cluster_cols = TRUE,
    filename = here(
      "output",
      "heatmaps",
      paste0(j, ".clustered.only_relevant_samples.pdf")),
    width = 12,
    height = 12)
})

# Gene set analyses ------------------------------------------------------------

# Load .gaf file with GO annotations.
gaf <- GOSemSim::read.gaf(
  here("data", "annotation", "PlasmoDB-66_Pfalciparum3D7_GO.gaf.gzip"))
idx <- ids2indices(
  split(gaf$TERM2GENE$Gene, gaf$TERM2GENE$GO),
  id = unlist(cfit$genes$GENEID))
gene.pathway <- data.frame(
  GeneID = gaf$TERM2GENE$Gene,
  PathwayID = gaf$TERM2GENE$GO)
pathway.names <- data.frame(
  PathwayID = gaf$TERM2NAME$GOID,
  Description = gaf$TERM2NAME$TERM)

# NOTE: kegga() performs over-representation analyses (i.e. hypergeometric
#       test) of gene lists.
#       I _think_ it is a competitive test.
# NOTE: This is a way to hack kegga() to perform goana()-style analysis
#       for non-supported species (such as plasmodium), which Alex Garnham
#       showed me on 2023-12-05.
lapply(colnames(cfit), function(j) {
  tt <- topTable(cfit, j, n = Inf, p.value = 0.05)
  keg_up <- kegga.default(
    de = unlist(tt$GENEID[tt$logFC > 0]),
    universe = unlist(cfit$genes$GENEID),
    gene.pathway = gene.pathway,
    pathway.names = pathway.names)
  keg_down <- kegga.default(
    de = unlist(tt$GENEID[tt$logFC < 0]),
    universe = unlist(cfit$genes$GENEID),
    gene.pathway = gene.pathway,
    pathway.names = pathway.names)
  write.csv(
    topKEGG(keg_up, p.value = 0.05, n = Inf),
    here("output", "DEGs", paste0(j, ".upregulated.kegga_with_GO.csv")))
  write.csv(
    topKEGG(keg_down, p.value = 0.05, n = Inf),
    here("output", "DEGs", paste0(j, ".downregulated.kegga_with_GO.csv")))
})

# NOTE: https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#gene-set-testing-with-camera
#       suggests using camera() rather than mroast() [or fry()] when doing
#       'fishing' analyses.

# NOTE: camera() is a competitive gene set test that accounts for inter-gene
#       correlation.
# NOTE: camera() can't incorporate random effect but cameraPR() can
#       (https://support.bioconductor.org/p/78299/).
lapply(colnames(cfit), function(j) {
  message(j)
  campr <- cameraPR(statistic = cfit$t[, j], index = idx)
  campr$GO <- rownames(campr)
  campr <- dplyr::left_join(campr, gaf$TERM2NAME, by = c("GO" = "GOID"))
  rownames(campr) <- campr$GO
  campr <- campr[, c("TERM", "NGenes", "Direction", "PValue", "FDR")]
  write.csv(
    campr,
    here("output", "DEGs", paste0(j, ".camera.csv")))
})

# Lasonder sex-linked gene expression ------------------------------------------

# See email from Danu on 2024-02-01:

# "A previous study isolated RNAs that were highly expressed in either males or
# females, and therefore based off their ratio in either sex classified genes
# as being male or female specific. I’ve highlighted the Top100 male or female
# genes in the attached file. Could we have a look at how these genes react in
# the GIDKO samples compared to the wild type? This file is the Lasonder 2016 set."

dir.create(here("output", "Lasonder"))

male_genes_tbl <- readxl::read_excel(
  here("data", "gene_lists", "Lasonder 2016_Top100 male and female genes.xlsx"),
  range = "A2:E99")
male_genes_tbl[!male_genes_tbl$`Male genes` %in% rownames(y), ] |>
  knitr::kable(caption = "Male genes not tested in DE analysis.")
male_genes_tbl[
  !male_genes_tbl$`Male genes` %in% rownames(y) &
    male_genes_tbl$`Male genes` %in% rownames(sce), ] |>
  knitr::kable(caption = "Male genes not tested in DE analysis but in dataset.")
female_genes_tbl <- readxl::read_excel(
  here("data", "gene_lists", "Lasonder 2016_Top100 male and female genes.xlsx"),
  range = "G2:K99")
female_genes_tbl[!female_genes_tbl$`Female genes` %in% rownames(y), ] |>
  knitr::kable(caption = "Female genes not tested in DE analysis.")
female_genes_tbl[
  !female_genes_tbl$`Female genes` %in% rownames(y) &
    female_genes_tbl$`Female genes` %in% rownames(sce), ] |>
  knitr::kable(caption = "Female genes not tested in DE analysis but in dataset.")

# Male genes
male_genes <- intersect(male_genes_tbl$`Male genes`, rownames(lcpm))
# 1. Samples ordered by `group`.
pheatmap(
  lcpm[male_genes, order(y$samples$group)],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = y$samples[, "group", drop = FALSE],
  main = "Lasonder male genes",
  annotation_colors = list(group = group_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  filename = here("output", "Lasonder", "Lasonder_male_genes.ordered.pdf"),
  width = 12,
  height = 12)
# 2. Samples clustered by expression pattern.
pheatmap(
  lcpm[male_genes, ],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = y$samples[, c("timepoint", "cell_line")],
  main = "Lasonder male genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    timepoint = timepoint_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = TRUE,
  filename = here("output", "Lasonder", "Lasonder_male_genes.clustered.pdf"),
  width = 12,
  height = 12)

# Female genes
female_genes <- intersect(female_genes_tbl$`Female genes`, rownames(lcpm))
# 1. Samples ordered by `group`.
pheatmap(
  lcpm[female_genes, ],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = y$samples[, "group", drop = FALSE],
  main = "Lasonder female genes",
  annotation_colors = list(group = group_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  filename = here("output", "Lasonder", "Lasonder_female_genes.ordered.pdf"),
  width = 12,
  height = 12)
# 2. Samples clustered by expression pattern.
pheatmap(
  lcpm[female_genes, order(y$samples$group)],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = y$samples[, c("timepoint", "cell_line")],
  main = "Lasonder female genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    timepoint = timepoint_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = TRUE,
  filename = here("output", "Lasonder", "Lasonder_female_genes.clustered.pdf"),
  width = 12,
  height = 12)

# Male and female
lasonder_genes <- c(male_genes, female_genes)
# 1. Samples ordered by `group`.
pheatmap(
  lcpm[lasonder_genes, order(y$samples$group)],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = y$samples[, "group", drop = FALSE],
  annotation_row = data.frame(
    sex = c(
      rep("Male", length(male_genes)),
      rep("Female", length(female_genes))),
    row.names = lasonder_genes),
  main = "Lasonder sex genes",
  annotation_colors = list(
    group = group_colours,
    sex = c(Male = "skyblue3", Female = "deeppink")),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  filename = here("output", "Lasonder", "Lasonder_sex_genes.ordered.pdf"),
  width = 12,
  height = 12)
# 2. Samples clustered by expression pattern.
pheatmap(
  lcpm[lasonder_genes, ],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = y$samples[, c("timepoint", "cell_line")],
  annotation_row = data.frame(
    sex = c(
      rep("Male", length(male_genes)),
      rep("Female", length(female_genes))),
    row.names = lasonder_genes),
  main = "Lasonder sex genes",
  annotation_colors = list(
    group = group_colours,
    sex = c(Male = "skyblue3", Female = "deeppink")),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = TRUE,
  filename = here("output", "Lasonder", "Lasonder_sex_genes.clustered.pdf"),
  width = 12,
  height = 12)

# Formal gene set analysis of Lasonder sex genes.
lasonder_index <- ids2indices(
  list(male = male_genes, female = female_genes, lasonder = lasonder_genes),
  rownames(cfit))
pdf(
  here("output", "Lasonder", "Lasonder.barcodeplots.pdf"),
  width = 5,
  height = 5)
l_of_fry <- lapply(colnames(cfit), function(j) {
  message(j)
  # TODO: Incorporate gene.weights? Something like log(Ratio) where Ratio comes
  #       from male_genes_tbl and female_genes_tbl?
  barcodeplot(
    statistics = cfit$t[, j],
    index = lasonder_index$female,
    index2 = lasonder_index$male,
    xlab = "t",
    col.bars = c("deeppink", "skyblue3"),
    main = j)
  barcodeplot(
    statistics = cfit$t[, j],
    index = lasonder_index$lasonder,
    xlab = "t",
    main = j,
    sub = "Lasander male and female genes")
  # NOTE: fry() is a self-contained gene set test.
  fry(
    y = cfit$EList,
    index = lasonder_index,
    design = design,
    contrast = cm[, j, drop = FALSE])
})
dev.off()
names(l_of_fry) <- colnames(cfit)
# TODO: Male genes start to 'flatline' in the barcode plots by the time we get
#       to the Day_12_vs_Day_9 comparisons. This is more obvious in some cell
#       lines (e.g., WTs) than others; discuss with Danu.
# TODO: The KO vs WT comparisons at Day 9 and Day 12 show downregulation of
#       sex-specific genes, particularly female-specific genes, in KOs;
#       discuss with Danu.

fry_female_df <- do.call(
  rbind,
  lapply(l_of_fry, function(fry) {
    fry["female", , drop = FALSE]
  }))
fry_male_df <- do.call(
  rbind,
  lapply(l_of_fry, function(fry) {
    fry["male", , drop = FALSE]
  }))
fry_lasonnder_df <- do.call(
  rbind,
  lapply(l_of_fry, function(fry) {
    fry["lasonder", , drop = FALSE]
  }))

# TODO: Look at the 'most significant' comparisons; discuss with Danu.
head(
  fry_lasonnder_df[
    order(
      fry_lasonnder_df$Direction,
      fry_lasonnder_df$FDR,
      decreasing = c(FALSE, FALSE),
      method = "radix"), ])
head(
  fry_lasonnder_df[
    order(
      fry_lasonnder_df$Direction,
      fry_lasonnder_df$FDR,
      decreasing = c(TRUE, FALSE),
      method = "radix"), ])

# Guerreiro transcriptional repression program gene expression -----------------

# See email from Danu on 2024-02-01:

# "We would like to look at genes that are known to be regulated by a
#  transcriptional repression programme in Plasmodium. I’ve attached the genes
#  of interest in here as well. Could we look at how they relate to the GIDKO
#  parasites at Day3/6/9/12? I’m sure they would be altered at Day 12 in the KO
#  samples compared to WT as those parasites are arrested at an earlier
#  developmental stage. This file is the Guerreiro 2014 set."

dir.create(here("output", "Guerreiro"))

# TODO: Is there a version of this dataset with logFCs?
transcriptional_repression_tbl <- readxl::read_excel(
  here(
    "data",
    "gene_lists",
    "Guerreiro 2014_Genes common to DOZI and CITH IP and down regulated in DOZI KO.xlsx"))
transcriptional_repression <- c(
  colnames(transcriptional_repression_tbl),
  unlist(transcriptional_repression_tbl))
names(transcriptional_repression) <- NULL
transcriptional_repression <- transcriptional_repression[
  !is.na(transcriptional_repression)]

data.frame(gene = setdiff(transcriptional_repression, rownames(y))) |>
  knitr::kable(
    caption = "Guerreiro transcriptional repression program not tested in DE analysis.")

data.frame(gene = intersect(
  setdiff(transcriptional_repression, rownames(y)),
  rownames(sce))) |>
  knitr::kable(
    caption = "Guerreiro transcriptional repression program not tested in DE analysis but in dataset.")

transcriptional_repression <- intersect(
  transcriptional_repression, rownames(y))

# 1. Samples ordered by `group`.
pheatmap(
  lcpm[transcriptional_repression, order(y$samples$group)],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = y$samples[, "group", drop = FALSE],
  main = "Guerreiro genes",
  annotation_colors = list(group = group_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  filename = here("output", "Guerreiro", "Guerreiro_genes.ordered.pdf"),
  width = 12,
  height = 12)
# 2. Samples clustered by expression pattern.
pheatmap(
  lcpm[transcriptional_repression, ],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = y$samples[, c("timepoint", "cell_line")],
  main = "Guerreiro genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    timepoint = timepoint_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = TRUE,
  filename = here("output", "Guerreiro", "Guerreiro_genes.clustered.pdf"),
  width = 12,
  height = 12)

# Formal gene set analysis of Lasonder sex genes.
guerreiro_index <- ids2indices(
  list(guerreiro = transcriptional_repression),
  rownames(cfit))
pdf(
  here("output", "Guerreiro", "Guerreiro.barcodeplots.pdf"),
  width = 5,
  height = 5)
l_of_fry <- lapply(colnames(cfit), function(j) {
  message(j)
  # TODO: Incorporate gene.weights? Don't have logFCs, but will ask.
  barcodeplot(
    statistics = cfit$t[, j],
    index = guerreiro_index[[1]],
    xlab = "t",
    main = j,
    sub = "Guerreiro transcriptional repression program genes")
  # NOTE: fry() is a self-contained gene set test.
  fry(
    y = cfit$EList,
    index = guerreiro_index,
    design = design,
    contrast = cm[, j, drop = FALSE])
})
dev.off()
names(l_of_fry) <- colnames(cfit)
fry_guerreiro_df <- do.call(
  rbind,
  lapply(l_of_fry, function(fry) {
    fry["guerreiro", , drop = FALSE]
  }))

# TODO: Look at the 'most significant' comparisons; discuss with Danu.
head(
  fry_guerreiro_df[
    order(
      fry_guerreiro_df$Direction,
      fry_lasonnder_df$PValue,
      decreasing = c(FALSE, FALSE),
      method = "radix"), ])
head(
  fry_guerreiro_df[
    order(
      fry_guerreiro_df$Direction,
      fry_lasonnder_df$PValue,
      decreasing = c(TRUE, FALSE),
      method = "radix"), ])

# Processing bodies (P-body), mRNA binding proteins & Translation gene sets ----

dir.create(here("output", "Pbody_mRNA_translation"))

# See email from Danu on 2024-11-25:

# "Here are the lists of the genes we are looking at. There are 3 gene sets all together [Translation, mRNA binding proteins, Processing bodies (P-body)].
# Could we assemble barcode plots for these sets similar to what you previously did for the Guerreiro (translational repression genes) and the Lasonder (male and female specific markers)?
# The data sets we’re most interested in comparing are:
#
# - All knockouts vs WT at Day 3
# - All knockouts vs WT at Day 6
# - All knockouts vs WT at Day 9
# - All knockouts vs WT at Day 12"

pbody_tbl <- readxl::read_excel(
  here("data", "gene_lists", "/Danu_Barcode_updated.xlsx"),
  sheet = "Pbody")
pbody_tbl[!pbody_tbl$Plasmodium_Pbody %in% rownames(y), ] |>
  knitr::kable(caption = "'Processing body' genes not tested in DE analysis.")
pbody_tbl[
  !pbody_tbl$Plasmodium_Pbody %in% rownames(y) &
    pbody_tbl$Plasmodium_Pbody %in% rownames(sce), ] |>
  knitr::kable(
    caption = "'Processing body' genes not tested in DE analysis but in dataset.")
mrna_tbl <- readxl::read_excel(
  here("data", "gene_lists", "/Danu_Barcode_updated.xlsx"),
  sheet = "mRNA")
mrna_tbl[!mrna_tbl$`Plasmodium Identifiers` %in% rownames(y), ] |>
  knitr::kable(
    caption = "'mRNA binding protein' genes not tested in DE analysis.")
mrna_tbl[
  !mrna_tbl$`Plasmodium Identifiers` %in% rownames(y) &
    mrna_tbl$`Plasmodium Identifiers` %in% rownames(sce), ] |>
  knitr::kable(
    caption = "'mRNA binding protein' genes not tested in DE analysis but in dataset.")
translation_tbl <- readxl::read_excel(
  here("data", "gene_lists", "/Danu_Barcode_updated.xlsx"),
  sheet = "Translation")
translation_tbl[!translation_tbl$`Gene ID` %in% rownames(y), ] |>
  knitr::kable(
    caption = "'Translation' genes not tested in DE analysis.")
translation_tbl[
  !translation_tbl$`Gene ID` %in% rownames(y) &
    translation_tbl$`Gene ID` %in% rownames(sce), ] |>
  knitr::kable(
    caption = "'Translation' genes not tested in DE analysis but in dataset.")

# Processing body genes
pbody_genes <- intersect(pbody_tbl$Plasmodium_Pbody, rownames(lcpm))
# 1. Samples ordered by `group`.
pheatmap(
  lcpm[pbody_genes, order(y$samples$group)],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = y$samples[, "group", drop = FALSE],
  main = "'Processing body' genes",
  annotation_colors = list(group = group_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  filename = here(
    "output",
    "Pbody_mRNA_translation",
    "Processing_body_genes.ordered.pdf"),
  width = 12,
  height = 12)
# 2. Samples clustered by expression pattern.
pheatmap(
  lcpm[pbody_genes, ],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = y$samples[, c("timepoint", "cell_line")],
  main = "Lasonder male genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    timepoint = timepoint_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = TRUE,
  filename = here(
    "output",
    "Pbody_mRNA_translation",
    "Processing_body_genes.clustered.pdf"),
  width = 12,
  height = 12)

# mRNA binding protein genes
mrna_genes <- intersect(mrna_tbl$`Plasmodium Identifiers`, rownames(lcpm))
# 1. Samples ordered by `group`.
pheatmap(
  lcpm[mrna_genes, order(y$samples$group)],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = y$samples[, "group", drop = FALSE],
  main = "'mRNA binding protein' genes",
  annotation_colors = list(group = group_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  filename = here(
    "output",
    "Pbody_mRNA_translation",
    "mRNA_binding_protein_genes.ordered.pdf"),
  width = 12,
  height = 12)
# 2. Samples clustered by expression pattern.
pheatmap(
  lcpm[mrna_genes, ],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = y$samples[, c("timepoint", "cell_line")],
  main = "'mRNA binding protein' genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    timepoint = timepoint_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = TRUE,
  filename = here(
    "output",
    "Pbody_mRNA_translation",
    "mRNA_binding_protein_genes.clustered.pdf"),
  width = 12,
  height = 12)

# Translation genes
translation_genes <- intersect(translation_tbl$`Gene ID`, rownames(lcpm))
# 1. Samples ordered by `group`.
pheatmap(
  lcpm[translation_genes, order(y$samples$group)],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = y$samples[, "group", drop = FALSE],
  main = "'Translation' genes",
  annotation_colors = list(group = group_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = FALSE,
  filename = here(
    "output",
    "Pbody_mRNA_translation",
    "Translation_genes.ordered.pdf"),
  width = 12,
  height = 12)
# 2. Samples clustered by expression pattern.
pheatmap(
  lcpm[translation_genes, ],
  scale = "row",
  color = colorRampPalette(c("blue","white","red"))(100),
  fontsize_row = 6,
  fontsize_col = 5,
  fontsize = 6,
  annotation_col = y$samples[, c("timepoint", "cell_line")],
  main = "'Translation' genes",
  annotation_colors = list(
    cell_line = cell_line_colours,
    timepoint = timepoint_colours),
  angle_col = 45,
  treeheight_row = 30,
  treeheight_col = 30,
  cluster_cols = TRUE,
  filename = here(
    "output",
    "Pbody_mRNA_translation",
    "Translation_genes.clustered.pdf"),
  width = 12,
  height = 12)

# Formal gene set analysis of P-body, mRNA binding protein, and translation
# genes
pbody_mrna_translation_index <- ids2indices(
  list(pbody = pbody_genes, mrna = mrna_genes, translation = translation_genes),
  rownames(cfit))
pdf(
  here(
    "output",
    "Pbody_mRNA_translation",
    "Pbody_mRNA_translation.barcodeplots.pdf"),
  width = 5,
  height = 5)
l_of_fry <- lapply(colnames(cfit), function(j) {
  message(j)
  # TODO: Incorporate gene.weights? Something like log(Ratio) where Ratio comes
  #       from male_genes_tbl and female_genes_tbl?
  barcodeplot(
    statistics = cfit$t[, j],
    index = pbody_mrna_translation_index$pbody,
    xlab = "t",
    main = j,
    sub = "'Processing bodies' genes")
  barcodeplot(
    statistics = cfit$t[, j],
    index = pbody_mrna_translation_index$mrna,
    xlab = "t",
    main = j,
    sub = "'mRNA binding protein' genes")
  barcodeplot(
    statistics = cfit$t[, j],
    index = pbody_mrna_translation_index$translation,
    xlab = "t",
    main = j,
    sub = "'Translation' genes")
  # NOTE: fry() is a self-contained gene set test.
  fry(
    y = cfit$EList,
    index = pbody_mrna_translation_index,
    design = design,
    contrast = cm[, j, drop = FALSE])
})
dev.off()
names(l_of_fry) <- colnames(cfit)

fry_pbody_df <- do.call(
  rbind,
  lapply(l_of_fry, function(fry) {
    fry["pbody", , drop = FALSE]
  }))
fry_mrna_df <- do.call(
  rbind,
  lapply(l_of_fry, function(fry) {
    fry["mrna", , drop = FALSE]
  }))
fry_translation_df <- do.call(
  rbind,
  lapply(l_of_fry, function(fry) {
    fry["translation", , drop = FALSE]
  }))

# TODO: Look at the 'most significant' comparisons; discuss with Danu.
head(
  fry_pbody_df[
    order(
      fry_pbody_df$Direction,
      fry_pbody_df$FDR,
      decreasing = c(FALSE, FALSE),
      method = "radix"), ])
head(
  fry_pbody_df[
    order(
      fry_pbody_df$Direction,
      fry_pbody_df$FDR,
      decreasing = c(TRUE, FALSE),
      method = "radix"), ])

head(
  fry_mrna_df[
    order(
      fry_mrna_df$Direction,
      fry_mrna_df$FDR,
      decreasing = c(FALSE, FALSE),
      method = "radix"), ])
head(
  fry_mrna_df[
    order(
      fry_mrna_df$Direction,
      fry_mrna_df$FDR,
      decreasing = c(TRUE, FALSE),
      method = "radix"), ])

head(
  fry_translation_df[
    order(
      fry_translation_df$Direction,
      fry_translation_df$FDR,
      decreasing = c(FALSE, FALSE),
      method = "radix"), ])
head(
  fry_translation_df[
    order(
      fry_translation_df$Direction,
      fry_translation_df$FDR,
      decreasing = c(TRUE, FALSE),
      method = "radix"), ])

# Time course analysis ---------------------------------------------------------

X <- ns(as.integer(y$samples$timepoint), df = 3)
Group <- relevel(y$samples$cell_line, "WT")
design_tc <- model.matrix(~Group * X)
colnames(design_tc) <- sub("Group", "", colnames(design_tc))

fit_tc <- voomLmFit(
  y,
  design_tc,
  block = y$samples$cell_line_rep,
  sample.weights = TRUE,
  plot = TRUE)
fit_tc <- eBayes(fit_tc)

# Outputs of time course analysis ----------------------------------------------

# Summarise the number of DEGs in each comparison
tt_gid1 <- topTable(
  fit_tc,
  coef = c("GID1KO:X1", "GID1KO:X2", "GID1KO:X3"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
tt_gid2 <- topTable(
  fit_tc,
  coef = c("GID2KO:X1", "GID2KO:X2", "GID2KO:X3"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
tt_gid7 <- topTable(
  fit_tc,
  coef = c("GID7KO:X1", "GID7KO:X2", "GID7KO:X3"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
tt_gid8 <- topTable(
  fit_tc,
  coef = c("GID8KO:X1", "GID8KO:X2", "GID8KO:X3"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
tt_gid9 <- topTable(
  fit_tc,
  coef = c("GID9KO:X1", "GID9KO:X2", "GID9KO:X3"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
# NOTE: This is an ANOVA-like test.
tt_any_ko <- topTable(
  fit_tc,
  coef = c(
    "GID1KO:X1", "GID1KO:X2", "GID1KO:X3",
    "GID2KO:X1", "GID2KO:X2", "GID2KO:X3",
    "GID7KO:X1", "GID7KO:X2", "GID7KO:X3",
    "GID8KO:X1", "GID8KO:X2", "GID8KO:X3",
    "GID9KO:X1", "GID9KO:X2", "GID9KO:X3"),
  number = Inf,
  p.value = 1,
  sort.by = "F")
data.frame(
  Sig = c(
    sum(tt_gid1$adj.P.Val < 0.05),
    sum(tt_gid2$adj.P.Val < 0.05),
    sum(tt_gid7$adj.P.Val < 0.05),
    sum(tt_gid8$adj.P.Val < 0.05),
    sum(tt_gid9$adj.P.Val < 0.05),
    sum(tt_any_ko$adj.P.Val < 0.05)),
  NotSig = c(
    sum(tt_gid1$adj.P.Val >= 0.05),
    sum(tt_gid2$adj.P.Val >= 0.05),
    sum(tt_gid7$adj.P.Val >= 0.05),
    sum(tt_gid8$adj.P.Val >= 0.05),
    sum(tt_gid9$adj.P.Val >= 0.05),
    sum(tt_any_ko$adj.P.Val >= 0.05)),
  row.names = c(
    "GID1KO_vs_WT", "GID2KO_vs_WT", "GID7KO_vs_WT",
    "GID8KO_vs_WT", "GID9KO_vs_WT", "any_KO_vs_WT"))

# Plots
lcpm <- cpm(y, log = TRUE)
dir.create(here("output", "timecourse"))
pdf(here("output", "timecourse", "any_KO_vs_WT.pdf"), width = 9, height = 3)
for (g in rownames(tt_any_ko)[tt_any_ko$adj.P.Val < 0.05]) {
  message(g)
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

# TODO: cutree_rows = ?
# TODO: More stringent FDR cutoff for gene selection?
timeCourseHeatmaps <- function(tt, n, cutree_rows = 6, FDR = 0.05) {
  g <- rownames(tt[tt$adj.P.Val < FDR, ])

  # 1. Samples ordered by `group`.
  hm <- pheatmap(
    lcpm[g, order(y$samples$group)],
    scale = "row",
    color = colorRampPalette(c("blue","white","red"))(100),
    fontsize_row = 6,
    fontsize_col = 5,
    fontsize = 6,
    annotation_col = y$samples[, "group", drop = FALSE],
    main = n,
    annotation_colors = list(group = group_colours),
    angle_col = 45,
    treeheight_row = 30,
    treeheight_col = 30,
    cluster_cols = FALSE,
    width = 12,
    height = 12,
    cutree_rows = cutree_rows,
    show_rownames = FALSE,
    silent = TRUE)
  gc <- cutree(hm$tree_row, k = cutree_rows)
  gc_colours <- setNames(
    palette.colors(cutree_rows, "Okabe-Ito"),
    seq_len(cutree_rows))
  pheatmap(
    lcpm[g, order(y$samples$group)],
    scale = "row",
    color = colorRampPalette(c("blue","white","red"))(100),
    fontsize_row = 6,
    fontsize_col = 5,
    fontsize = 6,
    annotation_col = y$samples[, "group", drop = FALSE],
    main = n,
    annotation_colors = list(group = group_colours, cluster = gc_colours),
    angle_col = 45,
    treeheight_row = 30,
    treeheight_col = 30,
    cluster_cols = FALSE,
    filename = here(
      "output",
      "timecourse",
      paste0(n, ".time_course_heatmap.samples_ordered.pdf")),
    width = 12,
    height = 12,
    cutree_rows = 6,
    show_rownames = FALSE,
    annotation_row = data.frame(cluster = factor(gc), row.names = names(gc)))

  # 2. Samples clustered by expression pattern.
  pheatmap(
    lcpm[g, ],
    scale = "row",
    color = colorRampPalette(c("blue","white","red"))(100),
    fontsize_row = 6,
    fontsize_col = 5,
    fontsize = 6,
    annotation_col = y$samples[, c("timepoint", "cell_line")],
    main = n,
    annotation_colors = list(
      cell_line = cell_line_colours,
      timepoint = timepoint_colours,
      cluster = gc_colours),
    angle_col = 45,
    treeheight_row = 30,
    treeheight_col = 30,
    cluster_cols = TRUE,
    filename = here(
      "output",
      "timecourse",
      paste0(n, ".time_course_heatmap.samples_clustered.pdf")),
    width = 12,
    height = 12,
    cutree_rows = 6,
    show_rownames = FALSE,
    annotation_row = data.frame(cluster = factor(gc), row.names = names(gc)))

  # 3. Samples subsetted to those involved in the comparison and then ordered
  #    by `group`.
  if (n != "any_KO_vs_WT") {
    cts <- strsplit(n, "_")[[1]][c(1, 3)]
    jj <- y$samples$cell_line %in% cts
    hm2 <- pheatmap(
      lcpm[g, jj][, order(y$samples$group[jj])],
      scale = "row",
      color = colorRampPalette(c("blue","white","red"))(100),
      fontsize_row = 6,
      fontsize_col = 5,
      fontsize = 6,
      annotation_col = droplevels(y$samples[jj, c("timepoint", "cell_line")]),
      main = gsub("\\_vs\\_", " vs. ", n),
      angle_col = 45,
      treeheight_row = 30,
      treeheight_col = 30,
      cluster_cols = FALSE,
      width = 12,
      height = 12,
      cutree_rows = 6,
      show_rownames = FALSE,
      silent = TRUE)
    gc2 <- cutree(hm2$tree_row, k = cutree_rows)
    gc2_colours <- setNames(
      palette.colors(cutree_rows, "Dark2"),
      seq_len(cutree_rows))
    pheatmap(
      lcpm[g, jj][, order(y$samples$group[jj])],
      scale = "row",
      color = colorRampPalette(c("blue","white","red"))(100),
      fontsize_row = 6,
      fontsize_col = 5,
      fontsize = 6,
      annotation_col = droplevels(y$samples[jj, c("timepoint", "cell_line")]),
      main = gsub("\\_vs\\_", " vs. ", n),
      annotation_colors = list(
        cell_line = cell_line_colours[
          levels(droplevels(y$samples$cell_line[jj]))],
        timepoint = timepoint_colours[
          levels(droplevels(y$samples$timepoint[jj]))],
        cluster2 = gc2_colours),
      angle_col = 45,
      treeheight_row = 30,
      treeheight_col = 30,
      cluster_cols = FALSE,
      filename = here(
        "output",
        "timecourse",
        paste0(
          n,
          ".time_course_heatmap.samples_ordered.only_relevant_samples.pdf")),
      width = 12,
      height = 12,
      cutree_rows = 6,
      show_rownames = FALSE,
      annotation_row = data.frame(
        cluster2 = factor(gc2),
        row.names = names(gc2)))

    val <- list(gc = gc, gc2 = gc2)
  } else {
    val <- list(gc = gc)
  }
  val
}
gcs_gid1 <- timeCourseHeatmaps(tt_gid1, "GID1KO_vs_WT")
gcs_gid2 <- timeCourseHeatmaps(tt_gid2, "GID2KO_vs_WT")
gcs_gid7 <- timeCourseHeatmaps(tt_gid7, "GID7KO_vs_WT")
gcs_gid8 <- timeCourseHeatmaps(tt_gid8, "GID8KO_vs_WT")
gcs_gid9 <- timeCourseHeatmaps(tt_gid9, "GID9KO_vs_WT")
gcs_any_ko <- timeCourseHeatmaps(tt_any_ko, "any_KO_vs_WT")

# Gene set analyses of timecourse data -----------------------------------------

# NOTE: Relies on output of timeCourseHeatmaps()
timeCourseGO <- function(gc) {
  val <- lapply(sort(unique(gc)), function(cl) {
    kegga.default(
      de = names(gc[gc == cl]),
      universe = unlist(fit_tc$genes$GENEID),
      gene.pathway = gene.pathway,
      pathway.names = pathway.names)
  })
  names(val) <- sort(unique(gc))
  val
}
# TODO: gc or gc2?
l_of_keg_gid1 <- timeCourseGO(gcs_gid1$gc2)
l_of_keg_gid2 <- timeCourseGO(gcs_gid2$gc2)
l_of_keg_gid7 <- timeCourseGO(gcs_gid7$gc2)
l_of_keg_gid8 <- timeCourseGO(gcs_gid8$gc2)
l_of_keg_gid9 <- timeCourseGO(gcs_gid9$gc2)
l_of_keg_any_ko <- timeCourseGO(gcs_any_ko$gc)

# An example of results
topKEGG(l_of_keg_any_ko[[1]], n = 5)
# TODO: Make CSVs of GO results; need to settle on a few other decisions (e.g,
#       clustering algorithm and number of clusters of gene patterns) before
#       its worth generating these.

# TODOs ------------------------------------------------------------------------

# - [ ] Is it possible to re-arrange the row dendograms (timecourse stuff) to
#       have the groups in the same order as the legend?
# - [ ] To get exemplars of each cluster of timecourse genes, can use gc in
#       conjunction with the gene list (and then re-order by significance, if
#       necessary).
# - [ ] If it's useful to highlight Lasonder and/or Guerreiro genes in
#       timecourse heatmap then need to include it in the above.
