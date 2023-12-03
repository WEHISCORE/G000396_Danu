# DE analysis of mini-bulk data for G000396_Danu
# Peter Hickey
# 2023-12-04

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
  c("cell_line", "timepoint", "biological_replicate", "technical_replicate", "sample", "group", "cell_line_rep")]

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

# voom using all replicates-----------------------------------------------------

y <- SE2DGEList(sce)

# Sum technical replicates
y <- sumTechReps(y, y$samples$sample)

keep <- filterByExpr(y)
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

glimmaMDS(y)

design <- model.matrix(~0 + group, y$samples)
colnames(design) <- sub("group", "", colnames(design))

v <- voom(y, design, plot = TRUE)
cor <- duplicateCorrelation(v, design, block = v$targets$cell_line_rep)
# NOTE: Very low (0.03)
cor$consensus
v <- voom(
  y,
  design,
  block = v$targets$sample,
  correlation = cor$consensus,
  plot = TRUE)
cor <- duplicateCorrelation(v, design, block = v$targets$cell_line_rep)
# NOTE: Basically doesn't change.
cor$consensus

fit <- lmFit(
  v,
  design,
  block = v$targets$cell_line_rep,
  correlation = cor$consensus)
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

  # Comparisons at Day 6 (cell lines vs. WT)
  GID1KO.Day_6_vs_WT.Day_6 = GID1KO.Day_6 - WT.Day_6,
  GID2KO.Day_6_vs_WT.Day_6 = GID2KO.Day_6 - WT.Day_6,
  GID7KO.Day_6_vs_WT.Day_6 = GID7KO.Day_6 - WT.Day_6,
  GID8KO.Day_6_vs_WT.Day_6 = GID8KO.Day_6 - WT.Day_6,
  GID9KO.Day_6_vs_WT.Day_6 = GID9KO.Day_6 - WT.Day_6,

  # Comparisons at Day 9 (cell lines vs. WT)
  GID1KO.Day_9_vs_WT.Day_9 = GID1KO.Day_9 - WT.Day_9,
  GID2KO.Day_9_vs_WT.Day_9 = GID2KO.Day_9 - WT.Day_9,
  GID7KO.Day_9_vs_WT.Day_9 = GID7KO.Day_9 - WT.Day_9,
  GID8KO.Day_9_vs_WT.Day_9 = GID8KO.Day_9 - WT.Day_9,
  GID9KO.Day_9_vs_WT.Day_9 = GID9KO.Day_9 - WT.Day_9,
  # Comparisons at Day 12 (cell lines vs. WT)
  GID1KO.Day_12_vs_WT.Day_12 = GID1KO.Day_12 - WT.Day_12,
  GID2KO.Day_12_vs_WT.Day_12 = GID2KO.Day_12 - WT.Day_12,
  GID7KO.Day_12_vs_WT.Day_12 = GID7KO.Day_12 - WT.Day_12,
  GID8KO.Day_12_vs_WT.Day_12 = GID8KO.Day_12 - WT.Day_12,
  GID9KO.Day_12_vs_WT.Day_12 = GID9KO.Day_12 - WT.Day_12,

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

  # TODO: Genes that are DE in an cell line vs. WT at start (Day3) via an F-test?
  levels = design)
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)
t(summary(decideTests(fit)))

y1 <- y
v1 <- v
fit1 <- fit

# voom filtering out low-quality replicates ------------------------------------

y <- SE2DGEList(sce)

# Subset to relevant samples
# NOTE: Not using criteria from preprocessing Rmd because it is stricter than
#       I want/need for DE analysis. Here, I really just want to remove those
#       libraries with crap library size.
libsize_drop <- isOutlier(colSums(y$counts), type = "lower", log = TRUE)
keep_rep <- !libsize_drop
summary(keep_rep)
y <- y[, keep_rep]
y$samples <- droplevels(y$samples)

# Sum technical replicates
y <- sumTechReps(y, y$samples$sample)

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

glimmaMDS(y)

design <- model.matrix(~0 + group, y$samples)
colnames(design) <- sub("group", "", colnames(design))

v <- voom(y, design, plot = TRUE)
cor <- duplicateCorrelation(v, design, block = v$targets$cell_line_rep)
# NOTE: Very low (0.04)
cor$consensus
v <- voom(
  y,
  design,
  block = v$targets$sample,
  correlation = cor$consensus,
  plot = TRUE)
cor <- duplicateCorrelation(v, design, block = v$targets$cell_line_rep)
# NOTE: Basically doesn't change.
cor$consensus

fit <- lmFit(
  v,
  design,
  block = v$targets$cell_line_rep,
  correlation = cor$consensus)
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

  # Comparisons at Day 6 (cell lines vs. WT)
  GID1KO.Day_6_vs_WT.Day_6 = GID1KO.Day_6 - WT.Day_6,
  GID2KO.Day_6_vs_WT.Day_6 = GID2KO.Day_6 - WT.Day_6,
  GID7KO.Day_6_vs_WT.Day_6 = GID7KO.Day_6 - WT.Day_6,
  GID8KO.Day_6_vs_WT.Day_6 = GID8KO.Day_6 - WT.Day_6,
  GID9KO.Day_6_vs_WT.Day_6 = GID9KO.Day_6 - WT.Day_6,

  # Comparisons at Day 9 (cell lines vs. WT)
  GID1KO.Day_9_vs_WT.Day_9 = GID1KO.Day_9 - WT.Day_9,
  GID2KO.Day_9_vs_WT.Day_9 = GID2KO.Day_9 - WT.Day_9,
  GID7KO.Day_9_vs_WT.Day_9 = GID7KO.Day_9 - WT.Day_9,
  GID8KO.Day_9_vs_WT.Day_9 = GID8KO.Day_9 - WT.Day_9,
  GID9KO.Day_9_vs_WT.Day_9 = GID9KO.Day_9 - WT.Day_9,
  # Comparisons at Day 12 (cell lines vs. WT)
  GID1KO.Day_12_vs_WT.Day_12 = GID1KO.Day_12 - WT.Day_12,
  GID2KO.Day_12_vs_WT.Day_12 = GID2KO.Day_12 - WT.Day_12,
  GID7KO.Day_12_vs_WT.Day_12 = GID7KO.Day_12 - WT.Day_12,
  GID8KO.Day_12_vs_WT.Day_12 = GID8KO.Day_12 - WT.Day_12,
  GID9KO.Day_12_vs_WT.Day_12 = GID9KO.Day_12 - WT.Day_12,

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

  # TODO: Genes that are DE in an cell line vs. WT at start (Day3) via an F-test?
  levels = design)
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)
t(summary(decideTests(fit)))

y2 <- y
v2 <- v
fit2 <- fit

# voomWithQualityWeights using all replicates ----------------------------------

y <- SE2DGEList(sce)

# Sum technical replicates
y <- sumTechReps(y, y$samples$sample)

keep <- filterByExpr(y)
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

glimmaMDS(y)

design <- model.matrix(~0 + group, y$samples)
colnames(design) <- sub("group", "", colnames(design))

v <- voomWithQualityWeights(y, design, plot = TRUE)
cor <- duplicateCorrelation(v, design, block = v$targets$cell_line_rep)
# NOTE: Very low (0.02)
cor$consensus
v <- voomWithQualityWeights(
  y,
  design,
  block = v$targets$sample,
  correlation = cor$consensus,
  plot = TRUE)
cor <- duplicateCorrelation(v, design, block = v$targets$cell_line_rep)
# NOTE: Basically doesn't change.
cor$consensus

fit <- lmFit(
  v,
  design,
  block = v$targets$cell_line_rep,
  correlation = cor$consensus)

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

  # Comparisons at Day 6 (cell lines vs. WT)
  GID1KO.Day_6_vs_WT.Day_6 = GID1KO.Day_6 - WT.Day_6,
  GID2KO.Day_6_vs_WT.Day_6 = GID2KO.Day_6 - WT.Day_6,
  GID7KO.Day_6_vs_WT.Day_6 = GID7KO.Day_6 - WT.Day_6,
  GID8KO.Day_6_vs_WT.Day_6 = GID8KO.Day_6 - WT.Day_6,
  GID9KO.Day_6_vs_WT.Day_6 = GID9KO.Day_6 - WT.Day_6,

  # Comparisons at Day 9 (cell lines vs. WT)
  GID1KO.Day_9_vs_WT.Day_9 = GID1KO.Day_9 - WT.Day_9,
  GID2KO.Day_9_vs_WT.Day_9 = GID2KO.Day_9 - WT.Day_9,
  GID7KO.Day_9_vs_WT.Day_9 = GID7KO.Day_9 - WT.Day_9,
  GID8KO.Day_9_vs_WT.Day_9 = GID8KO.Day_9 - WT.Day_9,
  GID9KO.Day_9_vs_WT.Day_9 = GID9KO.Day_9 - WT.Day_9,
  # Comparisons at Day 12 (cell lines vs. WT)
  GID1KO.Day_12_vs_WT.Day_12 = GID1KO.Day_12 - WT.Day_12,
  GID2KO.Day_12_vs_WT.Day_12 = GID2KO.Day_12 - WT.Day_12,
  GID7KO.Day_12_vs_WT.Day_12 = GID7KO.Day_12 - WT.Day_12,
  GID8KO.Day_12_vs_WT.Day_12 = GID8KO.Day_12 - WT.Day_12,
  GID9KO.Day_12_vs_WT.Day_12 = GID9KO.Day_12 - WT.Day_12,

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

  # TODO: Genes that are DE in an cell line vs. WT at start (Day3) via an F-test?
  levels = design)
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)
t(summary(decideTests(fit)))

y3 <- y
v3 <- v
fit3 <- fit

# voomWithQualityWeights filtering out low-quality replicates-------------------

y <- SE2DGEList(sce)

# Subset to relevant samples
# NOTE: Not using criteria from preprocessing Rmd because it is stricter than
#       I want/need for DE analysis. Here, I really just want to remove those
#       libraries with crap library size.
libsize_drop <- isOutlier(colSums(y$counts), type = "lower", log = TRUE)
keep_rep <- !libsize_drop
summary(keep_rep)
y <- y[, keep_rep]
y$samples <- droplevels(y$samples)

# Sum technical replicates
y <- sumTechReps(y, y$samples$sample)

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

glimmaMDS(y)

design <- model.matrix(~0 + group, y$samples)
colnames(design) <- sub("group", "", colnames(design))

v <- voomWithQualityWeights(y, design, plot = TRUE)
cor <- duplicateCorrelation(v, design, block = v$targets$cell_line_rep)
# NOTE: Very low (0.02)
cor$consensus
v <- voomWithQualityWeights(
  y,
  design,
  block = v$targets$sample,
  correlation = cor$consensus,
  plot = TRUE)
cor <- duplicateCorrelation(v, design, block = v$targets$cell_line_rep)
# NOTE: Basically doesn't change.
cor$consensus

fit <- lmFit(
  v,
  design,
  block = v$targets$cell_line_rep,
  correlation = cor$consensus)
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

  # Comparisons at Day 6 (cell lines vs. WT)
  GID1KO.Day_6_vs_WT.Day_6 = GID1KO.Day_6 - WT.Day_6,
  GID2KO.Day_6_vs_WT.Day_6 = GID2KO.Day_6 - WT.Day_6,
  GID7KO.Day_6_vs_WT.Day_6 = GID7KO.Day_6 - WT.Day_6,
  GID8KO.Day_6_vs_WT.Day_6 = GID8KO.Day_6 - WT.Day_6,
  GID9KO.Day_6_vs_WT.Day_6 = GID9KO.Day_6 - WT.Day_6,

  # Comparisons at Day 9 (cell lines vs. WT)
  GID1KO.Day_9_vs_WT.Day_9 = GID1KO.Day_9 - WT.Day_9,
  GID2KO.Day_9_vs_WT.Day_9 = GID2KO.Day_9 - WT.Day_9,
  GID7KO.Day_9_vs_WT.Day_9 = GID7KO.Day_9 - WT.Day_9,
  GID8KO.Day_9_vs_WT.Day_9 = GID8KO.Day_9 - WT.Day_9,
  GID9KO.Day_9_vs_WT.Day_9 = GID9KO.Day_9 - WT.Day_9,
  # Comparisons at Day 12 (cell lines vs. WT)
  GID1KO.Day_12_vs_WT.Day_12 = GID1KO.Day_12 - WT.Day_12,
  GID2KO.Day_12_vs_WT.Day_12 = GID2KO.Day_12 - WT.Day_12,
  GID7KO.Day_12_vs_WT.Day_12 = GID7KO.Day_12 - WT.Day_12,
  GID8KO.Day_12_vs_WT.Day_12 = GID8KO.Day_12 - WT.Day_12,
  GID9KO.Day_12_vs_WT.Day_12 = GID9KO.Day_12 - WT.Day_12,

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

  # TODO: Genes that are DE in an cell line vs. WT at start (Day3) via an F-test?
  levels = design)
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)
t(summary(decideTests(fit)))

y4 <- y
v4 <- v
fit4 <- fit

# voomLmFit using all replicates -----------------------------------------------

y <- SE2DGEList(sce)

# Sum technical replicates
y <- sumTechReps(y, y$samples$sample)

keep <- filterByExpr(y)
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

glimmaMDS(y)

design <- model.matrix(~0 + group, y$samples)
colnames(design) <- sub("group", "", colnames(design))

fit <- voomLmFit(
  y,
  design,
  block = y$samples$cell_line_rep,
  plot = TRUE,
  keep.EList = TRUE)
# NOTE: Very low correlation (0.03)
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

  # Comparisons at Day 6 (cell lines vs. WT)
  GID1KO.Day_6_vs_WT.Day_6 = GID1KO.Day_6 - WT.Day_6,
  GID2KO.Day_6_vs_WT.Day_6 = GID2KO.Day_6 - WT.Day_6,
  GID7KO.Day_6_vs_WT.Day_6 = GID7KO.Day_6 - WT.Day_6,
  GID8KO.Day_6_vs_WT.Day_6 = GID8KO.Day_6 - WT.Day_6,
  GID9KO.Day_6_vs_WT.Day_6 = GID9KO.Day_6 - WT.Day_6,

  # Comparisons at Day 9 (cell lines vs. WT)
  GID1KO.Day_9_vs_WT.Day_9 = GID1KO.Day_9 - WT.Day_9,
  GID2KO.Day_9_vs_WT.Day_9 = GID2KO.Day_9 - WT.Day_9,
  GID7KO.Day_9_vs_WT.Day_9 = GID7KO.Day_9 - WT.Day_9,
  GID8KO.Day_9_vs_WT.Day_9 = GID8KO.Day_9 - WT.Day_9,
  GID9KO.Day_9_vs_WT.Day_9 = GID9KO.Day_9 - WT.Day_9,
  # Comparisons at Day 12 (cell lines vs. WT)
  GID1KO.Day_12_vs_WT.Day_12 = GID1KO.Day_12 - WT.Day_12,
  GID2KO.Day_12_vs_WT.Day_12 = GID2KO.Day_12 - WT.Day_12,
  GID7KO.Day_12_vs_WT.Day_12 = GID7KO.Day_12 - WT.Day_12,
  GID8KO.Day_12_vs_WT.Day_12 = GID8KO.Day_12 - WT.Day_12,
  GID9KO.Day_12_vs_WT.Day_12 = GID9KO.Day_12 - WT.Day_12,

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

  # TODO: Genes that are DE in an cell line vs. WT at start (Day3) via an F-test?
  levels = design)
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)
t(summary(decideTests(fit)))

y5 <- y
fit5 <- fit

# voomLmFit filtering out low-quality replicates -------------------------------

y <- SE2DGEList(sce)

# Subset to relevant samples
# NOTE: Not using criteria from preprocessing Rmd because it is stricter than
#       I want/need for DE analysis. Here, I really just want to remove those
#       libraries with crap library size.
libsize_drop <- isOutlier(colSums(y$counts), type = "lower", log = TRUE)
keep_rep <- !libsize_drop
summary(keep_rep)
y <- y[, keep_rep]
y$samples <- droplevels(y$samples)

# Sum technical replicates
y <- sumTechReps(y, y$samples$sample)

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

glimmaMDS(y)

design <- model.matrix(~0 + group, y$samples)
colnames(design) <- sub("group", "", colnames(design))

fit <- voomLmFit(
  y,
  design,
  block = y$samples$cell_line_rep,
  plot = TRUE,
  keep.EList = TRUE)
# NOTE: Very low correlation (0.03)
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

  # Comparisons at Day 6 (cell lines vs. WT)
  GID1KO.Day_6_vs_WT.Day_6 = GID1KO.Day_6 - WT.Day_6,
  GID2KO.Day_6_vs_WT.Day_6 = GID2KO.Day_6 - WT.Day_6,
  GID7KO.Day_6_vs_WT.Day_6 = GID7KO.Day_6 - WT.Day_6,
  GID8KO.Day_6_vs_WT.Day_6 = GID8KO.Day_6 - WT.Day_6,
  GID9KO.Day_6_vs_WT.Day_6 = GID9KO.Day_6 - WT.Day_6,

  # Comparisons at Day 9 (cell lines vs. WT)
  GID1KO.Day_9_vs_WT.Day_9 = GID1KO.Day_9 - WT.Day_9,
  GID2KO.Day_9_vs_WT.Day_9 = GID2KO.Day_9 - WT.Day_9,
  GID7KO.Day_9_vs_WT.Day_9 = GID7KO.Day_9 - WT.Day_9,
  GID8KO.Day_9_vs_WT.Day_9 = GID8KO.Day_9 - WT.Day_9,
  GID9KO.Day_9_vs_WT.Day_9 = GID9KO.Day_9 - WT.Day_9,
  # Comparisons at Day 12 (cell lines vs. WT)
  GID1KO.Day_12_vs_WT.Day_12 = GID1KO.Day_12 - WT.Day_12,
  GID2KO.Day_12_vs_WT.Day_12 = GID2KO.Day_12 - WT.Day_12,
  GID7KO.Day_12_vs_WT.Day_12 = GID7KO.Day_12 - WT.Day_12,
  GID8KO.Day_12_vs_WT.Day_12 = GID8KO.Day_12 - WT.Day_12,
  GID9KO.Day_12_vs_WT.Day_12 = GID9KO.Day_12 - WT.Day_12,

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

  # TODO: Genes that are DE in an cell line vs. WT at start (Day3) via an F-test?
  levels = design)
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)
t(summary(decideTests(fit)))

y6 <- y
fit6 <- fit

# voomLmFit with quality weights using all replicates --------------------------

y <- SE2DGEList(sce)

# Sum technical replicates
y <- sumTechReps(y, y$samples$sample)

keep <- filterByExpr(y)
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

glimmaMDS(y)

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

  # Comparisons at Day 6 (cell lines vs. WT)
  GID1KO.Day_6_vs_WT.Day_6 = GID1KO.Day_6 - WT.Day_6,
  GID2KO.Day_6_vs_WT.Day_6 = GID2KO.Day_6 - WT.Day_6,
  GID7KO.Day_6_vs_WT.Day_6 = GID7KO.Day_6 - WT.Day_6,
  GID8KO.Day_6_vs_WT.Day_6 = GID8KO.Day_6 - WT.Day_6,
  GID9KO.Day_6_vs_WT.Day_6 = GID9KO.Day_6 - WT.Day_6,

  # Comparisons at Day 9 (cell lines vs. WT)
  GID1KO.Day_9_vs_WT.Day_9 = GID1KO.Day_9 - WT.Day_9,
  GID2KO.Day_9_vs_WT.Day_9 = GID2KO.Day_9 - WT.Day_9,
  GID7KO.Day_9_vs_WT.Day_9 = GID7KO.Day_9 - WT.Day_9,
  GID8KO.Day_9_vs_WT.Day_9 = GID8KO.Day_9 - WT.Day_9,
  GID9KO.Day_9_vs_WT.Day_9 = GID9KO.Day_9 - WT.Day_9,
  # Comparisons at Day 12 (cell lines vs. WT)
  GID1KO.Day_12_vs_WT.Day_12 = GID1KO.Day_12 - WT.Day_12,
  GID2KO.Day_12_vs_WT.Day_12 = GID2KO.Day_12 - WT.Day_12,
  GID7KO.Day_12_vs_WT.Day_12 = GID7KO.Day_12 - WT.Day_12,
  GID8KO.Day_12_vs_WT.Day_12 = GID8KO.Day_12 - WT.Day_12,
  GID9KO.Day_12_vs_WT.Day_12 = GID9KO.Day_12 - WT.Day_12,

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

  # TODO: Genes that are DE in an cell line vs. WT at start (Day3) via an F-test?
  levels = design)
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)
t(summary(decideTests(fit)))

y7 <- y
fit7 <- fit

# voomLmFit with quality weights filtering out low-quality replicates ----------

y <- SE2DGEList(sce)

# Subset to relevant samples
# NOTE: Not using criteria from preprocessing Rmd because it is stricter than
#       I want/need for DE analysis. Here, I really just want to remove those
#       libraries with crap library size.
libsize_drop <- isOutlier(colSums(y$counts), type = "lower", log = TRUE)
keep_rep <- !libsize_drop
summary(keep_rep)
y <- y[, keep_rep]
y$samples <- droplevels(y$samples)

# Sum technical replicates
y <- sumTechReps(y, y$samples$sample)

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

glimmaMDS(y)

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

  # Comparisons at Day 6 (cell lines vs. WT)
  GID1KO.Day_6_vs_WT.Day_6 = GID1KO.Day_6 - WT.Day_6,
  GID2KO.Day_6_vs_WT.Day_6 = GID2KO.Day_6 - WT.Day_6,
  GID7KO.Day_6_vs_WT.Day_6 = GID7KO.Day_6 - WT.Day_6,
  GID8KO.Day_6_vs_WT.Day_6 = GID8KO.Day_6 - WT.Day_6,
  GID9KO.Day_6_vs_WT.Day_6 = GID9KO.Day_6 - WT.Day_6,

  # Comparisons at Day 9 (cell lines vs. WT)
  GID1KO.Day_9_vs_WT.Day_9 = GID1KO.Day_9 - WT.Day_9,
  GID2KO.Day_9_vs_WT.Day_9 = GID2KO.Day_9 - WT.Day_9,
  GID7KO.Day_9_vs_WT.Day_9 = GID7KO.Day_9 - WT.Day_9,
  GID8KO.Day_9_vs_WT.Day_9 = GID8KO.Day_9 - WT.Day_9,
  GID9KO.Day_9_vs_WT.Day_9 = GID9KO.Day_9 - WT.Day_9,
  # Comparisons at Day 12 (cell lines vs. WT)
  GID1KO.Day_12_vs_WT.Day_12 = GID1KO.Day_12 - WT.Day_12,
  GID2KO.Day_12_vs_WT.Day_12 = GID2KO.Day_12 - WT.Day_12,
  GID7KO.Day_12_vs_WT.Day_12 = GID7KO.Day_12 - WT.Day_12,
  GID8KO.Day_12_vs_WT.Day_12 = GID8KO.Day_12 - WT.Day_12,
  GID9KO.Day_12_vs_WT.Day_12 = GID9KO.Day_12 - WT.Day_12,

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

  # TODO: Genes that are DE in an cell line vs. WT at start (Day3) via an F-test?
  levels = design)
fit <- contrasts.fit(fit, cm)
fit <- eBayes(fit)
t(summary(decideTests(fit)))

y8 <- y
fit8 <- fit

# Time course analysis ---------------------------------------------------------

# NOTE: Using the results from the
#       'voomLmFit with quality weights filtering out low-quality replicate'
#       analysis.

library(splines)
X <- ns(as.integer(y8$samples$timepoint), df = 3)
Group <- relevel(y8$samples$cell_line, "WT")
design <- model.matrix(~Group * X)
colnames(design) <- sub("Group", "", colnames(design))

fit <- voomLmFit(
  y8,
  design,
  block = y8$samples$cell_line_rep,
  sample.weights = TRUE,
  plot = TRUE)
fit <- eBayes(fit)

fit9 <- fit

# TODO: Could also extract genes that change over time within a given cell line
#       (rather than those that change where there change is different to the
#       change within WT).

# Outputs to share with Danu ---------------------------------------------------

# NOTE: Using the results from the
#       'voomLmFit with quality weights filtering out low-quality replicate'
#       analysis.

# Glimma plots
# TODO: Provide only necessary/useful columns in `anno`.
dir.create(here("tmp", "Glimma"))
lapply(colnames(fit8), function(j) {
  message(j)
  glimmaMA(
    x = fit8,
    dge = y8,
    coef = j,
    sample.cols = unname(group_colours[y8$samples$group]),
    html = here(
      "tmp",
      "Glimma",
      paste0(j, ".html")),
    main = j)
})

# Heatmaps
lcpm <- edgeR::cpm(y8, log = TRUE)
dir.create(here("tmp", "pheatmap"))

# TODO: Might need to use the `filename` argument of pheatmap() because
#       otherwise the last heatmap text gets overprinted.
pdf(here("tmp", "pheatmap", "pheatmap.ordered.pdf"), height = 12, width = 12)
lapply(colnames(fit8), function(j) {
  message(j)

  pheatmap::pheatmap(
    lcpm[
      rownames(topTable(fit8, coef = j, n = 100, p.value = 0.05)),
      order(y8$samples$group)],
    scale = "row",
    color = colorRampPalette(c("blue","white","red"))(100),
    fontsize_row = 6,
    fontsize_col = 5,
    fontsize = 6,
    # annotation_col = y8$samples[, c("cell_line", "timepoint")],
    annotation_col = y8$samples[, c("cell_line", "timepoint", "group")],
    main = gsub("_vs_", " vs. ", j),
    # annotation_colors = list(
    #   cell_line = cell_line_colours,
    #   timepoint = timepoint_colours),
    annotation_colors = list(
      cell_line = cell_line_colours,
      timepoint = timepoint_colours,
      group = group_colours),
    angle_col = 45,
    treeheight_row = 30,
    treeheight_col = 30,
    cluster_cols = FALSE)
})
dev.off()
# TODO: Might need to use the `filename` argument of pheatmap() because
#       otherwise the last heatmap text gets overprinted.
pdf(here("tmp", "pheatmap", "pheatmap.clustered.pdf"), height = 12, width = 12)
lapply(colnames(fit8), function(j) {
  message(j)

  pheatmap::pheatmap(
    lcpm[rownames(topTable(fit8, coef = j, n = 100, p.value = 0.05)), ],
    scale = "row",
    color = colorRampPalette(c("blue","white","red"))(100),
    fontsize_row = 6,
    fontsize_col = 5,
    fontsize = 6,
    # annotation_col = y8$samples[, c("timepoint", "cell_line")],
    annotation_col = y8$samples[, c("timepoint", "cell_line", "group")],
    main = gsub("_vs_", " vs. ", j),
    # annotation_colors = list(
    #   cell_line = cell_line_colours,
    #   timepoint = timepoint_colours),
    annotation_colors = list(
      cell_line = cell_line_colours,
      timepoint = timepoint_colours,
      group = group_colours),
    angle_col = 45,
    treeheight_row = 30,
    treeheight_col = 30)
})
dev.off()

# TODO: Only show samples involved in the contrast in the heatmap?

# Timecourse outputs
tt <- topTable(
  fit9,
  coef = c("GID1KO:X1", "GID1KO:X2", "GID1KO:X3"),
  number = Inf,
  p.value = 0.05)
nrow(tt)
dir.create(here("tmp", "timecourse"))
pdf(here("tmp", "timecourse", "GID1KO_vs_WT.pdf"), width = 9, height = 3)
for (g in rownames(tt)) {
  message(g)
  p <- ggplot(
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
    panel_border() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  print(p)
}
dev.off()

tt <- topTable(
  fit9,
  coef = c("GID2KO:X1", "GID2KO:X2", "GID2KO:X3"),
  number = Inf,
  p.value = 0.05)
nrow(tt)
pdf(here("tmp", "timecourse", "GID2KO_vs_WT.pdf"), width = 9, height = 3)
for (g in rownames(tt)) {
  message(g)
  p <- ggplot(
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
    panel_border() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  print(p)
}
dev.off()

tt <- topTable(
  fit9,
  coef = c("GID7KO:X1", "GID7KO:X2", "GID7KO:X3"),
  number = Inf,
  p.value = 0.05)
nrow(tt)
pdf(here("tmp", "timecourse", "GID7KO_vs_WT.pdf"), width = 9, height = 3)
for (g in rownames(tt)) {
  message(g)
  p <- ggplot(
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
    panel_border() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  print(p)
}
dev.off()

tt <- topTable(
  fit9,
  coef = c("GID8KO:X1", "GID8KO:X2", "GID8KO:X3"),
  number = Inf,
  p.value = 0.05)
nrow(tt)
pdf(here("tmp", "timecourse", "GID8KO_vs_WT.pdf"), width = 9, height = 3)
for (g in rownames(tt)) {
  message(g)
  p <- ggplot(
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
    panel_border() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  print(p)
}
dev.off()

tt <- topTable(
  fit9,
  coef = c("GID9KO:X1", "GID9KO:X2", "GID9KO:X3"),
  number = Inf,
  p.value = 0.05)
nrow(tt)
pdf(here("tmp", "timecourse", "GID9KO_vs_WT.pdf"), width = 9, height = 3)
for (g in rownames(tt)) {
  message(g)
  p <- ggplot(
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
    panel_border() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  print(p)
}
dev.off()

# TODOs ------------------------------------------------------------------------

# - [ ] TREAT?
