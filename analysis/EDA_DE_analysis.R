# DE analysis of mini-bulk data for G000396_Danu
# Peter Hickey
# 2023-11-23

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
# TODO: Experiment with UMI counts
counts(sce) <- assay(sce, "read_counts")
# NOTE: Only retain relevant colData columns.
cd <- colData(sce)
colData(sce) <- cd[
  ,
  c("cell_line", "timepoint", "biological_replicate", "technical_replicate", "sample")]
sce$group <- interaction(
  sce$cell_line,
  sce$timepoint,
  lex.order = TRUE,
  drop = TRUE)
sce$cell_line_rep <- interaction(
  sce$cell_line,
  sce$biological_replicate,
  lex.order = TRUE,
  drop = TRUE)

sample_colours <- setNames(
  palette.colors(nlevels(sce$group), "Alphabet"),
  levels(sce$group))
cell_line_colours <- setNames(
  palette.colors(nlevels(sce$cell_line), "Dark2"),
  levels(sce$cell_line))
timepoint_colour <- setNames(
  palette.colors(nlevels(sce$timepoint), "Set2"),
  levels(sce$timepoint))

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
  scale_colour_manual(values = sample_colours)

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
  scale_colour_manual(values = sample_colours)

glimmaMDS(y)

design <- model.matrix(~0 + group, y$samples)
colnames(design) <- sub("group", "", colnames(design))

v <- voom(y, design, plot = TRUE)
cor <- duplicateCorrelation(v, design, block = v$targets$cell_line_rep)
# NOTE: Very low (0.01)
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
  scale_colour_manual(values = sample_colours)

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
fit1 <- fit

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
  scale_colour_manual(values = sample_colours)

glimmaMDS(y)

design <- model.matrix(~0 + group, y$samples)
colnames(design) <- sub("group", "", colnames(design))

v <- voomWithQualityWeights(y, design, plot = TRUE)
cor <- duplicateCorrelation(v, design, block = v$targets$cell_line_rep)
# NOTE: Very low (0.01)
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
  scale_colour_manual(values = sample_colours)

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
  scale_colour_manual(values = sample_colours)

glimmaMDS(y)

design <- model.matrix(~0 + group, y$samples)
colnames(design) <- sub("group", "", colnames(design))

fit <- voomLmFit(
  y,
  design,
  block = y$samples$cell_line_rep,
  plot = TRUE,
  keep.EList = TRUE)
# NOTE: Very low correlation (0.01)
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
  scale_colour_manual(values = sample_colours)

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
  scale_colour_manual(values = sample_colours)

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
# NOTE: Very low correlation (0.01)
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

# Outputs to share with Danu ---------------------------------------------------


# TODOs ------------------------------------------------------------------------

# - [ ] TREAT
# - [ ] Timecourse analysis
