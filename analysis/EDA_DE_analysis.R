# DE analysis of mini-bulk data for G000396_Danu
# Peter Hickey
# 2023-11-22

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

# TODO: Other colours.
sample_colours <- setNames(
  palette.colors(nlevels(sce$group), "Alphabet"),
  levels(sce$group))
cell_line_colours <- setNames(
  palette.colors(nlevels(sce$cell_line), "Dark2"),
  levels(sce$cell_line))
timepoint_colour <- setNames(
  palette.colors(nlevels(sce$timepoint), "Set2"),
  levels(sce$timepoint))

# edgeR-QL using all replicates-------------------------------------------------

y <- SE2DGEList(sce)

# Sum technical replicates
y <- sumTechReps(y, y$samples$sample)

keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]

# TODO: TMM or TMMwsp?
y <- normLibSizes(y, method = "TMMwsp")
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

# TODO: Incorporate `biological_replicate` where appropriate (i.e. when
#       comparing across timepoints within a cell line).
# design <- model.matrix(~0 + group, y$samples)
# colnames(design) <- sub("group", "", colnames(design))
# design <- model.matrix(~cell_line_rep + timepoint, y$samples)
# colnames(design) <- sub("cell_line_rep|timepoint", "", colnames(design))
design <- model.matrix(~cell_line:biological_replicate + timepoint, y$samples)
colnames(design) <- sub(
  "cell_line|biological_replicate|timepoint",
  "",
  colnames(design))

y <- estimateDisp(y, design)
# TODO: BCV plot is really wonky if using UMI counts.
plotBCV(y)
# NOTE: Reasonably low dispersion (typical value for well-controlled
#       experiments in human is 0.4 and 0.1 for data on genetically identical
#       model organisms).
sqrt(y$common.dispersion)

fit <- glmQLFit(y, design)
# TODO: Bit of an odd uptick when using reads and a big uptick when using UMIs
plotQLDisp(fit)

# TODO: What contrasts?
contrasts <- makeContrasts(
  GID1KO.Day_3 - WT.Day_3,
  GID1KO.Day_6 - WT.Day_6,
  GID1KO.Day_9 - WT.Day_9,
  GID1KO.Day_12 - WT.Day_12,

  GID2KO.Day_3 - WT.Day_3,
  GID2KO.Day_6 - WT.Day_6,
  GID2KO.Day_9 - WT.Day_9,
  GID2KO.Day_12 - WT.Day_12,

  GID7KO.Day_3 - WT.Day_3,
  GID7KO.Day_6 - WT.Day_6,
  GID7KO.Day_9 - WT.Day_9,
  GID7KO.Day_12 - WT.Day_12,

  GID9KO.Day_3 - WT.Day_3,
  GID9KO.Day_6 - WT.Day_6,
  GID9KO.Day_9 - WT.Day_9,
  GID9KO.Day_12 - WT.Day_12,

  GID1KO.Day_6 - GID1KO.Day_3,
  GID1KO.Day_9 - GID1KO.Day_6,
  GID1KO.Day_12 - GID1KO.Day_9,

  GID2KO.Day_6 - GID2KO.Day_3,
  GID2KO.Day_9 - GID2KO.Day_6,
  GID2KO.Day_12 - GID2KO.Day_9,

  GID7KO.Day_6 - GID7KO.Day_3,
  GID7KO.Day_9 - GID7KO.Day_6,
  GID7KO.Day_12 - GID7KO.Day_9,

  GID9KO.Day_6 - GID9KO.Day_3,
  GID9KO.Day_9 - GID9KO.Day_6,
  GID9KO.Day_12 - GID9KO.Day_9,

  WT.Day_6 - WT.Day_3,
  WT.Day_9 - WT.Day_6,
  WT.Day_12 - WT.Day_9,

  levels = design)

# TODO: Need to loop over all contrasts.
l_qlf <- lapply(colnames(contrasts), function(j) {
  message(j)
  qlf <- glmQLFTest(fit, contrast = contrasts[, j, drop = FALSE])
  qlf
})
names(l_qlf) <- colnames(contrasts)
t(sapply(l_qlf, function(qlf) summary(decideTests(qlf))))

# TODO: Would need to loop over elements of l_qlf to make the relevant plots.
topTags(l_qlf[[1]])
glimmaMA(l_qlf[[1]], y)

y0 <- y
fit0 <- fit
l_qlf0 <- l_qlf

# edgeR-QL filtering out low-quality replicates --------------------------------

y <- SE2DGEList(sce)

# Subset to relevant samples
# Using criteria from preprocessing Rmd.
# TODO: This is probably stricter than I want. Really just want to remove those
#       libraries with crap library sizes.
# TODO: Would need to adjust if using UMI counts because QC metrics in sce are
#       based on read counts.
libsize_drop <- isOutlier(
  metric = cd$sum,
  nmads = 3,
  type = "lower",
  log = TRUE)
feature_drop <- isOutlier(
  metric = cd$detected,
  nmads = 3,
  type = "lower",
  log = TRUE)
mito_drop <- isOutlier(
  metric = cd$subsets_Mito_percent,
  nmads = 3,
  type = "higher")
ribo_drop <- isOutlier(
  metric = cd$subsets_Ribo_percent,
  nmads = 3,
  type = "higher")
keep_rep <- !(libsize_drop | feature_drop | mito_drop | ribo_drop)
summary(keep_rep)
y <- y[, keep_rep]
y$samples <- droplevels(y$samples)

# Sum technical replicates
y <- sumTechReps(y, y$samples$sample)

keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]

# TODO: TMM or TMMwsp?
y <- normLibSizes(y, method = "TMMwsp")
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

# TODO: Using one-way layout. Is there any reason to use something different?
design <- model.matrix(~0 + group, y$samples)
colnames(design) <- sub("group", "", colnames(design))

y <- estimateDisp(y, design)
# TODO: BCV plot is really wonky if using UMI counts.
plotBCV(y)
# NOTE: Reasonably low dispersion (typical value for well-controlled
#       experiments in human is 0.4 and 0.1 for data on genetically identical
#       model organisms).
sqrt(y$common.dispersion)

fit <- glmQLFit(y, design)
# TODO: Bit of an odd uptick when using reads and a big uptick when using UMIs
plotQLDisp(fit)

# TODO: What contrasts?
contrasts <- makeContrasts(
  GID1KO.Day_3 - WT.Day_3,
  GID1KO.Day_6 - WT.Day_6,
  GID1KO.Day_9 - WT.Day_9,
  GID1KO.Day_12 - WT.Day_12,

  GID2KO.Day_3 - WT.Day_3,
  GID2KO.Day_6 - WT.Day_6,
  GID2KO.Day_9 - WT.Day_9,
  GID2KO.Day_12 - WT.Day_12,

  GID7KO.Day_3 - WT.Day_3,
  GID7KO.Day_6 - WT.Day_6,
  GID7KO.Day_9 - WT.Day_9,
  GID7KO.Day_12 - WT.Day_12,

  GID9KO.Day_3 - WT.Day_3,
  GID9KO.Day_6 - WT.Day_6,
  GID9KO.Day_9 - WT.Day_9,
  GID9KO.Day_12 - WT.Day_12,

  GID1KO.Day_6 - GID1KO.Day_3,
  GID1KO.Day_9 - GID1KO.Day_6,
  GID1KO.Day_12 - GID1KO.Day_9,

  GID2KO.Day_6 - GID2KO.Day_3,
  GID2KO.Day_9 - GID2KO.Day_6,
  GID2KO.Day_12 - GID2KO.Day_9,

  GID7KO.Day_6 - GID7KO.Day_3,
  GID7KO.Day_9 - GID7KO.Day_6,
  GID7KO.Day_12 - GID7KO.Day_9,

  GID9KO.Day_6 - GID9KO.Day_3,
  GID9KO.Day_9 - GID9KO.Day_6,
  GID9KO.Day_12 - GID9KO.Day_9,

  WT.Day_6 - WT.Day_3,
  WT.Day_9 - WT.Day_6,
  WT.Day_12 - WT.Day_9,

  levels = design)

# TODO: Need to loop over all contrasts.
l_qlf <- lapply(colnames(contrasts), function(j) {
  message(j)
  qlf <- glmQLFTest(fit, contrast = contrasts[, j, drop = FALSE])
  qlf
})
names(l_qlf) <- colnames(contrasts)
t(sapply(l_qlf, function(qlf) summary(decideTests(qlf))))

# TODO: Would need to loop over elements of l_qlf to make the relevant plots.
topTags(l_qlf[[1]])
glimmaMA(l_qlf[[1]], y)

y1 <- y
fit1 <- fit
l_qlf1 <- l_qlf

# voom using all replicates-----------------------------------------------------

# voom filtering out low-quality replicates ------------------------------------

# voomLmFit using all replicates -----------------------------------------------

# voomLmFit filtering out low-quality replicates -------------------------------

# voomLmFit with quality weights using all replicates --------------------------

# voomLmFit with quality weights filtering out low-quality replicates ----------

# Outputs to share with Danu ---------------------------------------------------

# NOTE: Temporarily using the results from the
#       'edgeR-QL filtering out low-quality replicates' analysis.

# Glimma
dir.create(here("tmp", "Glimma"))
lapply(names(l_qlf1), function(j) {
  message(j)
  glimmaMA(
    l_qlf[[j]],
    y1,
    html = here(
      "tmp",
      "Glimma",
      paste0(sub(" - ", "_vs_", j), ".html")))
})

# Heatmaps
# TODO: Should heatmaps use v$E or some combination of v$E and
#       voom$weights if using voomWithQualityWeights? (fit$EList$E and
#       fit$EList$weights if using voomLmFit).
lcpm <- edgeR::cpm(y1, log = TRUE)
dir.create(here("tmp", "pheatmap"))
# TODO: Might need to use the `filename` argument of pheatmap() because
#       otherwise the last heatmap text gets overprinted.
pdf(here("tmp", "pheatmap", "pheatmap.pdf"), height = 12, width = 12)
lapply(names(l_qlf1), function(j) {
  message(j)
  pheatmap::pheatmap(
    # TODO: Show all samples or just those involved in the contrast?
    lcpm[
      rownames(topTags(l_qlf1[[j]], n = 100, p.value = 0.05)),
      order(y1$samples$group)],
    scale = "row",
    color = colorRampPalette(c("blue","white","red"))(100),
    fontsize_row = 6,
    fontsize_col = 4,
    fontsize = 6,
    annotation_col = y1$samples[, c("cell_line", "timepoint")],
    main = sub(" - ", "_vs_", j),
    # TODO: colours
    # annotation_colors = list(group = phenotype_of_interest_colours),
    angle_col = 45,
    treeheight_row = 30,
    treeheight_col = 30,
    # TODO: Cluster or not?
    cluster_cols = FALSE)
})
dev.off()

# Alternative heatmaps
dir.create(here("tmp", "coolmap"))
# TODO: Might need to use the `filename` argument of pheatmap() because
#       otherwise the last heatmap text gets overprinted.
pdf(here("tmp", "coolmap", "coolmap.pdf"), height = 12, width = 12)
lapply(names(l_qlf1), function(j) {
  message(j)
  coolmap(
    lcpm[
      rownames(topTags(l_qlf1[[j]], n = 100, p.value = 0.05)),
      order(y1$samples$group)],
    # TODO: Have to look at plot to figure out what these values should be; is
    #       there a better way?
    # colsep = c(4, 8),
    # TODO: Have to look at plot to figure out what these values should be; is
    #       there a better way?
    # margins = c(10, 5),
    srtCol = 45,
    cexCol = 0.8,
    # ColSideColors = phenotype_of_interest_colours[y10$samples$group],
    main = sub(" - ", "_vs_", j))
})
dev.off()
