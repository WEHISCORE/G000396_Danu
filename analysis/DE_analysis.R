# DE analysis of mini-bulk data for G000396_Danu
# Peter Hickey
# 2023-12-07

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

# Multi-level DE analysis ------------------------------------------------------

# NOTE: Using voomLmFit with quality weights after filtering out low-quality
#       replicates

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

dir.create(here("output", "Glimma"))
glimmaMDS(y, html = here("output", "Glimma", "overall.MDS.html"))

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

  levels = design)
cfit <- contrasts.fit(fit, cm)
cfit <- eBayes(cfit)

t(summary(decideTests(cfit)))

# TODO: Discuss with Gordon.
# TODO: UP TO HERE Ways to look for genes where the KOs are all (consistently)
#       different to the WTs
# Looking at genes that are different between each KO and WT at Day_3
coefs <- c("GID1KO.Day_3_vs_WT.Day_3", "GID2KO.Day_3_vs_WT.Day_3",
           "GID7KO.Day_3_vs_WT.Day_3", "GID8KO.Day_3_vs_WT.Day_3",
           "GID9KO.Day_3_vs_WT.Day_3")
# An F-test of whether any one of these contrasts is significant (i.e. an
# ANOVA-like test).
nrow(topTable(cfit, coef = coefs, n = Inf, p.value = 0.05))

# The next 2 give identical results because each contrast is `separate`-ly
summary(decideTests(cfit[, coefs], method = "separate"))
summary(decideTests(cfit, method = "separate"))[, coefs]

summary(decideTests(cfit[, coefs], method = "nestedF"))


classifyTestsF(cfit[, coefs], p.value = 0.05)

# Outputs of Multi-level DE analysis -------------------------------------------

# DEG lists as CSVs
dir.create(here("output", "DEGs"))
lapply(colnames(cfit), function(j) {
  message(j)
  tt <- topTable(cfit, coef = j, n = Inf, p.value = 1, sort.by = "P")
  write.csv(
    flattenDF(
      tt[, c("GENEID", "Name", "description", "logFC", "AveExpr", "t",
             "P.Value", "adj.P.Val", "B")]),
    here("output", "DEGs", paste0(j, ".DEGs.csv")))
})

# Glimma MA plots
dir.create(here("output", "Glimma"))
lapply(colnames(cfit), function(j) {
  message(j)
  glimmaMA(
    x = cfit,
    dge = y,
    coef = j,
    anno = y$genes[, c("GENEID", "Name", "description")],
    sample.cols = unname(group_colours[y$samples$group]),
    html = here(
      "output",
      "Glimma",
      paste0(j, ".html")),
    main = j)
})

# Heatmaps
lcpm <- edgeR::cpm(y, log = TRUE)
dir.create(here("output", "heatmaps"))
lapply(colnames(cfit), function(j) {
  message(j)

  # 1. Samples ordered by `group`.
  pheatmap(
    lcpm[
      rownames(topTable(cfit, coef = j, n = 100, p.value = 0.05)),
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
    lcpm[rownames(topTable(cfit, coef = j, n = 100, p.value = 0.05)), ],
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
    lcpm[rownames(topTable(cfit, coef = j, n = 100, p.value = 0.05)), jj],
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

# NOTE: fry() is a self-contained gene set test.
lapply(colnames(cfit), function(j) {
  message(j)
  f <- fry(cfit$EList, idx, design, contrast = cm[, j])
  f$GO <- rownames(f)
  f <- dplyr::left_join(f, gaf$TERM2NAME, by = c("GO" = "GOID"))
  rownames(f) <- f$GO
  f <- f[
    ,
    c("TERM", "NGenes", "Direction", "PValue", "FDR", "PValue.Mixed",
      "FDR.Mixed")]
  write.csv(
    f,
    here("output", "DEGs", paste0(j, ".fry.csv")))
})

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

tt_gid1 <- topTable(
  fit_tc,
  coef = c("GID1KO:X1", "GID1KO:X2", "GID1KO:X3"),
  number = Inf,
  p.value = 0.05)
nrow(tt_gid1)
tt_gid2 <- topTable(
  fit_tc,
  coef = c("GID2KO:X1", "GID2KO:X2", "GID2KO:X3"),
  number = Inf,
  p.value = 0.05)
nrow(tt_gid2)
tt_gid7 <- topTable(
  fit_tc,
  coef = c("GID7KO:X1", "GID7KO:X2", "GID7KO:X3"),
  number = Inf,
  p.value = 0.05)
nrow(tt_gid7)
tt_gid8 <- topTable(
  fit_tc,
  coef = c("GID8KO:X1", "GID8KO:X2", "GID8KO:X3"),
  number = Inf,
  p.value = 0.05)
nrow(tt_gid8)
tt_gid9 <- topTable(
  fit_tc,
  coef = c("GID9KO:X1", "GID9KO:X2", "GID9KO:X3"),
  number = Inf,
  p.value = 0.05)
nrow(tt_gid9)

# TODO: KOs vs. WT?
# NOTE: This is an ANOVA-like test.
tt_kos <-  topTable(
  fit_tc,
  coef = c(
    "GID1KO:X1", "GID1KO:X2", "GID1KO:X3",
    "GID2KO:X1", "GID2KO:X2", "GID2KO:X3",
    "GID7KO:X1", "GID7KO:X2", "GID7KO:X3",
    "GID8KO:X1", "GID8KO:X2", "GID8KO:X3",
    "GID9KO:X1", "GID9KO:X2", "GID9KO:X3"),
  number = Inf,
  p.value = 0.05)
nrow(tt_kos)

# TODO: Create outputs
