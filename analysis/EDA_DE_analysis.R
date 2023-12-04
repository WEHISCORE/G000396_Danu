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
library(splines)

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

# Heatmaps (perhaps only showing samples involved in contrast?)
# TODO

# Gene set analysis
# TODO: Competitive vs. self-contained test
#       https://support.bioconductor.org/p/107553/.

gaf <- GOSemSim::read.gaf(
  here("data", "annotation", "PlasmoDB-66_Pfalciparum3D7_GO.gaf.gzip"))
idx <- ids2indices(
  split(gaf$TERM2GENE$Gene, gaf$TERM2GENE$GO),
  id = unlist(cfit$genes$GENEID))
dir.create(here("output", "CameraGeneSetTests"))
lapply(colnames(cfit), function(j) {
  message(j)
  campr <- cameraPR(statistic = cfit$t[, j], index = idx)
  campr$GO <- rownames(campr)
  campr <- dplyr::left_join(campr, gaf$TERM2NAME, by = c("GO" = "GOID"))
  rownames(campr) <- campr$GO
  campr <- campr[, c("TERM", "NGenes", "Direction", "PValue", "FDR")]
  write.csv(
    campr,
    here("output", "CameraGeneSetTests", paste0(j, ".camera.csv")))
})

# NOTE: camera() can't incorporate random effect
#       (https://support.bioconductor.org/p/78299/)
cam <- camera(cfit$EList, idx, design, contrast = cm[, j])
cam$GO <- rownames(cam)
cam <- dplyr::left_join(cam, gaf$TERM2NAME, by = c("GO" = "GOID"))
rownames(cam) <- cam$GO
cam <- cam[, c("TERM", "NGenes", "Direction", "PValue", "FDR")]
# NOTE: cameraPR() gets around this limitation of camera()
campr <- cameraPR(cfit$t[, j], idx)
campr$GO <- rownames(campr)
campr <- dplyr::left_join(campr, gaf$TERM2NAME, by = c("GO" = "GOID"))
rownames(campr) <- campr$GO
campr <- campr[, c("TERM", "NGenes", "Direction", "PValue", "FDR")]
# NOTE: https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#gene-set-testing-with-camera suggests using
#       camera() rather than mroast() [or fry()] when doing 'fishing' analyses.
# f <- fry(cfit$EList, idx, design, contrast = cm[, j])
# f$GO <- rownames(f)
# f <- dplyr::left_join(f, gaf$TERM2NAME, by = c("GO" = "GOID"))
# rownames(f) <- f$GO
# f <- f[
#   ,
#   c("TERM", "NGenes", "Direction", "PValue", "FDR", "PValue.Mixed",
#     "FDR.Mixed")]
# NOTE: romer() can't work with voom.
# rom <- romer(
#   estimateDisp(y),
#   idx,
#   design,
#   cm[, j],
#   block = y$samples$cell_line_rep,
#   correlation = 0.02)

# Multi-level DE analysis: Aggregated KOs vs. WT -------------------------------

g <- factor(
  sub("GID[1|2|7|8|9]KO", "KO", y$samples$group),
  c(
    "KO.Day_3", "KO.Day_6", "KO.Day_9", "KO.Day_12",
    "WT.Day_3", "WT.Day_6", "WT.Day_9", "WT.Day_12"))
design_kos <- model.matrix(~0 + g)
colnames(design_kos) <- sub("g", "", colnames(design))

fit_kos <- voomLmFit(
  y,
  design_kos,
  block = y$samples$cell_line_rep,
  sample.weights = TRUE,
  plot = TRUE,
  keep.EList = TRUE)
# NOTE: Very low correlation (0.02)
fit$correlation

cm <- makeContrasts(
  # Day6 vs. Day3 (within 'cell line')
  KO.Day_6_vs_KO.Day_3 = KO.Day_6 - KO.Day_3,
  WT.Day_6_vs_WT.Day_3 = WT.Day_6 - WT.Day_3,

  # Day9 vs. Day6 (within 'cell line')
  KO.Day_9_vs_KO.Day_6 = KO.Day_9 - KO.Day_6,
  WT.Day_9_vs_WT.Day_6 = WT.Day_9 - WT.Day_6,

  # Day12 vs. Day9 (within 'cell line')
  KO.Day_12_vs_KO.Day_9 = KO.Day_12 - KO.Day_9,
  WT.Day_12_vs_WT.Day_9 = WT.Day_12 - WT.Day_9,

  # Comparisons at Day 3 (cell lines vs. WT)
  KO.Day_3_vs_WT.Day_3 = KO.Day_3 - WT.Day_3,

  # Comparisons at Day 6 (cell lines vs. WT)
  KO.Day_6_vs_WT.Day_6 = KO.Day_6 - WT.Day_6,

  # Comparisons at Day 9 (cell lines vs. WT)
  KO.Day_9_vs_WT.Day_9 = KO.Day_9 - WT.Day_9,

  # Comparisons at Day 12 (cell lines vs. WT)
  KO.Day_12_vs_WT.Day_12 = KO.Day_12 - WT.Day_12,

  # Interactions (Day6 vs. Day3; cell lines vs. WT)
  `(KO.Day_6_vs_KO.Day_3)_vs_(WT.Day_6_vs_WT.Day_3)` =
    (KO.Day_6 - KO.Day_3) - (WT.Day_6 - WT.Day_3),

  # Interactions (Day9 vs. Day6; cell lines vs. WT)
  `(KO.Day_9_vs_KO.Day_6)_vs_(WT.Day_9_vs_WT.Day_6)` =
    (KO.Day_9 - KO.Day_6) - (WT.Day_9 - WT.Day_6),

  # Interactions (Day12 vs. Day9; cell lines vs. WT)
  `(KO.Day_12_vs_KO.Day_9)_vs_(WT.Day_12_vs_WT.Day_9)` =
    (KO.Day_12 - KO.Day_9) - (WT.Day_12 - WT.Day_9),

  levels = design)
cfit <- contrasts.fit(fit, cm)
cfit <- eBayes(cfit)
t(summary(decideTests(cfit)))

# TODO: Create outputs

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

# Time course analysis: Aggregated KOs vs. WT ----------------------------------

X <- ns(as.integer(y$samples$timepoint), df = 3)
Group <- factor(
  sub("GID[1|2|7|8|9]KO", "KO", y$samples$cell_line),
  c("WT", "KO"))
design_tc <- model.matrix(~Group * X)
colnames(design_tc) <- sub("Group", "", colnames(design_tc))

fit_tc <- voomLmFit(
  y,
  design_tc,
  block = y$samples$cell_line_rep,
  sample.weights = TRUE,
  plot = TRUE)
fit_tc <- eBayes(fit_tc)

tt_ko <- topTable(
  fit_tc,
  coef = c("KO:X1", "KO:X2", "KO:X3"),
  number = Inf,
  p.value = 0.05)
nrow(tt_ko)

# TODO: KOs_except_GID1 vs. WT

# TODO: Create outputs
