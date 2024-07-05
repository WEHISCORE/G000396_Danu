# EDA integrating proteomics data with transcriptomics data.
# Peter Hickey
# 2024-06-27

# Notes ------------------------------------------------------------------------

# From Danu's email:
#
# - The columns of interest for us are the D6.GID7_D6.IGP comparisons (columns CC-CE).
# - I would like to look at the WT(IGP2) vs GID7KO samples at the Day6 timepoint.

# Setup ------------------------------------------------------------------------

library(here)

# Load data --------------------------------------------------------------------

all_proteomics_df <- data.table::fread(
  here("data", "proteomics", "full_dataset_P4762.csv"),
  data.table = FALSE)

rna_seq_df <- read.csv(
  here("output", "DEGs", "GID7KO.Day_6_vs_WT.Day_6.DEGs.csv"),
  row.names = 1)

# Munge data -------------------------------------------------------------------

coi_all <- colnames(all_proteomics_df) %in%
  c("ID_2",
    "Mass",
    "Protein.Group",
    "Protein.Ids",
    "Genes",
    "First.Protein.Description")
coi_comparison <- grepl("D6.GID7_D6.IGP", colnames(all_proteomics_df))
coi <- coi_all | coi_comparison

proteomics_df <- all_proteomics_df[, coi]

# Plots ------------------------------------------------------------------------

nrow(rna_seq_df)
length(unique(rna_seq_df$GENEID))
nrow(proteomics_df)
length(unique(proteomics_df$Genes))
common <- intersect(rna_seq_df$GENEID, proteomics_df$Genes)
length(common)

table(Sig = rna_seq_df$adj.P.Val < 0.05, Common = rna_seq_df$GENEID %in% common)
table(
  Sig = proteomics_df$D6.GID7_D6.IGP_adj.P.Val < 0.05,
  Common = proteomics_df$Genes %in% common)

x <- rna_seq_df[match(common, rna_seq_df$GENEID), ]
y <- proteomics_df[match(common, proteomics_df$Genes), ]
col <- dplyr::case_when(
  x$adj.P.Val < 0.05 & y$D6.GID7_D6.IGP_adj.P.Val > 0.05 ~ "dodgerblue",
  x$adj.P.Val > 0.05 & y$D6.GID7_D6.IGP_adj.P.Val < 0.05 ~ "orange",
  x$adj.P.Val < 0.05 & y$D6.GID7_D6.IGP_adj.P.Val < 0.05 ~ "red",
  TRUE ~ "black")

# Some randomization to mitigate overplotting.
set.seed(666)
i <- sample(nrow(x))
pdf(here("tmp/EDA_proteomics_integration.pdf"), 6, 6)
plot(
  x[i, "logFC"],
  y[i, "D6.GID7_D6.IGP_logFC"],
  xlab = "RNA-seq: logFC",
  ylab = "Proteomics: logFC",
  main = "GID7KO.Day_6 vs. WT.Day_6",
  col = col,
  pch = 16,
  cex = 0.8)
legend(
  "bottomright",
  legend = c("Non-sig", "RNA", "Protein", "Both"),
  pch = 16,
  col = c("black", "dodgerblue", "orange", "red"))
dev.off()
cor(x$logFC, y$D6.GID7_D6.IGP_logFC)

library(limma)
barcodeplot(
  x$logFC,
  index = y$D6.GID7_D6.IGP_adj.P.Val < 0.05,
  gene.weights = y$D6.GID7_D6.IGP_logFC[y$D6.GID7_D6.IGP_adj.P.Val < 0.05])

barcodeplot(
  y$D6.GID7_D6.IGP_logFC,
  index = x$adj.P.Val < 0.05,
  gene.weights = x$logFC[x$adj.P.Val < 0.05])

library(ggplot2)
library(plotly)
library(cowplot)
p <- ggplot(
  data.frame(
    x = x$logFC,
    y = y$D6.GID7_D6.IGP_logFC,
    col = dplyr::case_when(
      x$adj.P.Val < 0.05 & y$D6.GID7_D6.IGP_adj.P.Val > 0.05 ~ "RNA",
      x$adj.P.Val > 0.05 & y$D6.GID7_D6.IGP_adj.P.Val < 0.05 ~ "Protein",
      x$adj.P.Val < 0.05 & y$D6.GID7_D6.IGP_adj.P.Val < 0.05 ~ "Both",
      TRUE ~ "Non-sig"),
    gene = rownames(x))[i, ]) +
  geom_point(aes(x = x, y = y, colour = col, label = gene)) +
  scale_colour_manual(
    values = c(
      "RNA" = "dodgerblue",
      "Protein" = "orange",
      "Both" = "red",
      "Non-sig" = "black"),
    name = "DE status") +
  theme_cowplot() +
  xlab("RNA-seq: logFC") +
  ylab("Proteomics: logFC") +
  ggtitle("GID7KO.Day_6 vs. WT.Day_6")
ggplotly(p)

# TODOs ------------------------------------------------------------------------

# - [ ] What are the 'corrected' columns in proteomics data (e.g., 'corrected_D6.GID7_D6.IGP_adj.P.Val')?
# - [ ] The 'Significance' column seems to be based on the 'corrected' adjusted P-values (but these also include NA's which I don't understand).
