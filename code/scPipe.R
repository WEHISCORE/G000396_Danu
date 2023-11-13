# Process S000443 (G000396) with scPipe
# Peter Hickey
# 2023-11-13

# Setup ------------------------------------------------------------------------

library(here)
library(scPipe)
library(Rsubread)
library(BiocParallel)

register(MulticoreParam(workers = 8))

# Load sample sheet ------------------------------------------------------------

sample_sheet <- as(
  read.csv(
    here("data", "sample_sheets", "S000443_S000448.sample_sheet.csv"),
    row.names = 1),
  "DataFrame")

# Key variables ----------------------------------------------------------------

rpis <- unique(sample_sheet$illumina_index_index_number_separate_index_read)
names(rpis) <- rpis
sequencing_runs <- tapply(
  sample_sheet$sequencing_run,
  sample_sheet$illumina_index_index_number_separate_index_read,
  unique)
outdir <- here("data", "SCEs")
dir.create(outdir, recursive = TRUE)
extdir <- file.path(
  "/vast/scratch/users/hickey/G000396_Danu/scPipe",
  rpis)
names(extdir) <- rpis
sapply(extdir, dir.create, recursive = TRUE)
# NOTE: Only using first 7 nt of barcode.
read_structure <- get_read_str("CEL-Seq2")
read_structure$bl2 <- 7
# NOTE: There's no plasmodium falciparum in this list. I've just used human
#       because these variables aren't important for running scPipe.
# NOTE: Must be an element of biomaRt::listDatasets(), e.g.,
#       biomaRt::listDatasets(biomaRt::useEnsembl("ensembl"))[["dataset"]]
organism <- "hsapiens_gene_ensembl"
# NOTE: Must be an element of biomaRt::listAttributes(), e.g.,
#       biomaRt::listAttributes(biomaRt::useEnsembl("ensembl", organism))[["name"]]
gene_id_type <- "ensembl_gene_id"

# FASTQ files
r1_fq <- c(
  grep(
    "Undetermined",
    list.files(
      path = here("extdata", "S000443", "data", "fastq"),
      full.names = TRUE,
      pattern = glob2rx("*R1*.fastq.gz")),
    value = TRUE,
    invert = TRUE),
  grep(
    "Undetermined",
    list.files(
      path = here("extdata", "S000448", "data", "fastq"),
      full.names = TRUE,
      pattern = glob2rx("*R1*.fastq.gz")),
    value = TRUE,
    invert = TRUE))
r2_fq <- gsub("R1", "R2", r1_fq)
stopifnot(all(file.exists(r2_fq)))
tx_fq <- file.path(extdir, paste0(rpis, ".R2.fastq.gz"))
names(tx_fq) <- rpis
barcode_fq <- gsub("R2", "R1", tx_fq)

bplapply(rpis, function(rpi) {
  message(rpi)
  # NOTE: S000443 has RPI in file name whereas G000448 does not.
  cmd <- paste0(
    "cat ",
    grep(rpi, r1_fq, value = TRUE),
    " ",
    grep(rpi, r1_fq, value = TRUE, invert = TRUE),
    " > ",
    barcode_fq[[rpi]],
    "\n",
    "cat ",
    grep(rpi, r2_fq, value = TRUE),
    " ",
    grep(rpi, r2_fq, value = TRUE, invert = TRUE),
    " > ",
    tx_fq[[rpi]])
  system(cmd)
})

# Genome index
genome_index <- here(
  "extdata",
  "PlasmoDB-51_Pfalciparum3D7",
  "PlasmoDB-51_Pfalciparum3D7_with_ERCC")

# Genome annotation(s)
annofn <- c(
  here(
    "extdata",
    "PlasmoDB-51_Pfalciparum3D7",
    "PlasmoDB-51_Pfalciparum3D7.gff"),
  system.file("extdata", "ERCC92_anno.gff3", package = "scPipe"))

# Cell barcodes
bc_anno <- file.path(extdir, paste0(rpis, ".barcode_annotation.csv"))
names(bc_anno) <- rpis

for (rpi in rpis) {
  message(rpi)
  tmp <- sample_sheet[
    sample_sheet$illumina_index_index_number_separate_index_read == rpi, ]
  barcode_df <- data.frame(
    cell_id = row.names(tmp),
    barcode = strtrim(
      tmp$rd1_index_cell_index_index_sequence_as_in_c_rt1_primer,
      7))
  stopifnot(!anyDuplicated(barcode_df$barcode))
  write.csv(
    x = barcode_df,
    file = bc_anno[[rpi]],
    quote = FALSE,
    row.names = FALSE)
}

# Output files -----------------------------------------------------------------

combined_fq <- file.path(extdir, gsub("R[12]", "combined", basename(tx_fq)))
names(combined_fq) <- names(tx_fq)
subread_bam <- gsub("fastq.gz", "subread.bam", combined_fq, fixed = TRUE)
exon_bam <- gsub("subread", "exon", subread_bam)

# FASTQ reformatting -----------------------------------------------------------

filter_settings <- list(rmlow = TRUE, rmN = FALSE, minq = 20, numbq = 2)
# NOTE: Have to loop over files because sc_trim_barcode() is not vectorised.
bplapply(seq_along(tx_fq), function(i) {
  message(combined_fq[i])
  sc_trim_barcode(
    outfq = combined_fq[i],
    r1 = tx_fq[i],
    r2 = barcode_fq[i],
    read_structure = read_structure,
    filter_settings = filter_settings)
})

# Aligning reads to a reference genome -----------------------------------------

align(
  index = genome_index,
  readfile1 = combined_fq,
  output_file = subread_bam,
  nthreads = 8)

# Assigning reads to annotated exons -------------------------------------------

# NOTE: Have to make data.frame in SAF format of PlasmoDB-51_Pfalciparum3D7.gff
#       because scPipe::sc_exon_mapping() doesn't know how to parse this file.
plasmo_gr <- rtracklayer::import(annofn[[1]])
plasmo_exon_gr <- plasmo_gr[plasmo_gr$type == "exon"]
# NOTE: My understanding is SAF is 1-based (see
#       https://www.biostars.org/p/228636/#9470789)
plasmo_saf <- data.frame(
  # NOTE: Some exons have multiple `Parent`s. Don't know why or what this
  #       means, but for now just take the first one.
  GeneID = sapply(plasmo_exon_gr$Parent, "[[", 1),
  Chr = seqnames(plasmo_exon_gr),
  Start = start(plasmo_exon_gr),
  End = end(plasmo_exon_gr),
  Strand = strand(plasmo_exon_gr))
ercc_gr <- rtracklayer::import(annofn[[2]])
ercc_exon_gr <- ercc_gr[ercc_gr$type == "exon"]
ercc_saf <- data.frame(
  # NOTE: Some exons have multiple `Parent`s. Don't know why or what this
  #       means, but for now just take the first one.
  GeneID = ercc_exon_gr$Name,
  Chr = seqnames(ercc_exon_gr),
  Start = start(ercc_exon_gr),
  End = end(ercc_exon_gr),
  Strand = strand(ercc_exon_gr))
saf <- rbind(plasmo_saf, ercc_saf)

bam_tags <- list(am = "YE", ge = "GE", bc = "BC", mb = "OX")
bc_len <- read_structure$bl1 + read_structure$bl2
barcode_vector <- ""
UMI_len <- read_structure$ul
stnd <- TRUE
fix_chr <- FALSE
bplapply(seq_along(subread_bam), function(i) {
  message(i)
  sc_exon_mapping(
    inbam = subread_bam[i],
    outbam = exon_bam[i],
    annofn = saf,
    bam_tags = bam_tags,
    bc_len = bc_len,
    barcode_vector = barcode_vector,
    UMI_len = UMI_len,
    stnd = stnd,
    fix_chr = fix_chr)
})

# De-multiplexing data ---------------------------------------------------------

max_mis <- 1
has_UMI <- TRUE
mito <- "Pf3D7_MIT_v3"
bplapply(seq_along(exon_bam), function(i) {
  message(i)
  sc_demultiplex(
    inbam = exon_bam[i],
    outdir = extdir[i],
    bc_anno = bc_anno[i],
    max_mis = max_mis,
    bam_tags = bam_tags,
    mito = mito,
    has_UMI = has_UMI)
})

# Gene counting deduped data ---------------------------------------------------

UMI_cor <- 1
gene_fl <- FALSE
bplapply(seq_along(bc_anno), function(i) {
  message(i)
  sc_gene_counting(
    outdir = extdir[i],
    bc_anno = bc_anno[i],
    UMI_cor = UMI_cor,
    gene_fl = gene_fl)
})

# Create and save deduped SingleCellExperiment ---------------------------------

list_of_sce <- lapply(rpis, function(rpi) {
  message(rpi)
  create_sce_by_dir(
    datadir = extdir[[rpi]],
    organism = organism,
    gene_id_type = gene_id_type,
    pheno_data = sample_sheet[
      sample_sheet$illumina_index_index_number_separate_index_read == rpi, ],
    # NOTE: Create the report separately for more fine-grained control.
    report = FALSE)
})
sce <- do.call(combineCols, c(unname(list_of_sce), fill = 0L, delayed = FALSE))
sce$UMI_deduped <- TRUE
# NOTE: Need to tidy up metadata to support use with scPipe::QC_metrics() on
#       the returned object. This assumes that all elements of `list_of_sce`
#       have the same metadata.
metadata(sce) <- list(
  scPipe = list(
    version = metadata(sce)[[1]][["version"]],
    QC_cols = metadata(sce)[[1]][["QC_cols"]],
    demultiplex_info = data.frame(
      status = metadata(sce)[[1]][["demultiplex_info"]][["status"]],
      count = rowSums(
        do.call(
          cbind,
          lapply(
            metadata(sce)[seq(1, length(metadata(sce)), 2)],
            function(x) {
              x[["demultiplex_info"]][["count"]]
            })))),
    UMI_dup_info = data.frame(
      duplication.number =
        metadata(sce)[[1]][["UMI_dup_info"]][["duplication.number"]],
      count = rowSums(
        do.call(
          cbind,
          lapply(
            metadata(sce)[seq(1, length(metadata(sce)), 2)],
            function(x) {
              x[["UMI_dup_info"]][["count"]]
            }))))),
  Biomart = metadata(sce)[[2]])

assay(sce, withDimnames = FALSE) <- as(
  assay(sce, withDimnames = FALSE),
  "dgCMatrix")
sce <- splitAltExps(
  sce,
  ifelse(grepl("^ERCC", rownames(sce)), "ERCC", "Endogenous"))

saveRDS(
  sce,
  file.path(outdir, "G000396_Danu.UMI_deduped.scPipe.SCE.rds"),
  compress = "xz")

# Create QC report of deduped data----------------------------------------------

library(readr)
library(plotly)
library(DT)
library(scran)
library(scater)
library(Rtsne)
dir.create(here("output", "scPipe"), recursive = TRUE)
# NOTE: Needs a fix for https://github.com/LuyiTian/scPipe/issues/100.
lapply(rpis, function(rpi) {
  message(rpi)
  create_report(
    sample_name = rpi,
    outdir = extdir[[rpi]],
    r1 = tx_fq[[rpi]],
    r2 = barcode_fq[[rpi]],
    outfq = combined_fq[[rpi]],
    read_structure = read_structure,
    filter_settings = filter_settings,
    align_bam = subread_bam[[rpi]],
    genome_index = genome_index,
    map_bam = exon_bam[[rpi]],
    exon_anno = annofn,
    stnd = stnd,
    fix_chr = fix_chr,
    barcode_anno = bc_anno[[rpi]],
    max_mis = max_mis,
    UMI_cor = UMI_cor,
    gene_fl = gene_fl,
    organism = organism,
    gene_id_type = gene_id_type)

  # NOTE: Copy the QC report to the repository.
  file.copy(
    from = file.path(extdir[[rpi]], "report.nb.html"),
    to = here(
      "output",
      "scPipe",
      paste0(rpi, ".scPipe_QC_report.nb.html")),
    overwrite = TRUE)
})

# Gene counting non-deduped data -----------------------------------------------

# NOTE: Need to use my own gene counting function because not using UMI
#       deduplication.
geneCountingNoUMIDedup <- function(outdir, bc_anno) {
  files <- list.files(file.path(outdir, "count"), full.names = TRUE)
  names(files) <- sub("\\.csv", "", basename(files))
  counts <- lapply(files, function(file) {
    message(basename(file))
    data.table::fread(file, select = 1)[, table(gene_id)]
  })
  genes <- Reduce(union, lapply(counts, names))
  x <- matrix(
    0L,
    nrow = length(genes),
    ncol = length(files),
    dimnames = list(genes, names(counts)))
  for (j in names(counts)) {
    xx <- counts[[j]]
    x[names(xx), j] <- xx
  }
  z <- cbind(
    data.frame(gene_id = rownames(x)),
    as.data.frame(x))
  data.table::fwrite(
    x = z,
    file = file.path(paste0(outdir, "_no_dedup"), "gene_count.csv"),
    row.names = FALSE,
    nThread = 1)
}

sapply(paste0(extdir, "_no_dedup"), dir.create)
bplapply(names(bc_anno), function(n) {
  message(n)
  geneCountingNoUMIDedup(
    outdir = extdir[[n]],
    bc_anno = bc_anno[[n]])
})

# Create and save non-deduped SingleCellExperiment -----------------------------

list_of_no_dedup_sce <- lapply(rpis, function(rpi) {
  x <- data.table::fread(
    file.path(paste0(extdir[[rpi]], "_no_dedup"), "gene_count.csv"))
  counts <- as.matrix(x[, -1])
  sce <- SingleCellExperiment(
    list(counts = counts),
    colData = sample_sheet[colnames(counts), ])
  rownames(sce) <- x$gene_id
  sce
})

no_dedup_sce <- do.call(
  combineCols,
  c(unname(list_of_no_dedup_sce), fill = 0L, delayed = FALSE))
no_dedup_sce$UMI_deduped <- FALSE
colnames(no_dedup_sce) <- paste0(colnames(no_dedup_sce), ".not_UMI_deduped")
assay(no_dedup_sce, withDimnames = FALSE) <- unname(
  as(assay(no_dedup_sce, withDimnames = FALSE), "dgCMatrix"))
no_dedup_sce <- splitAltExps(
  no_dedup_sce,
  ifelse(grepl("^ERCC", rownames(no_dedup_sce)), "ERCC", "Endogenous"))
# NOTE: Order genes and samples as in `sce`.
no_dedup_sce <- no_dedup_sce[
  rownames(sce),
  paste0(colnames(sce), ".not_UMI_deduped")]

saveRDS(
  no_dedup_sce,
  file.path(outdir, "G000396_Danu.not_UMI_deduped.scPipe.SCE.rds"),
  compress = "xz")

# Copy outputs -----------------------------------------------------------------

# NOTE: This is very fragile and will likely require modification for each
#       project.
cmd <- paste0(
  "rsync --verbose --human-readable --recursive --progress --archive \ ",
  unique(dirname(extdir)), "\ ",
  here("extdata"))
system(cmd)
