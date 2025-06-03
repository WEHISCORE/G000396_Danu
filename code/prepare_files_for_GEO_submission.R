# Prepare G000396_Danu data for GEO submission
# Peter Hickey
# 2025-06-03

library(here)
library(SingleCellExperiment)

outdir <- here("GEO")
dir.create(outdir, recursive = TRUE)

# FASTQs -----------------------------------------------------------------------

# NOTE: These plate-level FASTQ files are created by code/scPipe.R
dir.create(file.path(outdir, "FASTQ"))

# NOTE: There was only 1 RPI
file.copy(
  from = here("extdata", "scPipe", "RPI-43", "RPI-43.R1.fastq.gz"),
  to = file.path(outdir, "FASTQ"),
  recursive = FALSE,
  overwrite = FALSE)
file.copy(
  from = here("extdata", "scPipe", "RPI-43", "RPI-43.R2.fastq.gz"),
  to = file.path(outdir, "FASTQ"),
  recursive = FALSE,
  overwrite = FALSE)

# Mini-bulk SCE ----------------------------------------------------------------

sce <- readRDS(here("data", "SCEs", "G000396_Danu.preprocessed.SCE.rds"))
stopifnot(all(lengths(rowData(sce)$GENEID)) == 1)
rownames(sce) <- unlist(rowData(sce)$GENEID)
rowData(sce) <- S4Vectors::make_zero_col_DFrame(nrow(sce))
colnames(sce) <- paste0(sce$sample, ".", sce$technical_replicate)
cns <- c(
  # NOTE: Commented out columns that are technical factors rather than
  # "plate_number", "well_position", "sample_type", "original_plate",
  # "c_rt1_primer_name", "rd1_index_cell_index_index_sequence_as_in_c_rt1_primer",
  # "illumina_index_index_number_separate_index_read",
  # "illumina_index_index_sequence_separate_index_read",
  # "rt1_index_primer_sequences", "sequencing_run",
  "cell_line", "timepoint", "biological_replicate", "technical_replicate",
  "sample", "group", "cell_line_rep")
colData(sce) <- colData(sce)[, cns]

# Export data
# NOTE: There is only 1 non-zero ERCC count, so not exporting these.
dir.create(file.path(outdir, "SCE"))
# Gene counts
# NOTE: Exporting read counts.
write.csv(
  x = as.data.frame(as.matrix(assay(sce, "read_counts"))),
  file = gzfile(
    file.path(outdir, "SCE", "G000396_Danu.gene_counts.csv.gz")),
  row.names = TRUE)
# colData
write.csv(
  x = as.data.frame(colData(sce)),
  file = gzfile(
    file.path(outdir, "SCE", "G000396_Danu.sample_sheet.csv.gz")),
  row.names = TRUE)

# MD5sums ----------------------------------------------------------------------

tools::md5sum(
  c(
    file.path(outdir, "FASTQ", "RPI-43.R1.fastq.gz"),
    file.path(outdir, "FASTQ", "RPI-43.R2.fastq.gz"),
    file.path(outdir, "SCE", "G000396_Danu.gene_counts.csv.gz"),
    file.path(outdir, "SCE", "G000396_Danu.sample_sheet.csv.gz")))
# /vast/scratch/users/hickey/G000396_Danu/GEO/FASTQ/RPI-43.R1.fastq.gz
# "9a45c6edba03941da03e30e53f62dbf3"
# /vast/scratch/users/hickey/G000396_Danu/GEO/FASTQ/RPI-43.R2.fastq.gz
# "f5f02bec105c308b2a1ad56dfc09c68a"
# /vast/scratch/users/hickey/G000396_Danu/GEO/SCE/G000396_Danu.gene_counts.csv.gz
# "f633346fc05dc1dd16b8b3175f9a7ce8"
# /vast/scratch/users/hickey/G000396_Danu/GEO/SCE/G000396_Danu.sample_sheet.csv.gz
# "4ae27927d693d9ee8a16132bd8be3ba6"
