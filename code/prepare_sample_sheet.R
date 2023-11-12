# Prepare sample sheet (S000443) for use with scPipe.R.
# Peter Hickey
# 2023-11-13

# Setup ------------------------------------------------------------------------

library(here)
library(readxl)
library(dplyr)
library(janitor)

# Construct S000443/S000448 sample sheet ---------------------------------------

# NOTE: Same library was sequenced over 2 runs:
#       S000443: MiSeq Nano
#       S000448: NextSeq

file_S000443_S000448 <- here(
  "data",
  "sample_sheets",
  "G000396_Danu_MB_S000443_SeqprimerSept23.xlsx")

# NOTE: Header row is split across 2 lines, which I combine into 1 before
#       reading in the rest of the spreadsheet.
header_row <- read_excel(
  path = file_S000443_S000448,
  sheet = "Samples",
  skip = 2,
  n_max = 1)

# NOTE: No FACS data
header_row <- paste0(colnames(header_row), header_row[1, ])
header_row <- gsub("^\\.\\.\\.[0-9]+", "", header_row)
sample_sheet <- read_excel(
  path = file_S000443_S000448,
  sheet = "Samples",
  # NOTE: Manually specifying range because there are some elements in column J
  #       that shouldn't be there or included.
  range = "A5:K388",
  col_names = header_row,
  # NOTE: Setting the max guess_max value avoids problems with incorrectly
  #       guessed columns
  #       (https://github.com/tidyverse/readxl/issues/414#issuecomment-352437730)
  guess_max = 1048576)

# Tidy up names and empty rows/columns.
sample_sheet <- clean_names(sample_sheet)
sample_sheet <- remove_empty(
  sample_sheet,
  which = c("rows", "cols"))

# Filter out those wells without a cell index sequence.
# NOTE: There are some empty wells, specifically those with
#       - `sample_type` == "empty"
#       - `original_plate` == "N/A"
#       - `day` == "N/A",
#       - `sample` == "N/A"
#       - `technical_replicate` == "N/A"
#       - `illumina_index_index_number_separate_index_read` == "N/A"
#       - `illumina_index_index_sequence_separate_index_read` == "N/A"
sample_sheet <- filter(
  sample_sheet,
  !rd1_index_cell_index_index_sequence_as_in_c_rt1_primer %in%
    c("N/A", "removed"))

# Some final tidying.
sample_sheet <- mutate(
  sample_sheet,
  plate_number = "Danu_MB_1",
  # NOTE: There are some wonky well_positions (e.g., 'I19=A1') that need to
  #       be fixed (these occur because it means well I19 with primer A1,
  #       in SCORE's terminology. I've asked for this to be avoided going
  #       forward.).
  well_position = gsub(" ", "", well_position),
  well_position = sapply(strsplit(well_position, "="), "[[", 1),
  well_position = factor(
    x = well_position,
    levels = unlist(
      lapply(LETTERS[1:16], function(x) paste0(x, 1:24)),
      use.names = TRUE)),
  # NOTE: Making everything RPI-43.
  illumina_index_index_number_separate_index_read = "RPI-43",
  sequencing_run = "S000443_S000448") |>
  arrange(plate_number, well_position) |>
  select(plate_number, everything())

# Re-naming of samples ---------------------------------------------------------

# NOTE: The sample IDs in the original 96-well plate sample sheets are an
#       absolute mess. There are duplicate sample IDs and non-consecutive
#       numbering of biological replicates within each experimental condition.
#       I discussed this with Danu in an email on 2023-11-10. Ultimately, the
#       resolution was to manually identify problematic samples and manually
#       re-label them. Sigh.

# NOTE: This data frame contains the samples that need modification/correction.
# NOTE: This data frame contains the samples that need modification/correction.
to_modify_df <- data.frame(
  # This is the well position in the 384-well CEL-Seq2 plate.
  well_position = c(
    "D13", "D15", "D17", "D19", "D21",
    "E13", "E15", "E17", "E19", "E21",
    "F1", "F3", "F5", "F7", "F9",
    "G1", "G3", "G5", "G7", "G9",
    "F13", "F15", "F17", "F19", "F21",
    "G13", "G15", "G17", "G19", "G21",
    "L13", "L15", "L17", "L19", "L21",
    "M13", "M15", "M17", "M19", "M21",
    "N1", "N3", "N5", "N7", "N9",
    "O1", "O3", "O5", "O7", "O9",
    "N13", "N15", "N17", "N19", "N21",
    "O13", "O15", "O17", "O19", "O21",
    "D14", "D16", "D18", "D20", "D22",
    "E14", "E16", "E18", "E20", "E22",
    "F2", "F4", "F6", "F8", "F10",
    "G2", "G4", "G6", "G8", "G10",
    "F14", "F16", "F18", "F20", "F22",
    "G14", "G16", "G18", "G20", "G22",
    "L14", "L16", "L18", "L20", "L22",
    "M14", "M16", "M18", "M20", "M22",
    "N2", "N4", "N6", "N8", "N10",
    "O2", "O4", "O6", "O8", "O10",
    "N14", "N16", "N18", "N20", "N22",
    "O14", "O16", "O18", "O20", "O22"),
  # This is the old sample ID used in the original 96-well plates.
  sample = c(
    rep(
      c("GID7KO_BioRep 2", "GID7KO_BioRep 3", "GID7KO_BioRep 4", "GID7KO_BioRep 5", "GID7KO_BioRep 6"),
      2),
    rep(
      c("GID8KO_BioRep 2", "GID8KO_BioRep 3", "GID8KO_BioRep 4", "GID8KO_BioRep 5", "GID8KO_BioRep 6"),
      2),
    rep(
      c("GID9KO_BioRep 3", "GID9KO_BioRep 4", "GID9KO_BioRep 5", "GID9KO_BioRep 6", "GID9KO_BioRep 7"),
      2),
    rep(
      c("GID7KO_BioRep 2", "GID7KO_BioRep 3", "GID7KO_BioRep 4", "GID7KO_BioRep 5", "GID7KO_BioRep 6"),
      2),
    rep(
      c("GID8KO_BioRep 2", "GID8KO_BioRep 3", "GID8KO_BioRep 4", "GID8KO_BioRep 5", "GID8KO_BioRep 6"),
      2),
    rep(
      c("GID9KO_BioRep 3", "GID9KO_BioRep 4", "GID9KO_BioRep 5", "GID9KO_BioRep 6", "GID9KO_BioRep 7"),
      2),
    rep(
      c("GID7KO_BioRep 2", "GID7KO_BioRep 9", "GID7KO_BioRep 3", "GID7KO_BioRep 5", "GID7KO_BioRep 6"),
      2),
    rep(
      c("GID8KO_BioRep 2", "GID8KO_BioRep 9", "GID8KO_BioRep 3", "GID8KO_BioRep 5", "GID8KO_BioRep 6"),
      2),
    rep(
      c("GID9KO_BioRep 9", "GID9KO_BioRep 4", "GID9KO_BioRep 3", "GID9KO_BioRep 6", "GID9KO_BioRep 7"),
      2),
    rep(
      c("GID7KO_BioRep 2", "GID7KO_BioRep 12", "GID7KO_BioRep 3", "GID7KO_BioRep 5", "GID7KO_BioRep 12"),
      2),
    rep(
      c("GID8KO_BioRep 2", "GID8KO_BioRep 12", "GID8KO_BioRep 3", "GID8KO_BioRep 5", "GID8KO_BioRep 12"),
      2),
    rep(
      c("GID9KO_BioRep 12", "GID9KO_BioRep 4", "GID9KO_BioRep 3", "GID9KO_BioRep 12", "GID9KO_BioRep 7"),
      2)),
  # This is the corrected sample ID for the original 96-well plate layout.
  corrected_sample = c(
    rep(
      c("GID7KO_BioRep 1", "GID7KO_BioRep 2", "GID7KO_BioRep 3", "GID7KO_BioRep 4", "GID7KO_BioRep 5"),
      2),
    rep(
      c("GID8KO_BioRep 1", "GID8KO_BioRep 2", "GID8KO_BioRep 3", "GID8KO_BioRep 4", "GID8KO_BioRep 5"),
      2),
    rep(
      c("GID9KO_BioRep 1", "GID9KO_BioRep 2", "GID9KO_BioRep 3", "GID9KO_BioRep 4", "GID9KO_BioRep 5"),
      2),
    rep(
      c("GID7KO_BioRep 1", "GID7KO_BioRep 2", "GID7KO_BioRep 3", "GID7KO_BioRep 4", "GID7KO_BioRep 5"),
      2),
    rep(
      c("GID8KO_BioRep 1", "GID8KO_BioRep 2", "GID8KO_BioRep 3", "GID8KO_BioRep 4", "GID8KO_BioRep 5"),
      2),
    rep(
      c("GID9KO_BioRep 1", "GID9KO_BioRep 2", "GID9KO_BioRep 3", "GID9KO_BioRep 4", "GID9KO_BioRep 5"),
      2),
    rep(
      c("GID7KO_BioRep 1", "GID7KO_BioRep 2", "GID7KO_BioRep 3", "GID7KO_BioRep 4", "GID7KO_BioRep 5"),
      2),
    rep(
      c("GID8KO_BioRep 1", "GID8KO_BioRep 2", "GID8KO_BioRep 3", "GID8KO_BioRep 4", "GID8KO_BioRep 5"),
      2),
    rep(
      c("GID9KO_BioRep 1", "GID9KO_BioRep 2", "GID9KO_BioRep 3", "GID9KO_BioRep 4", "GID9KO_BioRep 5"),
      2),
    rep(
      c("GID7KO_BioRep 1", "GID7KO_BioRep 2", "GID7KO_BioRep 3", "GID7KO_BioRep 4", "GID7KO_BioRep 5"),
      2),
    rep(
      c("GID8KO_BioRep 1", "GID8KO_BioRep 2", "GID8KO_BioRep 3", "GID8KO_BioRep 4", "GID8KO_BioRep 5"),
      2),
    rep(
      c("GID9KO_BioRep 1", "GID9KO_BioRep 2", "GID9KO_BioRep 3", "GID9KO_BioRep 4", "GID9KO_BioRep 5"),
      2)))

# Rename bad samples and recombine with good samples.
ok <- !sample_sheet$well_position %in% to_modify_df$well_position
sample_sheet <- bind_rows(
  left_join(sample_sheet[!ok, ], to_modify_df) |>
    mutate(
      sample = corrected_sample,
      corrected_sample = NULL),
  mutate(sample_sheet[ok, ])) |>
  mutate(
    well_position =
      factor(well_position, levels(sample_sheet$well_position))) |>
  arrange(well_position)

# Split sample ID into its components.
sample_sheet$cell_line <- sapply(strsplit(sample_sheet$sample, "_"), "[[", 1)
sample_sheet$timepoint <- sub(" ", "_", sample_sheet$day)
sample_sheet$day <- NULL
sample_sheet$biological_replicate <- sub(
  ".* ",
  "",
  as.character(sample_sheet$sample))
sample_sheet$sample <- NULL

sample_sheet <- select(
  sample_sheet,
  plate_number,
  well_position,
  sample_type,
  original_plate,
  cell_line,
  timepoint,
  biological_replicate,
  technical_replicate,
  everything())

# Construct final sample sheet -------------------------------------------------

sample_sheet <- sample_sheet |>
  mutate(rowname = paste0(plate_number, "_", well_position)) |>
  tibble::column_to_rownames("rowname")

# NOTE: Check that there aren't any malformed well_positions (e.g., 'I19=A1').
stopifnot(!anyNA(sample_sheet$well_position))

# Write sample sheet to disk ---------------------------------------------------

write.csv(
  sample_sheet,
  file = here("data", "sample_sheets", "S000443_S000448.sample_sheet.csv"))
