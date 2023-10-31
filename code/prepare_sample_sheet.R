# Prepare sample sheet (S000443) for use with scPipe.R
# Peter Hickey
# 2023-10-31

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
