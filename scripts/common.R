# Shared helpers: read a facility table into a standard shape
suppressPackageStartupMessages(library(tidyverse))

# Raw annotation column -> clean name. Extend when a new table format arrives.
ANNOTATION <- c(
  protein_id    = "Protein.Group",
  gene          = "Genes",
  protein_name  = "Protein.Names",
  description   = "First.Protein.Description",
  n_sequences   = "N.Sequences",
  n_proteotypic = "N.Proteotypic.Sequences"
)
# Numeric columns that are not samples (facility statistics)
NOT_SAMPLES <- "^(FC|logp|logq)_|_Omni$"

read_table <- function(path) {
  read_tsv(path, show_col_types = FALSE) |>
    rename(any_of(c(Protein.Group = "Protein"))) |>
    rename(any_of(ANNOTATION))
}

sample_columns <- function(raw) {
  raw |>
    select(where(is.numeric), -any_of(names(ANNOTATION)), -matches(NOT_SAMPLES)) |>
    names()
}
