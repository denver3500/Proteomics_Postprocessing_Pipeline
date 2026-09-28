# Raw table + approved sample sheet -> clean wide intensity table + sample sheet
source("scripts/common.R")

raw <- read_table(snakemake@input$table)
sheet_path <- snakemake@input$sheet
sheet <- read_tsv(sheet_path, col_types = cols(.default = "c"))

# Verify the approved sheet against the table
stopifnot(
  "sample sheet needs columns: column, condition, replicate" =
    all(c("column", "condition", "replicate") %in% names(sheet)),
  "sample sheet lists columns missing from the table" = all(sheet$column %in% names(raw))
)
unlisted <- setdiff(sample_columns(raw), sheet$column)
if (length(unlisted)) {
  stop(
    length(unlisted), " sample(s) in the table are not in ", sheet_path, ": ",
    str_flatten_comma(head(unlisted, 3)), if (length(unlisted) > 3) ", ...",
    ". Add them (see results/", snakemake@wildcards$dataset, "/samples_draft.tsv) ",
    "or list them with a reason in the exclude column.",
    call. = FALSE
  )
}

samples <- sheet |>
  filter(if_all(any_of("exclude"), is.na)) |>  # keep rows without an exclude reason
  select(-any_of("exclude")) |>
  mutate(replicate = as.integer(replicate)) |>
  arrange(fct_inorder(condition), replicate) |>
  mutate(sample = str_c(condition, "_R", replicate), .before = 1) |>
  relocate(condition, replicate, .after = sample) |>
  relocate(column, .after = last_col())

stopifnot(
  "every sample needs a condition and an integer replicate" =
    !anyNA(samples$condition) && !anyNA(samples$replicate),
  "sample names (condition_R#) are not unique" = !anyDuplicated(samples$sample)
)

intensities <- raw |>
  select(any_of(names(ANNOTATION)), all_of(set_names(samples$column, samples$sample))) |>
  filter(!str_starts(protein_id, "cRAP-")) |>
  mutate(
    across(any_of(c("protein_id", "gene", "protein_name")), \(x) str_remove_all(x, ";cRAP-[^;]*")),
    gene = coalesce(gene, protein_name, protein_id)
  )

write_tsv(intensities, snakemake@output$intensities, na = "")
write_tsv(samples, snakemake@output$samples, na = "")
