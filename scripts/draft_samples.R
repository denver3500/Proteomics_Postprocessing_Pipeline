# Guess a sample sheet from facility names: {date}_{operator}_{condition}[_{run}]_R{replicate}
source("scripts/common.R")

cols <- sample_columns(read_table(snakemake@input$table))

prefix <- cols[1]
while (!all(startsWith(cols, prefix))) prefix <- str_sub(prefix, 1, -2)
prefix <- coalesce(str_extract(prefix, "^.*_"), "")
name <- str_sub(cols, nchar(prefix) + 1)

draft <- tibble(
  column    = cols,
  condition = str_remove(name, "(_\\d+)?_R\\d+$"),
  label     = default_label(condition),  # readable name for plots; edit freely
  replicate = as.integer(str_match(name, "_R(\\d+)$")[, 2]),
  exclude   = NA  # fill in a reason (e.g. "failed QC") to drop a sample
) |>
  mutate(replicate = coalesce(replicate, row_number()), .by = condition) |>
  arrange(fct_inorder(condition), replicate)

singletons <- draft |> count(condition) |> filter(n == 1) |> pull(condition)
if (length(singletons)) {
  warning(
    length(singletons), " condition(s) with a single sample, e.g. ", singletons[1],
    ". Sample names probably contain subject IDs: fix the condition column before approving.",
    call. = FALSE
  )
}

write_tsv(draft, snakemake@output[[1]], na = "")
