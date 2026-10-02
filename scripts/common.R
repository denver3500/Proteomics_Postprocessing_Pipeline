# Shared helpers: read a facility table into a standard shape; contrasts and significance for DE; figure folders
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

# Readable condition name for plots and tables, used when the sample sheet has no label
default_label <- function(condition) str_replace_all(condition, "_", " ")

condition_labels <- function(samples) {
  distinct(samples, condition, label) |> deframe()
}

# Every pair of conditions, named A_vs_B. A is listed first in the sample sheet, so log2fc > 0 means higher in A
condition_pairs <- function(samples) {
  pairs <- combn(unique(samples$condition), 2, simplify = FALSE)
  set_names(pairs, map_chr(pairs, \(p) str_c(p[1], "_vs_", p[2])))
}

contrast_labels <- function(samples) {
  labels <- condition_labels(samples)
  map_chr(condition_pairs(samples), \(p) str_c(labels[[p[1]]], " vs ", labels[[p[2]]]))
}

read_de <- function(dir, samples) {
  names(condition_pairs(samples)) |>
    set_names() |>
    map(\(name) read_tsv(file.path(dir, str_c(name, ".tsv")), show_col_types = FALSE))
}

# Significant = below every p-value cutoff given and at least the given |log2fc|; a missing cutoff does not filter
call_significance <- function(stats, criterion) {
  passes <- with(stats, padj < (criterion$padj %||% Inf) & pvalue < (criterion$pvalue %||% Inf) &
    abs(log2fc) >= (criterion$log2fc %||% 0))
  mutate(stats, direction = case_when(passes & log2fc > 0 ~ "up", passes ~ "down", .default = "ns"))
}

describe_criterion <- function(criterion) {
  c(
    if (!is.null(criterion$padj)) str_glue("padj < {criterion$padj}"),
    if (!is.null(criterion$pvalue)) str_glue("p < {criterion$pvalue}"),
    if (!is.null(criterion$log2fc)) str_glue("|log2FC| >= {criterion$log2fc}")
  ) |>
    str_flatten(" & ")
}

count_significant <- function(de) {
  tibble(
    contrast = names(de),
    up       = map_int(de, \(x) sum(x$direction == "up")),
    down     = map_int(de, \(x) sum(x$direction == "down")),
    total    = up + down
  )
}

# Constants the reports quote. QC flags: a sample is flagged when it stands out from the rest of the dataset
OUTLIER_MADS     <- 3   # proteins, contaminants, replicate correlation: > 3 MADs worse than the median
MAX_MEDIAN_SHIFT <- 1   # median intensity: > 1 log2 (2-fold) away from the typical sample
N_TOP_HITS       <- 50  # proteins in the DE top-hits heatmap
MIN_SET_SIZE     <- 10  # enrichment: gene sets are tested when 10-500 of their genes were detected
MAX_SET_SIZE     <- 500
ENRICH_PADJ      <- 0.05  # a gene set counts as enriched below this padj

# Enrichment works on genes: a protein group counts as its first gene (that of the leading protein)
first_gene <- function(gene) str_remove(gene, ";.*")

# Species from UniProt entry names (S10A8_HUMAN), as msigdbr names it and its MSigDB (HS/MM)
detect_species <- function(proteins) {
  suffix <- str_match(proteins$protein_name %||% character(), "_([A-Z0-9]+)(;|$)")[, 2]
  common <- names(which.max(table(suffix)))
  switch(common %||% "",
    HUMAN = list(name = "Homo sapiens", db = "HS"),
    MOUSE = list(name = "Mus musculus", db = "MM"),
    stop("enrichment supports human and mouse; cannot tell the species from protein names (S10A8_HUMAN), found: ",
         common %||% "none", call. = FALSE)
  )
}

# MSigDB serving a collection (Snakefile COLLECTIONS entry) for this species: its own, or human mapped to orthologs
msigdb_source <- function(msigdb, species) {
  db <- if (is.null(msigdb[[species$db]])) "HS" else species$db
  list(
    db = db,
    code = msigdb[[db]],
    note = str_glue("MSigDB {msigdb[[db]]} (msigdbr {packageVersion('msigdbr')})",
                    if (db != species$db) ", human sets mapped to {species$name} orthologs" else "")
  )
}

# MSigDB set name -> readable label: HALLMARK_MTORC1_SIGNALING -> MTORC1 SIGNALING
term_label <- function(term) str_replace_all(str_remove(term, "^[A-Z]+_"), "_", " ")

# Figure folder, one per picture: tables (data.tsv, ...) + plot.R, a copy of scripts/figures/<script>.R, run in the
# folder to draw the picture. The folder works on its own, so it can be copied anywhere and restyled
make_figure <- function(folder, ..., script = basename(folder)) {
  dir.create(folder, recursive = TRUE)
  iwalk(list(...), \(table, file) {
    table |>
      mutate(across(where(is.double), \(x) signif(x, 6))) |>  # plenty for drawing; full statistics are in de/all_proteins
      write_tsv(file.path(folder, str_c(file, ".tsv")), na = "")
  })
  file.copy(file.path("scripts/figures", str_c(script, ".R")), file.path(folder, "plot.R"))
  source(file.path(folder, "plot.R"), local = new.env(), chdir = TRUE)
  invisible(folder)
}

# Markdown image of a figure folder's PNG for the reports, shown at 72 px per inch like plots drawn by knitr
figure_image <- function(...) {
  path <- normalizePath(list.files(file.path(...), "\\.png$", full.names = TRUE))
  info <- attr(png::readPNG(path, info = TRUE), "info")
  str_glue("![]({path}){{width={round(info$dim[1] * 72 / info$dpi[1])}px}}")
}
