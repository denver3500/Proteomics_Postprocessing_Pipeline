# Shared helpers: read a facility table into a standard shape; contrasts, significance and volcano plots for DE
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

DIRECTION_COLOURS <- c(up = "#B2182B", down = "#2166AC", ns = "grey75")

# Volcano of one contrast with up/down counts in the top corners; expects the direction column from call_significance()
volcano_plot <- function(stats, criterion, title, n_labels = 10) {
  labelled <- stats |> filter(direction != "ns") |> slice_min(pvalue, n = n_labels, with_ties = FALSE)
  fc <- criterion$log2fc %||% 0

  # Significance starts at the largest p that still passes the p-value cutoffs (for padj this varies per contrast)
  passing <- filter(stats, padj < (criterion$padj %||% Inf), pvalue < (criterion$pvalue %||% Inf))
  p_line  <- if (nrow(passing)) geom_hline(yintercept = -log10(max(passing$pvalue)), linetype = "dashed", colour = "grey40")
  fc_line <- if (fc > 0) geom_vline(xintercept = c(-fc, fc), linetype = "dashed", colour = "grey40")

  ggplot(arrange(stats, direction != "ns"), aes(log2fc, -log10(pvalue), colour = direction)) +
    geom_point(size = 0.8, alpha = 0.7) +
    p_line +
    fc_line +
    ggrepel::geom_text_repel(aes(label = gene), data = labelled, size = 2.5, max.overlaps = Inf, show.legend = FALSE) +
    annotate("label", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, label = str_glue("{sum(stats$direction == 'down')} down"),
             colour = DIRECTION_COLOURS[["down"]], fontface = "bold", label.size = 0) +
    annotate("label", x = Inf, y = Inf, hjust = 1.1, vjust = 1.3, label = str_glue("{sum(stats$direction == 'up')} up"),
             colour = DIRECTION_COLOURS[["up"]], fontface = "bold", label.size = 0) +
    scale_colour_manual(values = DIRECTION_COLOURS, breaks = c("up", "down"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(title = title, x = "log2 fold change", y = "-log10 p-value")
}
