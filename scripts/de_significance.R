# limma statistics + one significance criterion -> significant proteins per contrast, counts and DE figures
source("scripts/common.R")

criterion <- snakemake@params$criterion
samples   <- read_tsv(snakemake@input$samples, show_col_types = FALSE)
proteins  <- read_tsv(snakemake@input$intensities, show_col_types = FALSE)
pairs     <- condition_pairs(samples)
labels    <- condition_labels(samples)
titles    <- contrast_labels(samples)
de        <- read_de(snakemake@input$stats, samples) |> map(\(x) call_significance(x, criterion))
figures   <- snakemake@output$figures

dir.create(snakemake@output$significant)
iwalk(de, \(stats, name) {
  stats |>
    filter(direction != "ns") |>
    relocate(direction) |>
    write_tsv(file.path(snakemake@output$significant, str_c(name, ".tsv")), na = "")
})

counts <- count_significant(de) |> mutate(label = titles[contrast], .after = contrast)
counts |>
  mutate(criterion = describe_criterion(criterion)) |>
  write_tsv(snakemake@output$summary)

make_figure(file.path(figures, "counts"), data = select(counts, contrast, label, up, down))

# Significance starts at the largest p that still passes the p-value cutoffs (for padj this varies per contrast)
p_cutoff <- function(stats) {
  passing <- filter(stats, padj < (criterion$padj %||% Inf), pvalue < (criterion$pvalue %||% Inf))
  if (nrow(passing)) max(passing$pvalue) else NA
}

# One volcano folder per contrast: its proteins + a one-row table with the title and cutoff lines
iwalk(de, \(stats, name) {
  make_figure(
    file.path(figures, "volcano", name), script = "volcano",
    data = select(stats, gene, log2fc, pvalue, direction),
    contrast = tibble(
      contrast  = name,
      title     = titles[[name]],
      subtitle  = str_c(describe_criterion(criterion), "; up = higher in ", labels[[pairs[[name]][1]]]),
      p_cutoff  = p_cutoff(stats),
      fc_cutoff = na_if(criterion$log2fc %||% 0, 0)
    )
  )
})

make_figure(
  file.path(figures, "pvalues"),
  data = de |>
    map(\(x) select(x, protein_id, pvalue)) |>
    list_rbind(names_to = "contrast") |>
    mutate(contrast = titles[contrast]) |>
    pivot_wider(names_from = contrast, values_from = pvalue)
)

# Top hits: the significant proteins with the smallest p-value in any contrast
top_ids <- de |>
  map(\(x) filter(x, direction != "ns")) |>
  list_rbind() |>
  summarise(pvalue = min(pvalue), .by = protein_id) |>
  slice_min(pvalue, n = N_TOP_HITS, with_ties = FALSE) |>
  pull(protein_id)

if (length(top_ids) >= 2) {
  make_figure(
    file.path(figures, "top_hits"),
    data = proteins |>
      filter(protein_id %in% top_ids) |>
      select(gene, all_of(samples$sample)) |>
      mutate(gene = make.unique(gene)),
    samples = tibble(sample = samples$sample, condition = unname(labels[samples$condition]))
  )
}
