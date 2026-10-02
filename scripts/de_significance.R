# limma statistics + one significance criterion -> significant proteins and a volcano plot per contrast, counts
source("scripts/common.R")

criterion <- snakemake@params$criterion
samples   <- read_tsv(snakemake@input$samples, show_col_types = FALSE)
pairs     <- condition_pairs(samples)
labels    <- condition_labels(samples)
titles    <- contrast_labels(samples)
de        <- read_de(snakemake@input$stats, samples) |> map(\(x) call_significance(x, criterion))

dir.create(snakemake@output$significant)
dir.create(snakemake@output$volcano)

iwalk(de, \(stats, name) {
  stats |>
    filter(direction != "ns") |>
    relocate(direction) |>
    write_tsv(file.path(snakemake@output$significant, str_c(name, ".tsv")), na = "")

  plot <- volcano_plot(stats, criterion, titles[[name]]) +
    labs(subtitle = str_glue("{describe_criterion(criterion)}; up = higher in {labels[[pairs[[name]][1]]]}"))
  ggsave(file.path(snakemake@output$volcano, str_c(name, ".pdf")), plot, width = 6, height = 5)
})

count_significant(de) |>
  mutate(label = titles[contrast], .after = contrast) |>
  mutate(criterion = describe_criterion(criterion)) |>
  write_tsv(snakemake@output$summary)
