# Significant proteins of one criterion + one gene set collection -> over-representation of up and down proteins
# per contrast, figures
source("scripts/common.R")

criterion  <- snakemake@params$criterion
samples    <- read_tsv(snakemake@input$samples, show_col_types = FALSE)
sets       <- read_tsv(snakemake@input$gene_sets, show_col_types = FALSE)
collection <- snakemake@wildcards$collection
pathways   <- split(sets$gene, sets$term)
ids        <- distinct(sets, term, id)
universe   <- unique(sets$gene)  # detected genes annotated in the collection; unannotated genes cannot be enriched
pairs      <- condition_pairs(samples)
labels     <- condition_labels(samples)
titles     <- contrast_labels(samples)
de         <- read_de(snakemake@input$stats, samples) |> map(\(x) call_significance(x, criterion))
out        <- snakemake@output[[1]]

direction_genes <- function(stats, direction) {
  intersect(first_gene(stats$gene[stats$direction == direction]), universe)
}

over_represented <- function(genes) {
  fgsea::fora(pathways, genes, universe, minSize = MIN_SET_SIZE, maxSize = MAX_SET_SIZE) |>
    as_tibble() |>
    filter(overlap > 0) |>
    transmute(term = pathway, overlap, set_size = size, fold_enrichment = foldEnrichment, pvalue = pval, padj,
              genes = map_chr(overlapGenes, str_flatten, ";"))
}

ora <- map(de, \(stats) {
  c(up = "up", down = "down") |>
    map(\(direction) over_represented(direction_genes(stats, direction))) |>
    list_rbind(names_to = "direction") |>
    left_join(ids, by = "term") |>
    relocate(id, .after = term) |>
    arrange(pvalue)
})

dir.create(out, recursive = TRUE)
iwalk(ora, \(result, name) write_tsv(result, file.path(out, str_c(name, ".tsv")), na = ""))

summary <- tibble(
  contrast   = names(ora),
  label      = titles[contrast],
  up         = map_int(ora, \(x) sum(x$padj < ENRICH_PADJ & x$direction == "up")),
  down       = map_int(ora, \(x) sum(x$padj < ENRICH_PADJ & x$direction == "down")),
  genes_up   = map_int(de, \(x) length(direction_genes(x, "up"))),
  genes_down = map_int(de, \(x) length(direction_genes(x, "down")))
)
write_tsv(summary, file.path(out, "summary.tsv"))

iwalk(ora, \(result, name) {
  make_figure(
    file.path(out, "figures", name), script = "ora",
    data = transmute(result, direction, label = term_label(term), overlap, fold_enrichment, pvalue, padj,
                     enriched = padj < ENRICH_PADJ),
    contrast = tibble(
      contrast = name,
      title    = titles[[name]],
      subtitle = str_glue("ORA, {collection}; proteins with {describe_criterion(criterion)}; ",
                          "enriched = padj < {ENRICH_PADJ}; up = higher in {labels[[pairs[[name]][1]]]}")
    )
  )
})
