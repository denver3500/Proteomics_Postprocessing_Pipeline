# Significant proteins of one criterion + one gene set collection -> over-representation of up, down and all
# significant proteins per contrast, figures
source("scripts/common.R")

criterion  <- snakemake@params$criterion
samples    <- read_tsv(snakemake@input$samples, show_col_types = FALSE)
sets       <- read_tsv(snakemake@input$gene_sets, show_col_types = FALSE)
collection <- snakemake@wildcards$collection
pathways   <- split(sets$gene, sets$term)
ids        <- distinct(sets, term, id)
universe   <- unique(sets$gene)  # detected genes annotated in the collection; unannotated genes cannot be enriched
titles     <- contrast_labels(samples)
de         <- read_de(snakemake@input$stats, samples) |> map(\(x) call_significance(x, criterion))
out        <- snakemake@output[[1]]

# Significant genes per direction; "any" (both together) catches sets whose proteins go both ways
significant_genes <- function(stats) {
  genes <- map(c(up = "up", down = "down"), \(d) intersect(first_gene(stats$gene[stats$direction == d]), universe))
  c(genes, list(any = union(genes$up, genes$down)))
}

over_represented <- function(genes) {
  fgsea::fora(pathways, genes, universe, minSize = MIN_SET_SIZE, maxSize = MAX_SET_SIZE) |>
    as_tibble() |>
    filter(overlap > 0) |>
    transmute(term = pathway, overlap, set_size = size, fold_enrichment = foldEnrichment, pvalue = pval, padj,
              overlap_genes = overlapGenes)
}

genes <- map(de, significant_genes)
ora <- map(genes, \(lists) {
  map(lists, over_represented) |>
    list_rbind(names_to = "direction") |>
    mutate(
      n_up   = map_int(overlap_genes, \(g) sum(g %in% lists$up)),
      n_down = map_int(overlap_genes, \(g) sum(g %in% lists$down)),
      genes  = map_chr(overlap_genes, str_flatten, ";"),
      .keep  = "unused"
    ) |>
    left_join(ids, by = "term") |>
    relocate(id, .after = term) |>
    arrange(pvalue)
})

# Sets for the figure: enriched up or down, and mixed = enriched only with both directions, with proteins in both
enriched_sets <- function(result) {
  enriched <- filter(result, padj < ENRICH_PADJ)
  directional <- filter(enriched, direction != "any")
  mixed <- filter(enriched, direction == "any", !term %in% directional$term, n_up > 0, n_down > 0)
  bind_rows(directional, mutate(mixed, direction = "mixed")) |>
    transmute(
      panel = direction,
      label = term_label(term),
      padj, n_up, n_down,
      note = if_else(panel == "mixed", str_c(n_up, " up, ", n_down, " down"), str_c(overlap, " proteins")) |>
        str_c(sprintf(", %.1fx", fold_enrichment))
    )
}
enriched <- map(ora, enriched_sets)

dir.create(out, recursive = TRUE)
iwalk(ora, \(result, name) write_tsv(result, file.path(out, str_c(name, ".tsv")), na = ""))

summary <- tibble(
  contrast   = names(ora),
  label      = titles[contrast],
  up         = map_int(enriched, \(x) sum(x$panel == "up")),
  down       = map_int(enriched, \(x) sum(x$panel == "down")),
  mixed      = map_int(enriched, \(x) sum(x$panel == "mixed")),
  genes_up   = map_int(genes, \(x) length(x$up)),
  genes_down = map_int(genes, \(x) length(x$down))
)
write_tsv(summary, file.path(out, "summary.tsv"))

iwalk(enriched, \(sets, name) {
  counts <- filter(summary, contrast == name)
  enrichment_figure(
    file.path(out, "figures", name), sets, samples, name, method = "ora",
    subtitle = str_glue("ORA, {collection}; proteins with {describe_criterion(criterion)}; enriched = padj < {ENRICH_PADJ}"),
    message = str_glue("No enriched gene sets: {counts$genes_up} up and {counts$genes_down} down significant proteins ",
                       "are in {collection} sets.")
  )
})
