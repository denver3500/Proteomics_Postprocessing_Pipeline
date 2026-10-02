# All proteins ranked by limma t + one gene set collection -> GSEA per contrast (positive NES = higher in A), figures
source("scripts/common.R")

samples    <- read_tsv(snakemake@input$samples, show_col_types = FALSE)
sets       <- read_tsv(snakemake@input$gene_sets, show_col_types = FALSE)
collection <- snakemake@wildcards$collection
pathways   <- split(sets$gene, sets$term)
ids        <- distinct(sets, term, id)
pairs      <- condition_pairs(samples)
labels     <- condition_labels(samples)
titles     <- contrast_labels(samples)
de         <- read_de(snakemake@input$stats, samples)
out        <- snakemake@output[[1]]

# One value per gene: proteins sharing a first gene keep the strongest t
rank_genes <- function(stats) {
  stats |>
    filter(!is.na(t)) |>
    mutate(gene = first_gene(gene)) |>
    slice_max(abs(t), by = gene, n = 1, with_ties = FALSE) |>
    select(gene, t) |>
    deframe()
}

gsea <- map(de, \(stats) {
  set.seed(1)  # fgsea estimates small p-values by sampling; a fixed seed makes reruns identical
  fgsea::fgsea(pathways, rank_genes(stats), minSize = MIN_SET_SIZE, maxSize = MAX_SET_SIZE, nPermSimple = 10000) |>
    as_tibble() |>
    transmute(term = pathway, NES, ES, size, pvalue = pval, padj, leading_edge = map_chr(leadingEdge, str_flatten, ";")) |>
    left_join(ids, by = "term") |>
    relocate(id, .after = term) |>
    arrange(pvalue)
})

dir.create(out, recursive = TRUE)
iwalk(gsea, \(result, name) write_tsv(result, file.path(out, str_c(name, ".tsv")), na = ""))

summary <- tibble(
  contrast = names(gsea),
  label    = titles[contrast],
  up       = map_int(gsea, \(x) sum(x$padj < ENRICH_PADJ & x$NES > 0, na.rm = TRUE)),
  down     = map_int(gsea, \(x) sum(x$padj < ENRICH_PADJ & x$NES < 0, na.rm = TRUE)),
  tested   = map_int(gsea, nrow)
)
write_tsv(summary, file.path(out, "summary.tsv"))

iwalk(gsea, \(result, name) {
  make_figure(
    file.path(out, "figures", name), script = "gsea",
    data = result |>
      filter(!is.na(NES)) |>
      transmute(label = term_label(term), NES, pvalue, padj, enriched = padj < ENRICH_PADJ),
    contrast = tibble(
      contrast = name,
      title    = titles[[name]],
      subtitle = str_glue("GSEA, {collection}; enriched = padj < {ENRICH_PADJ}; positive NES = higher in {labels[[pairs[[name]][1]]]}")
    )
  )
})

# Overview: every set enriched in at least one contrast, across all contrasts
all_contrasts <- list_rbind(gsea, names_to = "contrast")
enriched_terms <- unique(all_contrasts$term[all_contrasts$padj < ENRICH_PADJ & !is.na(all_contrasts$padj)])
if (length(gsea) >= 2 && length(enriched_terms)) {
  make_figure(
    file.path(out, "figures", "overview"), script = "gsea_overview",
    data = all_contrasts |>
      filter(term %in% enriched_terms) |>
      transmute(label = term_label(term), contrast = titles[contrast], NES, padj, enriched = padj < ENRICH_PADJ)
  )
}
