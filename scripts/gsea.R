# All proteins ranked by limma t + one gene set collection -> GSEA per contrast (positive NES = higher in A), plus a
# second run on |t| that finds sets whose proteins go both ways; figures
source("scripts/common.R")

samples    <- read_tsv(snakemake@input$samples, show_col_types = FALSE)
sets       <- read_tsv(snakemake@input$gene_sets, show_col_types = FALSE)
collection <- snakemake@wildcards$collection
pathways   <- split(sets$gene, sets$term)
ids        <- distinct(sets, term, id)
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

run_fgsea <- function(ranks, ...) {
  set.seed(1)  # fgsea estimates small p-values by sampling; a fixed seed makes reruns identical
  fgsea::fgsea(pathways, ranks, minSize = MIN_SET_SIZE, maxSize = MAX_SET_SIZE, nPermSimple = 10000, ...) |>
    as_tibble()
}

# call: up / down = enriched by signed t; mixed = enriched only by |t|, with leading-edge proteins in both directions.
# n_up / n_down count the leading edge that made the call (signed for up/down, |t| otherwise)
gsea <- map(de, \(stats) {
  ranks <- rank_genes(stats)
  any <- run_fgsea(abs(ranks), scoreType = "pos") |>
    transmute(term = pathway, any_NES = NES, any_pvalue = pval, any_padj = padj, any_leading_edge = leadingEdge)
  run_fgsea(ranks) |>
    transmute(term = pathway, NES, pvalue = pval, padj, size, leading_edge = leadingEdge) |>
    left_join(any, by = "term") |>
    mutate(
      directional = coalesce(padj < ENRICH_PADJ, FALSE),
      edge   = if_else(directional, leading_edge, any_leading_edge),
      n_up   = map_int(edge, \(g) sum(ranks[g] > 0)),
      n_down = map_int(edge, \(g) sum(ranks[g] < 0)),
      call   = case_when(
        directional ~ if_else(NES > 0, "up", "down"),
        any_padj < ENRICH_PADJ & n_up > 0 & n_down > 0 ~ "mixed",
        .default = ""
      ),
      across(c(leading_edge, any_leading_edge), \(x) map_chr(x, str_flatten, ";"))
    ) |>
    select(-directional, -edge) |>
    left_join(ids, by = "term") |>
    relocate(id, .after = term) |>
    arrange(pvalue)
})

dir.create(out, recursive = TRUE)
iwalk(gsea, \(result, name) write_tsv(result, file.path(out, str_c(name, ".tsv")), na = ""))

summary <- tibble(
  contrast = names(gsea),
  label    = titles[contrast],
  up       = map_int(gsea, \(x) sum(x$call == "up")),
  down     = map_int(gsea, \(x) sum(x$call == "down")),
  mixed    = map_int(gsea, \(x) sum(x$call == "mixed")),
  tested   = map_int(gsea, nrow)
)
write_tsv(summary, file.path(out, "summary.tsv"))

iwalk(gsea, \(result, name) {
  enrichment_figure(
    file.path(out, "figures", name),
    sets = result |>
      filter(call != "") |>
      transmute(panel = call, label = term_label(term), padj = if_else(call == "mixed", any_padj, padj), n_up, n_down,
                note = if_else(call == "mixed", str_c(n_up, " up, ", n_down, " down"), sprintf("NES %.2f", NES))),
    samples, name, method = "gsea",
    subtitle = str_glue("GSEA, {collection}; enriched = padj < {ENRICH_PADJ}"),
    message = str_glue("No {collection} gene set is enriched (padj < {ENRICH_PADJ}) in this contrast.")
  )
})

# Overview: every set called in at least one contrast, across all contrasts; padj of the test behind the call
all_contrasts <- list_rbind(gsea, names_to = "contrast") |> mutate(padj = if_else(call == "mixed", any_padj, padj))
called_terms <- unique(all_contrasts$term[all_contrasts$call != ""])
if (length(gsea) >= 2 && length(called_terms)) {
  make_figure(
    file.path(out, "figures", "overview"), script = "gsea_overview",
    data = all_contrasts |>
      filter(term %in% called_terms) |>
      transmute(label = term_label(term), contrast = titles[contrast], NES, padj, call)
  )
}
