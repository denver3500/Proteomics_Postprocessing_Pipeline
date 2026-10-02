# Clean intensities + sample sheet -> limma statistics for every pair of conditions (+ a model summary)
source("scripts/common.R")
suppressPackageStartupMessages(library(limma))

samples  <- read_tsv(snakemake@input$samples, show_col_types = FALSE)
proteins <- read_tsv(snakemake@input$intensities, show_col_types = FALSE)

conditions <- unique(samples$condition)
if (length(conditions) < 2) stop("differential expression needs at least two conditions", call. = FALSE)
pairs <- condition_pairs(samples)
log2i <- as.matrix(proteins[samples$sample])

# One coefficient per condition (its mean). Contrasts are built as numbers (+1 for A, -1 for B),
# so condition names never have to be valid R names
design <- model.matrix(~ 0 + factor(samples$condition, levels = conditions))
colnames(design) <- conditions
contrasts <- sapply(pairs, \(p) (conditions == p[1]) - (conditions == p[2]))
rownames(contrasts) <- conditions

# Samples sharing a block (e.g. the same mouse) are correlated; limma estimates one correlation for all proteins.
# A sample with an empty block is its own block
paired <- "block" %in% names(samples)
if (paired) {
  block <- coalesce(as.character(samples$block), samples$sample)
  correlation <- duplicateCorrelation(log2i, design, block = block)$consensus.correlation
  fit <- lmFit(log2i, design, block = block, correlation = correlation)
} else {
  fit <- lmFit(log2i, design)
}
fit <- eBayes(contrasts.fit(fit, contrasts), trend = TRUE, robust = TRUE)

means <- map(set_names(conditions), \(c) rowMeans(log2i[, samples$condition == c, drop = FALSE], na.rm = TRUE))

dir.create(snakemake@output$stats, recursive = TRUE)
iwalk(pairs, \(p, name) {
  stats <- topTable(fit, coef = name, number = Inf, sort.by = "none")
  proteins |>
    select(any_of(names(ANNOTATION))) |>
    mutate(
      log2fc = stats$logFC,
      "mean_{p[1]}" := means[[p[1]]],
      "mean_{p[2]}" := means[[p[2]]],
      ave_log2 = stats$AveExpr,
      t = stats$t,
      pvalue = stats$P.Value,
      padj = stats$adj.P.Val
    ) |>
    arrange(pvalue) |>
    write_tsv(file.path(snakemake@output$stats, str_c(name, ".tsv")), na = "")
})

tribble(
  ~item, ~value,
  "Method", "limma: one mean per condition, moderated t-test (eBayes with intensity trend, robust)",
  "Paired by", if (paired) str_glue("block column ({n_distinct(block)} blocks)") else "none",
  "Within-block correlation", if (paired) as.character(round(correlation, 3)) else "not used",
  "Proteins", as.character(nrow(proteins)),
  "Conditions", str_flatten_comma(conditions),
  "Contrasts", as.character(length(pairs))
) |>
  write_tsv(snakemake@output$model)
