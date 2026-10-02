# Clean intensities + sample sheet (+ raw table for contaminants) -> QC tables and figures
source("scripts/common.R")

# side = +1: high values are bad; side = -1: low values are bad
is_outlier <- function(x, side) side * (x - median(x, na.rm = TRUE)) > OUTLIER_MADS * mad(x, na.rm = TRUE)

samples  <- read_tsv(snakemake@input$samples, show_col_types = FALSE) |>
  mutate(condition = fct_inorder(label))  # readable labels (one per condition), in sample-sheet order
proteins <- read_tsv(snakemake@input$intensities, show_col_types = FALSE)
raw      <- read_table(snakemake@input$table)
sheet    <- read_tsv(snakemake@input$sheet, col_types = cols(.default = "c"))

log2i <- as.matrix(proteins[samples$sample])  # log2 intensities, proteins x samples

# Contaminants: share of each sample's total intensity that comes from cRAP proteins
crap <- str_starts(raw$protein_id, "cRAP-")
raw_linear <- 2^as.matrix(raw[samples$column])
contaminants <- colSums(raw_linear[crap, , drop = FALSE], na.rm = TRUE) / colSums(raw_linear, na.rm = TRUE)

# Sample-to-sample correlation; `others` blanks the diagonal so a sample is not compared to itself
cors <- cor(log2i, use = "pairwise.complete.obs")
others <- cors
diag(others) <- NA
replicate_cors <- ifelse(outer(samples$condition, samples$condition, "=="), others, NA)

# Extreme low values (below the dataset's lower boxplot whisker) usually mark imputed, i.e. undetected, proteins
quartiles <- quantile(log2i, c(0.25, 0.75), na.rm = TRUE)
low_cutoff <- quartiles[[1]] - 1.5 * diff(quartiles)

qc <- samples |>
  select(sample, condition, replicate) |>
  mutate(
    quantified       = colSums(!is.na(log2i)),
    missing_pct      = 100 * colMeans(is.na(log2i)),
    low_values_pct   = 100 * colMeans(log2i < low_cutoff, na.rm = TRUE),
    contaminants_pct = 100 * contaminants,
    median_shift     = apply(log2i - apply(log2i, 1, median, na.rm = TRUE), 2, median, na.rm = TRUE),
    cor_within       = apply(replicate_cors, 2, median, na.rm = TRUE),
    cor_all          = apply(others, 2, median, na.rm = TRUE),
    flags = str_c(
      if_else(is_outlier(quantified, -1), "few proteins; ", "", missing = ""),
      if_else(is_outlier(contaminants_pct, +1), "high contaminants; ", "", missing = ""),
      if_else(abs(median_shift) > MAX_MEDIAN_SHIFT, "intensity shift; ", "", missing = ""),
      if_else(is_outlier(cor_within, -1), "low replicate correlation; ", "", missing = "")
    ) |> str_remove("; $")
  )

# Per-condition CV of each protein across replicates
cv <- proteins |>
  select(protein_id, all_of(samples$sample)) |>
  pivot_longer(-protein_id, names_to = "sample", values_to = "log2") |>
  left_join(select(samples, sample, condition), by = "sample") |>
  summarise(cv = 100 * sd(2^log2, na.rm = TRUE) / mean(2^log2, na.rm = TRUE), .by = c(condition, protein_id))

conditions <- qc |>
  summarise(replicates = n(), median_cor_within = median(cor_within), .by = condition) |>
  left_join(summarise(cv, median_cv_pct = median(cv, na.rm = TRUE), .by = condition), by = "condition")

# Dataset summary
contaminant_names <- str_remove(if ("protein_name" %in% names(raw)) raw$protein_name else raw$protein_id, "^cRAP-")
top_contaminants <- contaminant_names[crap][order(-rowSums(raw_linear[crap, , drop = FALSE], na.rm = TRUE))]
excluded <- sheet |> filter(if_any(any_of("exclude"), \(x) !is.na(x)))
complete <- log2i[complete.cases(log2i), , drop = FALSE]

overview <- tribble(
  ~item, ~value,
  "Dataset", snakemake@wildcards$dataset,
  "Samples", str_glue("{nrow(samples)} in {n_distinct(samples$condition)} conditions"),
  "Proteins", as.character(nrow(proteins)),
  "Proteins in every sample", as.character(nrow(complete)),
  "Missing values", str_glue("{round(100 * mean(is.na(log2i)), 1)}%"),
  "Single-peptide proteins", if ("n_sequences" %in% names(proteins)) {
    str_glue("{round(100 * mean(proteins$n_sequences == 1, na.rm = TRUE), 1)}%")
  } else "not available",
  "Ambiguous protein groups", as.character(sum(str_detect(proteins$protein_id, ";"))),
  "Contaminants removed", str_glue("{sum(crap)} cRAP proteins (top: {str_flatten_comma(head(top_contaminants, 5))})"),
  "Excluded samples", if (nrow(excluded)) str_flatten_comma(str_glue("{excluded$column} ({excluded$exclude})")) else "none",
  "Flagged samples", str_glue("{sum(qc$flags != '')} of {nrow(qc)}")
)

# PCA on the proteins measured in every sample of x: PC1/PC2 per sample + % variance of each
pca_scores <- function(x) {
  pca <- prcomp(t(x[complete.cases(x), , drop = FALSE]))
  variance <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  tibble(sample = colnames(x), PC1 = pca$x[, 1], PC2 = pca$x[, 2], var1 = variance[1], var2 = variance[2])
}

write_tsv(overview, snakemake@output$overview)
write_tsv(qc, snakemake@output$sample_metrics, na = "")
write_tsv(conditions, snakemake@output$condition_metrics, na = "")

figures <- snakemake@output$figures
sample_conditions <- select(samples, sample, condition)

make_figure(file.path(figures, "intensities"), data = select(proteins, protein_id, all_of(samples$sample)), samples = sample_conditions)
make_figure(file.path(figures, "contaminants"), data = select(qc, sample, condition, contaminants_pct))
make_figure(file.path(figures, "correlation"), data = bind_cols(sample_conditions, as_tibble(cors)))
make_figure(file.path(figures, "cv"), data = pivot_wider(cv, names_from = condition, values_from = cv))

if (nrow(complete) >= 2) {
  make_figure(file.path(figures, "pca"), data = left_join(select(qc, sample, condition, replicate, flags), pca_scores(log2i), by = "sample"))
}

# Each pair of conditions gets its own PCA, so a group that dominates the PCA above cannot hide their differences
if (nlevels(samples$condition) >= 3) {
  pairs <- combn(levels(samples$condition), 2, simplify = FALSE) |>
    map(\(pair) {
      pca_scores(log2i[, samples$condition %in% pair, drop = FALSE]) |>
        mutate(pair = str_c(pair[1], " vs ", pair[2]), .before = 1)
    }) |>
    list_rbind() |>
    left_join(select(samples, sample, condition, replicate), by = "sample") |>
    relocate(condition, replicate, .after = sample)
  make_figure(file.path(figures, "pca_pairs"), data = pairs)
}
