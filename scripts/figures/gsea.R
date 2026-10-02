# GSEA of one contrast: normalized enrichment score (NES) of the most significant gene sets per direction.
# Faded bars are not enriched (padj above the cutoff in the subtitle).
# data.tsv: one row per gene set; contrast.tsv: title. Edit the settings, then run in this folder: Rscript plot.R
suppressPackageStartupMessages(library(tidyverse))

WIDTH     <- 8     # inches
HEIGHT    <- NA    # inches; NA = grows with the number of sets
BASE_SIZE <- 11    # font size
COLOURS   <- c(up = "#B2182B", down = "#2166AC")
N_SETS    <- 10    # most significant sets shown per direction
MAX_LABEL <- 60    # longer set names are shortened
DPI       <- 300   # PNG resolution

data <- read_tsv("data.tsv", show_col_types = FALSE)
info <- read_tsv("contrast.tsv", show_col_types = FALSE)

shown <- data |>
  mutate(direction = if_else(NES > 0, "up", "down")) |>
  slice_min(pvalue, n = N_SETS, by = direction, with_ties = FALSE) |>
  mutate(label = fct_reorder(str_trunc(label, MAX_LABEL), NES))

plot <- ggplot(shown, aes(NES, label, fill = direction, alpha = enriched)) +
  geom_col() +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = COLOURS, guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3), labels = c(`TRUE` = "enriched", `FALSE` = "not enriched"),
                     name = NULL) +
  labs(title = info$title, subtitle = str_wrap(info$subtitle, 90), x = "normalized enrichment score (NES)", y = NULL) +
  theme_grey(BASE_SIZE) +
  theme(legend.position = "bottom", legend.box = "vertical", plot.title.position = "plot")

height <- coalesce(HEIGHT, max(3, 0.22 * nrow(shown) + 1.8))
walk(c("svg", "png"), \(ext) ggsave(str_c("gsea_", info$contrast, ".", ext), plot, width = WIDTH, height = height, dpi = DPI))
