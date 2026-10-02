# GSEA across contrasts: NES of the gene sets enriched in at least one contrast (red = higher in the first condition).
# A dot marks the contrasts where the set is enriched. Sets with similar patterns are placed next to each other.
# data.tsv: one row per gene set and contrast. Edit the settings, then run in this folder: Rscript plot.R
suppressPackageStartupMessages(library(tidyverse))

WIDTH     <- NA    # inches; NA = grows with the number of contrasts
HEIGHT    <- NA    # inches; NA = grows with the number of sets
BASE_SIZE <- 10    # font size
GRADIENT  <- c("#2166AC", "white", "#B2182B")  # negative, 0, positive NES
N_SETS    <- 30    # sets with the smallest padj in any contrast
MAX_LABEL <- 60    # longer set names are shortened
DPI       <- 300   # PNG resolution

data <- read_tsv("data.tsv", show_col_types = FALSE)

top <- data |>
  summarise(padj = min(padj, na.rm = TRUE), .by = label) |>
  slice_min(padj, n = N_SETS, with_ties = FALSE) |>
  pull(label)
shown <- filter(data, label %in% top)

nes <- shown |>
  select(label, contrast, NES) |>
  pivot_wider(names_from = contrast, values_from = NES, values_fill = 0) |>
  column_to_rownames("label") |>
  as.matrix()
rows <- if (nrow(nes) > 1) rownames(nes)[hclust(dist(nes))$order] else rownames(nes)
shown <- mutate(shown, label = factor(label, rows), contrast = fct_inorder(contrast))
limit <- max(abs(shown$NES), na.rm = TRUE)

plot <- ggplot(shown, aes(contrast, label, fill = NES)) +
  geom_tile(colour = "white") +
  geom_point(data = filter(shown, enriched), size = 0.8) +
  scale_fill_gradient2(low = GRADIENT[1], mid = GRADIENT[2], high = GRADIENT[3], limits = c(-limit, limit)) +
  scale_y_discrete(labels = \(label) str_trunc(label, MAX_LABEL)) +
  labs(x = NULL, y = NULL, caption = "dot = enriched in that contrast") +
  theme_minimal(BASE_SIZE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

width  <- coalesce(WIDTH, max(6, 0.3 * nlevels(shown$contrast) + 5))
height <- coalesce(HEIGHT, max(3, 0.18 * length(rows) + 2.5))
walk(c("svg", "png"), \(ext) ggsave(str_c("overview.", ext), plot, width = width, height = height, dpi = DPI))
