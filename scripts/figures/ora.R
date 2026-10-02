# Over-representation in one contrast: the most significant gene sets among up and among down proteins.
# x = fold enrichment, dot size = significant proteins in the set; faded dots are not enriched (padj above the cutoff).
# data.tsv: one row per gene set and direction; contrast.tsv: title. Edit the settings, then run in this folder: Rscript plot.R
suppressPackageStartupMessages(library(tidyverse))

WIDTH     <- 9     # inches
HEIGHT    <- NA    # inches; NA = grows with the number of sets
BASE_SIZE <- 11    # font size
COLOURS   <- c(up = "#B2182B", down = "#2166AC")
N_SETS    <- 10    # most significant sets shown per direction
MAX_LABEL <- 60    # longer set names are shortened
DPI       <- 300   # PNG resolution

data <- read_tsv("data.tsv", show_col_types = FALSE, col_types = cols(direction = "c", label = "c"))
info <- read_tsv("contrast.tsv", show_col_types = FALSE)

# Each panel orders its own sets (most significant on top), so a set found in both directions gets one row in each
shown <- data |>
  slice_min(pvalue, n = N_SETS, by = direction, with_ties = FALSE) |>
  mutate(
    direction = factor(direction, c("up", "down")),
    row = fct_reorder(str_c(str_trunc(label, MAX_LABEL), "@", direction), -pvalue)
  )

plot <- ggplot(shown, aes(fold_enrichment, row, colour = direction, alpha = enriched, size = overlap)) +
  geom_point() +
  facet_wrap(vars(direction), ncol = 1, scales = "free_y", drop = TRUE) +
  scale_y_discrete(labels = \(row) str_remove(row, "@[^@]*$")) +
  scale_x_continuous(expand = expansion(mult = 0.08)) +
  scale_colour_manual(values = COLOURS, guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3), labels = c(`TRUE` = "enriched", `FALSE` = "not enriched"),
                     name = NULL, guide = guide_legend(override.aes = list(size = 3))) +
  scale_size_area(max_size = 5, name = "proteins") +
  labs(title = info$title, subtitle = str_wrap(info$subtitle, 90), x = "fold enrichment", y = NULL) +
  theme_grey(BASE_SIZE) +
  theme(legend.position = "bottom", legend.box = "vertical", plot.title.position = "plot")

if (!nrow(shown)) {
  plot <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No gene set shares a significant protein in this contrast") +
    labs(title = info$title, subtitle = str_wrap(info$subtitle, 90)) +
    theme_void(BASE_SIZE)
}

height <- coalesce(HEIGHT, max(3, 0.22 * nrow(shown) + 2.2))
walk(c("svg", "png"), \(ext) ggsave(str_c("ora_", info$contrast, ".", ext), plot, width = WIDTH, height = height, dpi = DPI))
