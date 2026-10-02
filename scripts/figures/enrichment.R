# Enriched gene sets of one contrast (GSEA or ORA), most significant on top, in up, down and mixed panels.
# Mixed sets are enriched only when up and down proteins are combined. Each bar is split by the share of the set's
# proteins higher in either condition. data.tsv: one row per enriched set; contrast.tsv: titles and labels.
# Edit the settings, then run in this folder: Rscript plot.R
suppressPackageStartupMessages(library(tidyverse))

WIDTH     <- 9     # inches
HEIGHT    <- NA    # inches; NA = grows with the number of sets
BASE_SIZE <- 11    # font size
COLOURS   <- c(up = "#B2182B", down = "#2166AC")
N_SETS    <- 10    # most significant sets shown per panel
MAX_LABEL <- 60    # longer set names are shortened
DPI       <- 300   # PNG resolution

data <- read_tsv("data.tsv", show_col_types = FALSE, col_types = cols(panel = "c", label = "c", note = "c", .default = "d"))
info <- read_tsv("contrast.tsv", show_col_types = FALSE)

shown <- data |>
  mutate(panel = factor(panel, c("up", "down", "mixed"))) |>
  slice_min(padj, n = N_SETS, by = panel, with_ties = FALSE) |>
  mutate(score = -log10(padj), row = fct_reorder(str_c(str_trunc(label, MAX_LABEL), "@", panel), score))

segments <- shown |>
  pivot_longer(c(n_up, n_down), names_to = "direction", names_prefix = "n_", values_to = "n") |>
  mutate(width = score * n / sum(n), .by = row)

plot <- ggplot(shown, aes(y = row)) +
  geom_col(aes(x = width, fill = direction), data = segments, width = 0.75) +
  geom_text(aes(x = score, label = note), hjust = -0.08, size = 0.8 * BASE_SIZE / .pt, colour = "grey25") +
  facet_grid(rows = vars(panel), scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = \(row) str_remove(row, "@[^@]*$")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  scale_fill_manual(values = COLOURS, breaks = c("up", "down"), name = "proteins",
                    labels = c(up = str_c("higher in ", info$higher_up), down = str_c("higher in ", info$higher_down))) +
  labs(title = info$title, subtitle = str_wrap(info$subtitle, 90), x = "-log10 padj", y = NULL,
       caption = if ("mixed" %in% shown$panel) "mixed: enriched only when up and down proteins are combined") +
  theme_grey(BASE_SIZE) +
  theme(legend.position = "bottom", plot.title.position = "plot", strip.text.y = element_text(angle = 0))

if (!nrow(shown)) {
  plot <- ggplot() +
    annotate("text", x = 0, y = 0, label = str_wrap(info$message, 60), size = BASE_SIZE / .pt, colour = "grey25") +
    labs(title = info$title, subtitle = str_wrap(info$subtitle, 90)) +
    theme_void(BASE_SIZE) +
    theme(plot.title.position = "plot", plot.margin = margin(10, 10, 10, 10))
}

height <- coalesce(HEIGHT, if (nrow(shown)) 2 + 0.25 * nrow(shown) + 0.2 * n_distinct(shown$panel) else 2.5)
walk(c("svg", "png"), \(ext) ggsave(str_c(info$method, "_", info$contrast, ".", ext), plot, width = WIDTH, height = height, dpi = DPI))
