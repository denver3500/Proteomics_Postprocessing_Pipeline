# Volcano plot of one contrast; dashed lines mark the significance cutoffs.
# data.tsv: one row per protein; contrast.tsv: title and cutoffs. Edit the settings, then run in this folder: Rscript plot.R
suppressPackageStartupMessages(library(tidyverse))

WIDTH      <- 6     # inches
HEIGHT     <- 5     # inches
BASE_SIZE  <- 11    # font size
POINT_SIZE <- 0.8
COLOURS    <- c(up = "#B2182B", down = "#2166AC", ns = "grey75")
N_LABELS   <- 10    # most significant proteins labelled
HIGHLIGHT  <- c()   # genes always labelled and circled, e.g. c("TRAP1", "GLS")
DPI        <- 300   # PNG resolution

data <- read_tsv("data.tsv", show_col_types = FALSE)
info <- read_tsv("contrast.tsv", show_col_types = FALSE)

labelled <- data |>
  filter(direction != "ns") |>
  slice_min(pvalue, n = N_LABELS, with_ties = FALSE) |>
  bind_rows(filter(data, gene %in% HIGHLIGHT)) |>
  distinct() |>
  mutate(text_colour = if_else(gene %in% HIGHLIGHT, "black", COLOURS[direction]))

plot <- ggplot(arrange(data, direction != "ns"), aes(log2fc, -log10(pvalue), colour = direction)) +
  geom_point(size = POINT_SIZE, alpha = 0.7) +
  geom_point(data = filter(labelled, gene %in% HIGHLIGHT), shape = 21, size = 2.5, colour = "black") +
  (if (!is.na(info$p_cutoff)) geom_hline(yintercept = -log10(info$p_cutoff), linetype = "dashed", colour = "grey40")) +
  (if (!is.na(info$fc_cutoff)) geom_vline(xintercept = c(-1, 1) * info$fc_cutoff, linetype = "dashed", colour = "grey40")) +
  (if (nrow(labelled)) ggrepel::geom_text_repel(aes(label = gene), data = labelled, colour = labelled$text_colour,
                                                size = 2.5, max.overlaps = Inf, show.legend = FALSE)) +
  annotate("label", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, label = str_glue("{sum(data$direction == 'down')} down"),
           colour = COLOURS[["down"]], fontface = "bold", linewidth = 0) +
  annotate("label", x = Inf, y = Inf, hjust = 1.1, vjust = 1.3, label = str_glue("{sum(data$direction == 'up')} up"),
           colour = COLOURS[["up"]], fontface = "bold", linewidth = 0) +
  scale_colour_manual(values = COLOURS, breaks = c("up", "down"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(title = info$title, subtitle = info$subtitle, x = "log2 fold change", y = "-log10 p-value") +
  theme_grey(BASE_SIZE)

walk(c("svg", "png"), \(ext) ggsave(str_c("volcano_", info$contrast, ".", ext), plot, width = WIDTH, height = HEIGHT, dpi = DPI))
