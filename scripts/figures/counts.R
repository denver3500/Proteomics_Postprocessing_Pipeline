# Significant proteins per contrast, down to the left and up to the right.
# Edit the settings, then run in this folder: Rscript plot.R
suppressPackageStartupMessages(library(tidyverse))

WIDTH     <- 9    # inches
HEIGHT    <- NA   # inches; NA = grows with the number of contrasts
BASE_SIZE <- 11   # font size
COLOURS   <- c(up = "#B2182B", down = "#2166AC")
DPI       <- 300  # PNG resolution

data <- read_tsv("data.tsv", show_col_types = FALSE)

plot <- data |>
  pivot_longer(c(up, down), names_to = "direction", values_to = "n") |>
  mutate(signed = if_else(direction == "down", -n, n), label = fct_rev(fct_inorder(label))) |>
  ggplot(aes(signed, label, fill = direction)) +
  geom_col() +
  geom_text(aes(label = n, hjust = if_else(direction == "down", 1.15, -0.15)), size = 3) +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = COLOURS) +
  scale_x_continuous(labels = abs, expand = expansion(mult = 0.12)) +
  labs(x = "significant proteins (down | up)", y = NULL, fill = NULL) +
  theme_grey(BASE_SIZE)

height <- coalesce(HEIGHT, max(2.5, 0.3 * nrow(data)))
walk(c("svg", "png"), \(ext) ggsave(str_c("counts.", ext), plot, width = WIDTH, height = height, dpi = DPI))
