# p-value histogram per contrast: flat with a peak near 0 is healthy.
# data.tsv: one row per protein, one column per contrast. Edit the settings, then run in this folder: Rscript plot.R
suppressPackageStartupMessages(library(tidyverse))

WIDTH     <- 9         # inches
HEIGHT    <- NA        # inches; NA = grows with the number of contrasts
BASE_SIZE <- 11        # font size
COLUMNS   <- 4         # panels per row
BIN_WIDTH <- 0.05
FILL      <- "grey50"
DPI       <- 300       # PNG resolution

data <- read_tsv("data.tsv", show_col_types = FALSE) |>
  pivot_longer(-protein_id, names_to = "contrast", values_to = "pvalue") |>
  mutate(contrast = fct_inorder(contrast))

plot <- ggplot(data, aes(pvalue)) +
  geom_histogram(breaks = seq(0, 1, BIN_WIDTH), fill = FILL, na.rm = TRUE) +
  facet_wrap(vars(contrast), ncol = COLUMNS, scales = "free_y") +
  labs(x = "p-value", y = "proteins") +
  theme_grey(BASE_SIZE) +
  theme(strip.text = element_text(size = 7))

height <- coalesce(HEIGHT, max(2.5, 1.8 * ceiling(nlevels(data$contrast) / COLUMNS)))
walk(c("svg", "png"), \(ext) ggsave(str_c("pvalues.", ext), plot, width = WIDTH, height = height, dpi = DPI))
