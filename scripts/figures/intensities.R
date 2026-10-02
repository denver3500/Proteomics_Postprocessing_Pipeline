# Intensity distribution per sample. Edit the settings, then run in this folder: Rscript plot.R
suppressPackageStartupMessages(library(tidyverse))

WIDTH     <- 9     # inches
HEIGHT    <- NA    # inches; NA = grows with the number of samples
BASE_SIZE <- 11    # font size
COLOURS   <- NULL  # one per condition, e.g. c("M0 WT" = "#1B9E77", "M0 KO" = "#D95F02"); NULL = default palette
DPI       <- 300   # PNG resolution

samples <- read_tsv("samples.tsv", show_col_types = FALSE)
data <- read_tsv("data.tsv", show_col_types = FALSE) |>
  pivot_longer(-protein_id, names_to = "sample", values_to = "log2") |>
  left_join(samples, by = "sample") |>
  mutate(sample = factor(sample, samples$sample), condition = factor(condition, unique(samples$condition)))

plot <- ggplot(data, aes(log2, fct_rev(sample), fill = condition)) +
  geom_boxplot(outlier.size = 0.3, na.rm = TRUE) +
  (if (!is.null(COLOURS)) scale_fill_manual(values = COLOURS)) +
  labs(x = "log2 intensity", y = NULL) +
  theme_grey(BASE_SIZE)

height <- coalesce(HEIGHT, max(4, 0.18 * nrow(samples)))
walk(c("svg", "png"), \(ext) ggsave(str_c("intensities.", ext), plot, width = WIDTH, height = height, dpi = DPI))
