# Share of each sample's intensity from cRAP contaminants; dashed line = dataset median.
# Edit the settings, then run in this folder: Rscript plot.R
suppressPackageStartupMessages(library(tidyverse))

WIDTH     <- 9     # inches
HEIGHT    <- NA    # inches; NA = grows with the number of samples
BASE_SIZE <- 11    # font size
COLOURS   <- NULL  # one per condition, e.g. c("M0 WT" = "#1B9E77", "M0 KO" = "#D95F02"); NULL = default palette
DPI       <- 300   # PNG resolution

data <- read_tsv("data.tsv", show_col_types = FALSE) |>
  mutate(sample = fct_inorder(sample), condition = fct_inorder(condition))

plot <- ggplot(data, aes(contaminants_pct, fct_rev(sample), fill = condition)) +
  geom_col() +
  geom_vline(xintercept = median(data$contaminants_pct), linetype = "dashed") +
  (if (!is.null(COLOURS)) scale_fill_manual(values = COLOURS)) +
  labs(x = "% of total intensity from cRAP contaminants", y = NULL) +
  theme_grey(BASE_SIZE)

height <- coalesce(HEIGHT, max(4, 0.18 * nrow(data)))
walk(c("svg", "png"), \(ext) ggsave(str_c("contaminants.", ext), plot, width = WIDTH, height = height, dpi = DPI))
