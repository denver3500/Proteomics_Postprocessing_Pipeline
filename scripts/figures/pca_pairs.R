# A separate PCA for every pair of conditions, one panel each; points are labelled with replicate numbers.
# Edit the settings, then run in this folder: Rscript plot.R
suppressPackageStartupMessages(library(tidyverse))

WIDTH      <- 9     # inches
HEIGHT     <- NA    # inches; NA = grows with the number of pairs
BASE_SIZE  <- 11    # font size
COLUMNS    <- 3     # panels per row
POINT_SIZE <- 1.8
COLOURS    <- NULL  # one per condition, e.g. c("M0 WT" = "#1B9E77", "M0 KO" = "#D95F02"); NULL = default palette
DPI        <- 300   # PNG resolution

data <- read_tsv("data.tsv", show_col_types = FALSE) |>
  mutate(
    condition = fct_inorder(condition),
    panel = fct_inorder(str_c(pair, "\nPC1 ", var1, "%, PC2 ", var2, "%"))
  )

plot <- ggplot(data, aes(PC1, PC2, colour = condition)) +
  geom_point(size = POINT_SIZE) +
  geom_text(aes(label = str_c("R", replicate)), size = 2.5, vjust = -0.8, show.legend = FALSE) +
  facet_wrap(vars(panel), scales = "free", ncol = COLUMNS) +
  (if (!is.null(COLOURS)) scale_colour_manual(values = COLOURS)) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  theme_grey(BASE_SIZE) +
  theme(legend.position = "top", axis.text = element_text(size = 6))

height <- coalesce(HEIGHT, 3 * ceiling(nlevels(data$panel) / COLUMNS))
walk(c("svg", "png"), \(ext) ggsave(str_c("pca_pairs.", ext), plot, width = WIDTH, height = height, dpi = DPI))
