# Coefficient of variation of each protein across the replicates of a condition (lower = more reproducible).
# Edit the settings, then run in this folder: Rscript plot.R
suppressPackageStartupMessages(library(tidyverse))

WIDTH     <- 9     # inches
HEIGHT    <- NA    # inches; NA = grows with the number of conditions
BASE_SIZE <- 11    # font size
MAX_CV    <- 100   # x axis limit, %
COLOURS   <- NULL  # one per condition, e.g. c("M0 WT" = "#1B9E77", "M0 KO" = "#D95F02"); NULL = default palette
DPI       <- 300   # PNG resolution

data <- read_tsv("data.tsv", show_col_types = FALSE) |>
  pivot_longer(-protein_id, names_to = "condition", values_to = "cv") |>
  mutate(condition = fct_inorder(condition))

plot <- ggplot(data, aes(cv, fct_rev(condition), fill = condition)) +
  geom_boxplot(outlier.shape = NA, na.rm = TRUE, show.legend = FALSE) +
  coord_cartesian(xlim = c(0, MAX_CV)) +
  (if (!is.null(COLOURS)) scale_fill_manual(values = COLOURS)) +
  labs(x = "CV % per protein", y = NULL) +
  theme_grey(BASE_SIZE)

height <- coalesce(HEIGHT, max(2.5, 0.45 * nlevels(data$condition)))
walk(c("svg", "png"), \(ext) ggsave(str_c("cv.", ext), plot, width = WIDTH, height = height, dpi = DPI))
