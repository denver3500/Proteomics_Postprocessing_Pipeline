# PCA of all samples on the proteins measured in every sample; flagged samples are labelled with their replicate.
# Edit the settings, then run in this folder: Rscript plot.R
suppressPackageStartupMessages(library(tidyverse))

WIDTH      <- 9     # inches
HEIGHT     <- 6     # inches
BASE_SIZE  <- 11    # font size
POINT_SIZE <- 2.5
COLOURS    <- NULL  # one per condition, e.g. c("M0 WT" = "#1B9E77", "M0 KO" = "#D95F02"); NULL = default palette
DPI        <- 300   # PNG resolution

data <- read_tsv("data.tsv", show_col_types = FALSE) |>
  mutate(condition = fct_inorder(condition))

plot <- ggplot(data, aes(PC1, PC2, colour = condition)) +
  geom_point(size = POINT_SIZE) +
  geom_text(aes(label = str_c("R", replicate)), data = filter(data, !is.na(flags)), size = 3, vjust = -0.9, show.legend = FALSE) +
  (if (!is.null(COLOURS)) scale_colour_manual(values = COLOURS)) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  labs(x = str_glue("PC1 ({data$var1[1]}%)"), y = str_glue("PC2 ({data$var2[1]}%)")) +
  theme_grey(BASE_SIZE)

walk(c("svg", "png"), \(ext) ggsave(str_c("pca.", ext), plot, width = WIDTH, height = HEIGHT, dpi = DPI))
