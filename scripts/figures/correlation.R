# Sample-to-sample Pearson correlation, clustered on 1 - r; the colour bar marks each sample's condition.
# Edit the settings, then run in this folder: Rscript plot.R
suppressPackageStartupMessages(library(tidyverse))

WIDTH     <- 9          # inches
HEIGHT    <- NA         # inches; NA = grows with the number of samples
FONT_SIZE <- 7
PALETTE   <- "viridis"  # any name from hcl.pals()
COLOURS   <- NULL       # one per condition, e.g. c("M0 WT" = "#1B9E77", "M0 KO" = "#D95F02"); NULL = default palette
DPI       <- 300        # PNG resolution

data <- read_tsv("data.tsv", show_col_types = FALSE)
cors <- as.matrix(data[data$sample])
rownames(cors) <- data$sample
conditions <- fct_inorder(data$condition)
annotation <- data.frame(condition = conditions, row.names = data$sample)
distance <- as.dist(1 - cors)

plot <- pheatmap::pheatmap(
  cors,
  color = hcl.colors(100, PALETTE),
  clustering_distance_rows = distance,
  clustering_distance_cols = distance,
  clustering_method = "average",
  annotation_row = annotation,
  annotation_col = annotation,
  annotation_colors = list(condition = COLOURS %||% set_names(scales::hue_pal()(nlevels(conditions)), levels(conditions))),
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  border_color = NA,
  fontsize = FONT_SIZE,
  silent = TRUE
)

height <- coalesce(HEIGHT, max(4, 0.18 * nrow(data)) + 1.5)
walk(c("svg", "png"), \(ext) ggsave(str_c("correlation.", ext), plot$gtable, width = WIDTH, height = height, dpi = DPI))
