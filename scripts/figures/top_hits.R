# Intensities of the most significant proteins, z-scored per protein (red = above the protein's mean).
# Samples are grouped by condition. Edit the settings, then run in this folder: Rscript plot.R
suppressPackageStartupMessages(library(tidyverse))

WIDTH     <- 9    # inches
HEIGHT    <- NA   # inches; NA = grows with the number of proteins
FONT_SIZE <- 7
GRADIENT  <- c("#2166AC", "white", "#B2182B")  # low, mean, high
COLOURS   <- NULL # one per condition, e.g. c("M0 WT" = "#1B9E77", "M0 KO" = "#D95F02"); NULL = default palette
DPI       <- 300  # PNG resolution

samples <- read_tsv("samples.tsv", show_col_types = FALSE)
data    <- read_tsv("data.tsv", show_col_types = FALSE)
heatmap <- as.matrix(data[samples$sample])
rownames(heatmap) <- data$gene
conditions <- fct_inorder(samples$condition)

plot <- pheatmap::pheatmap(
  heatmap,
  scale = "row",
  color = colorRampPalette(GRADIENT)(100),
  cluster_cols = FALSE,
  gaps_col = head(cumsum(rle(samples$condition)$lengths), -1),
  annotation_col = data.frame(condition = conditions, row.names = samples$sample),
  annotation_colors = list(condition = COLOURS %||% set_names(scales::hue_pal()(nlevels(conditions)), levels(conditions))),
  annotation_names_col = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  fontsize = FONT_SIZE,
  silent = TRUE
)

height <- coalesce(HEIGHT, max(4, 0.14 * nrow(data) + 2))
walk(c("svg", "png"), \(ext) ggsave(str_c("top_hits.", ext), plot$gtable, width = WIDTH, height = height, dpi = DPI))
