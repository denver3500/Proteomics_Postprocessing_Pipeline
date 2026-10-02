# One MSigDB collection (msigdbr) -> gene sets for the dataset's species, restricted to genes detected in the dataset
source("scripts/common.R")

proteins <- read_tsv(snakemake@input$intensities, show_col_types = FALSE)
species  <- detect_species(proteins)
msigdb   <- msigdb_source(snakemake@params$msigdb, species)
code     <- str_split(msigdb$code, ":", n = 2)[[1]]  # C5:GO:BP -> collection C5, subcollection GO:BP

sets <- msigdbr::msigdbr(
  species = species$name,
  db_species = msigdb$db,
  collection = code[1],
  subcollection = if (length(code) == 2) code[2]
)

# id: the source database's identifier (GO:0006096, hsa00010, R-HSA-...) when there is one, else MSigDB's own
sets |>
  transmute(
    term = gs_name,
    id = if_else(str_detect(coalesce(gs_exact_source, ""), "^\\S+$"), gs_exact_source, gs_id),
    gene = gene_symbol
  ) |>
  filter(gene %in% first_gene(proteins$gene)) |>
  distinct() |>
  write_tsv(snakemake@output[[1]])
