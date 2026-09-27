# =============================================================================
# 11c_amr_matrix.R — AMR gene presence/absence heatmap
# Reads : 05_amr_virulence/03_amr_matrix_annotated.tsv (serovar, ST, isolate, genes)
#         06_plasmid_mge/03_amr_gene_location.tsv (molecule_type per gene)
# Writes: 11_figures/fig3_amr_matrix.png
# Notes : column annotation = dominant location per gene (chromosome/plasmid/both).
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(pheatmap)
})

amr_dir  <- Sys.getenv("AMR_DIR")
plas_dir <- Sys.getenv("PLASMID_DIR")
fig_dir  <- Sys.getenv("FIGURES_DIR")

am <- read_tsv(file.path(amr_dir, "03_amr_matrix_annotated.tsv"),
               show_col_types = FALSE)

# matrix: rows = isolates, cols = genes
gene_cols <- setdiff(names(am), c("serovar", "ST", "isolate"))
m <- as.matrix(am[, gene_cols])
rownames(m) <- am$isolate
storage.mode(m) <- "numeric"

# row annotation: serovar (ordered so clades group together)
am <- am %>% arrange(serovar, isolate)
m  <- m[am$isolate, , drop = FALSE]
row_anno <- data.frame(serovar = am$serovar)
rownames(row_anno) <- am$isolate

# column annotation: dominant location per gene
loc <- read_tsv(file.path(plas_dir, "03_amr_gene_location.tsv"),
                show_col_types = FALSE) %>%
  count(gene, molecule_type) %>%
  group_by(gene) %>%
  summarise(location = if (n_distinct(molecule_type) > 1) "both"
                       else molecule_type[which.max(n)],
            .groups = "drop")

col_anno <- tibble(gene = gene_cols) %>%
  left_join(loc, by = "gene") %>%
  mutate(location = replace_na(location, "unknown")) %>%
  column_to_rownames("gene")

ann_colours <- list(location = c(chromosome = "#1b7837",
                                 plasmid    = "#762a83",
                                 both       = "#e08214",
                                 unknown    = "grey80"))

png(file.path(fig_dir, "fig3_amr_matrix.png"),
    width = 2600, height = 1500, res = 200)
pheatmap(m,
         color = c("grey92", "#2c7fb8"),
         legend_breaks = c(0, 1), legend_labels = c("absent", "present"),
         cluster_rows = FALSE, cluster_cols = TRUE,
         annotation_row = row_anno,
         annotation_col = col_anno,
         annotation_colors = ann_colours,
         fontsize_row = 7, fontsize_col = 8,
         main = "AMR gene carriage by isolate (chromosomal vs plasmid-borne)")
dev.off()
cat("wrote fig3_amr_matrix.png\n")

