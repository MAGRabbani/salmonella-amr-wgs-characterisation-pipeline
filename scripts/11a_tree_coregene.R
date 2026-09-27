# =============================================================================
# 11a_tree_coregene.R — whole-set core-gene tree, tips coloured by serovar
# Reads : 09_phylogeny/01_core_gene_tree/core_tree.treefile
#         04_serotype_mlst/typing_summary.tsv (col 5 SISTR_serovar)
# Writes: 11_figures/fig1_coregene_tree.png
# =============================================================================
suppressPackageStartupMessages({
  library(ggtree); library(treeio); library(tidyverse)
})

phy_dir  <- Sys.getenv("PHYLO_DIR")
sero_dir <- Sys.getenv("SEROTYPE_DIR")
fig_dir  <- Sys.getenv("FIGURES_DIR")

tree <- read.tree(file.path(phy_dir, "01_core_gene_tree", "core_tree.treefile"))

meta <- read_tsv(file.path(sero_dir, "typing_summary.tsv"),
                 show_col_types = FALSE) %>%
  select(isolate, serovar = SISTR_serovar)

p <- ggtree(tree, layout = "rectangular") %<+% meta +
  geom_tippoint(aes(colour = serovar), size = 2.6) +
  geom_tiplab(size = 2.4, align = TRUE, linesize = 0.2) +
  scale_colour_viridis_d(end = 0.9, name = "Serovar") +
  theme_tree2() +
  ggtitle("Salmonella core-gene phylogeny (32 isolates)") +
  theme(plot.title = element_text(size = 12)) +
  hexpand(0.30) +                       # room on the right for full labels
  coord_cartesian(clip = "off")

ggsave(file.path(fig_dir, "fig1_coregene_tree.png"),
       p, width = 10, height = 8, dpi = 300)
cat("wrote fig1_coregene_tree.png\n")

