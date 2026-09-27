# =============================================================================
# 11b_tree_perserovar.R — per-serovar SNP trees (snippy/Gubbins/IQ-TREE)
# Reads : 09_phylogeny/02_snippy/<Serovar>/iqtree/<Serovar>_snp.treefile
# Writes: 11_figures/fig2_perserovar_trees.png
# Notes : trees include a 'Reference' tip (snippy-core implicit outgroup).
#         Flagged tips: Infantis 746/758 (pESI- sublineage),
#         Enteritidis ERR3843447 (outlier), Kentucky ERR9714962 (outlier).
# =============================================================================
suppressPackageStartupMessages({
  library(ggtree); library(treeio); library(tidyverse); library(patchwork)
})

phy_dir <- Sys.getenv("PHYLO_DIR")
fig_dir <- Sys.getenv("FIGURES_DIR")

serovars <- c("Typhimurium", "Enteritidis", "Infantis", "Kentucky")
flag <- list(Infantis    = c("SRR30149746", "SRR30149758"),
             Enteritidis = c("ERR3843447"),
             Kentucky    = c("ERR9714962"))

make_tree <- function(s) {
  tf <- file.path(phy_dir, "02_snippy", s, "iqtree", paste0(s, "_snp.treefile"))
  if (!file.exists(tf)) { message("missing ", tf); return(NULL) }
  tr <- read.tree(tf)
  hl <- flag[[s]]
  d  <- tibble(label = tr$tip.label,
               kind = case_when(label == "Reference" ~ "reference",
                                label %in% hl        ~ "flagged",
                                TRUE                 ~ "isolate"))
  ggtree(tr) %<+% d +
    geom_tippoint(aes(colour = kind), size = 2, show.legend = FALSE) +
    scale_colour_manual(values = c(reference = "grey50",
                                   flagged   = "firebrick",
                                   isolate   = "black")) +
    geom_tiplab(size = 2.2) +
    theme_tree2() +
    ggtitle(s) +
    theme(plot.title = element_text(size = 11)) +
    hexpand(0.35) +                     # room for right-edge labels
    coord_cartesian(clip = "off")
}

plots <- Filter(Negate(is.null), lapply(serovars, make_tree))

ggsave(file.path(fig_dir, "fig2_perserovar_trees.png"),
       wrap_plots(plots, ncol = 2), width = 14, height = 10, dpi = 300)
cat("wrote fig2_perserovar_trees.png\n")

