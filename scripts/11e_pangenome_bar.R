# =============================================================================
# 11e_pangenome_bar.R — pan-genome structure + serovar-specific gene counts
# Reads : 08_pangenome/gene_presence_absence.Rtab (categories computed here)
#         08_pangenome/02_serovar_specific_genes.tsv
# Writes: 11_figures/fig5_pangenome.png
# Notes : category thresholds (Panaroo): core >=99%, soft 95-99%,
#         shell 15-95%, cloud <15% of isolates.
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(patchwork)
})

pan_dir <- Sys.getenv("PANGENOME_DIR")
fig_dir <- Sys.getenv("FIGURES_DIR")

rtab <- read_tsv(file.path(pan_dir, "gene_presence_absence.Rtab"),
                 show_col_types = FALSE)
iso_cols <- setdiff(names(rtab), names(rtab)[1])   # first col = Gene
n_iso <- length(iso_cols)

freq <- rowSums(rtab[, iso_cols]) / n_iso
cat_counts <- tibble(freq) %>%
  mutate(category = case_when(freq >= 0.99 ~ "Core",
                              freq >= 0.95 ~ "Soft-core",
                              freq >= 0.15 ~ "Shell",
                              TRUE         ~ "Cloud")) %>%
  count(category) %>%
  mutate(category = factor(category,
                           levels = c("Core","Soft-core","Shell","Cloud")))

p1 <- ggplot(cat_counts, aes(category, n, fill = category)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.2) +
  scale_fill_viridis_d(end = 0.9) +
  labs(x = NULL, y = "Gene clusters", title = "Pan-genome structure") +
  theme_bw() + theme(plot.title = element_text(size = 12))

spec <- read_tsv(file.path(pan_dir, "02_serovar_specific_genes.tsv"),
                 show_col_types = FALSE) %>%
  arrange(desc(n_specific_genes)) %>%
  mutate(serovar = factor(serovar, levels = serovar))

p2 <- ggplot(spec, aes(serovar, n_specific_genes, fill = serovar)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = n_specific_genes), vjust = -0.3, size = 3.2) +
  scale_fill_viridis_d(end = 0.9) +
  labs(x = NULL, y = "Specific genes", title = "Serovar-specific genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(size = 12))

ggsave(file.path(fig_dir, "fig5_pangenome.png"),
       p1 + p2, width = 11, height = 4.5, dpi = 300)
cat("wrote fig5_pangenome.png\n")

