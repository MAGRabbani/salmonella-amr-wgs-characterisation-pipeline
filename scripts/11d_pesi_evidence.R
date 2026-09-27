# =============================================================================
# 11d_pesi_evidence.R — pESI carriage evidence panel (3 long-read Infantis)
# Reads : 07_annotation/01_annotation_summary.tsv (size_Mb, CDS)
#         10_longread_hybrid/01_assembly_overview.tsv (total_plasmid_bp)
# Writes: 11_figures/fig4_pesi_evidence.png
# Notes : 748 = pESI+, 746 & 758 = pESI-negative.
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(patchwork)
})

annot_dir <- Sys.getenv("ANNOT_DIR")
lr_dir    <- Sys.getenv("LONGREAD_DIR")
fig_dir   <- Sys.getenv("FIGURES_DIR")

trio <- c("SRR30149748", "SRR30149746", "SRR30149758")

ann <- read_tsv(file.path(annot_dir, "01_annotation_summary.tsv"),
                show_col_types = FALSE) %>%
  filter(isolate %in% trio) %>%
  select(isolate, size_Mb, CDS)

lr <- read_tsv(file.path(lr_dir, "01_assembly_overview.tsv"),
               show_col_types = FALSE) %>%
  filter(isolate %in% trio) %>%
  mutate(plasmid_Mb = total_plasmid_bp / 1e6) %>%
  select(isolate, plasmid_Mb)

d <- ann %>%
  left_join(lr, by = "isolate") %>%
  mutate(pESI = if_else(isolate == "SRR30149748", "pESI+", "pESI-"),
         isolate = factor(isolate, levels = trio)) %>%
  pivot_longer(c(size_Mb, CDS, plasmid_Mb),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         size_Mb    = "Genome size (Mb)",
                         CDS        = "CDS count",
                         plasmid_Mb = "Plasmid content (Mb)"))

p <- ggplot(d, aes(isolate, value, fill = pESI)) +
  geom_col(width = 0.65) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = c("pESI+" = "#762a83", "pESI-" = "grey60")) +
  labs(x = NULL, y = NULL, fill = NULL,
       title = "pESI carriage: SRR30149748 vs pESI-negative pair") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(size = 12))

ggsave(file.path(fig_dir, "fig4_pesi_evidence.png"),
       p, width = 10, height = 4.2, dpi = 300)
cat("wrote fig4_pesi_evidence.png\n")

