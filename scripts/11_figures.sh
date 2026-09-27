#!/bin/bash
# =============================================================================
# Script       : 11_figures.sh
# Description  : Render presentation figures from steps 04-10 outputs.
#                Runs five R scripts: whole-set core-gene tree, per-serovar
#                SNP trees, AMR gene matrix (chromosome vs plasmid), pESI
#                five-way evidence panel, and pan-genome structure bars.
#                Login-node job (no scheduler): R is single-threaded here.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Steps 04, 05a, 06a, 07a, 08a, 09, 10a complete; env apha_figures.
# Usage        : Run from ${PROJECT_DIR} :  bash scripts/11_figures.sh
# Dependencies : conda env ENV_FIG (ggtree, treeio, tidyverse, pheatmap, patchwork)
# =============================================================================

set -euo pipefail

# ---- Load config ------------------------------------------------------------
source config/config.sh

# ---- Load environment -------------------------------------------------------
. /etc/profile.d/modules.sh
module load igmm/apps/anaconda/2023.03
source activate "${ENV_FIG}"

# ---- Output directory -------------------------------------------------------
mkdir -p "${FIGURES_DIR}"

# config already exports the vars the R scripts read via Sys.getenv()

for R in 11a_tree_coregene.R \
         11b_tree_perserovar.R \
         11c_amr_matrix.R \
         11d_pesi_evidence.R \
         11e_pangenome_bar.R ; do
    echo "---- running ${R} on $(date) ----"
    Rscript "scripts/${R}" || echo "WARN: ${R} failed"
done

echo "step 11 finished on $(date); figures in ${FIGURES_DIR}"

# End of script: 11_figures.sh
