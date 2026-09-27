#!/bin/bash
# =============================================================================
# Script       : 09_phylo_core.sh
# Description  : Whole-set maximum-likelihood phylogeny from the Panaroo
#                core-gene alignment (all 32 isolates). Shows serovar-level
#                structure across the full panel. IQ-TREE with ModelFinder and
#                1000 ultrafast bootstraps.
#                Single SGE job.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 08_pangenome.sh first (core_gene_alignment.aln).
# Usage        : Run from ${PROJECT_DIR} :  qsub scripts/09_phylo_core.sh
# Dependencies : conda env apha_phylo (iqtree)
# Last updated : 26/09/2026
# =============================================================================
#$ -N sc_09_phylo_core
#$ -cwd
#$ -M m.a.g.rabbani@roslin.ed.ac.uk
#$ -m as
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -P roslin_smith_grp
#$ -o logs/
#$ -e logs/

set -euo pipefail
source config/config.sh
THREADS=4

. /etc/profile.d/modules.sh
module load igmm/apps/anaconda/2023.03
source activate apha_phylo

ALN="${PANGENOME_DIR}/core_gene_alignment.aln"
OUT="${PHYLO_DIR}/01_core_gene_tree"
mkdir -p "${OUT}"

[[ -s "${ALN}" ]] || { echo "ERROR: core alignment not found: ${ALN}"; exit 1; }

echo "Core-gene phylogeny started on $(date)"

# -m MFP     : ModelFinder picks the best substitution model
# -B 1000    : 1000 ultrafast bootstraps
# -T         : threads (AUTO caps at THREADS)
iqtree2 \
    -s "${ALN}" \
    -m MFP \
    -B 1000 \
    -T "${THREADS}" \
    --prefix "${OUT}/core_tree" \
    -redo

echo "Core-gene phylogeny completed on $(date)"
echo "Tree: ${OUT}/core_tree.treefile"

# End of script: 09_phylo_core.sh

