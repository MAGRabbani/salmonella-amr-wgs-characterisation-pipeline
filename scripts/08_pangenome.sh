#!/bin/bash
# =============================================================================
# Script       : 08_pangenome.sh
# Description  : Pan-genome analysis of all isolates with Panaroo.
#                Consumes Bakta GFF3 (with sequence) from step 07, builds the
#                gene presence/absence matrix, pan-genome graph, summary stats,
#                and a core-gene alignment (MAFFT) for downstream phylogeny.
#                Single job (all isolates together), not an array.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 07_annotation.sh first (needs *.gff3 with ##FASTA block).
# Usage        : Run from ${PROJECT_DIR} :  qsub scripts/08_pangenome.sh
# Dependencies : conda env apha_pangenome (panaroo, mafft)
# Last updated : 26/09/2026
# =============================================================================

# =============================================================================
# Grid Engine options
#   single job (not an array); Panaroo runs on all 32 genomes together
# =============================================================================
#$ -N sc_08_pangenome
#$ -cwd
#$ -M m.a.g.rabbani@roslin.ed.ac.uk
#$ -m as
#$ -l h_vmem=8G
#$ -pe sharedmem 8
#$ -P roslin_smith_grp
#$ -o logs/
#$ -e logs/

set -euo pipefail

# ---- Load config ------------------------------------------------------------
source config/config.sh
THREADS=8

# ---- Load environment -------------------------------------------------------
. /etc/profile.d/modules.sh
module load igmm/apps/anaconda/2023.03
source activate apha_pangenome

# ---- Output directory -------------------------------------------------------
mkdir -p "${PANGENOME_DIR}"

# ---- Gather GFF3 inputs -----------------------------------------------------
# Bakta writes one GFF3 per isolate (with ##FASTA block) under ANNOT_DIR/<acc>/
GFFS=( "${ANNOT_DIR}"/*/*.gff3 )
N_GFF=${#GFFS[@]}
[[ "${N_GFF}" -gt 0 ]] || { echo "ERROR: no GFF3 files found under ${ANNOT_DIR}"; exit 1; }
echo "Pan-genome started on $(date) with ${N_GFF} genomes"

# ---- Panaroo ----------------------------------------------------------------
# --clean-mode moderate : balanced error correction (recommended default)
# -a core               : build core-gene alignment (feeds step 09 phylogeny)
# --aligner mafft       : per-gene alignment tool
# --core_threshold 0.95 : gene in >=95% isolates counted as core (default)

panaroo \
    -i "${GFFS[@]}" \
    -o "${PANGENOME_DIR}" \
    --clean-mode moderate \
    --remove-invalid-genes \
    -a core \
    --aligner mafft \
    --core_threshold 0.95 \
    -t "${THREADS}"

echo "Pan-genome completed on $(date)"

# End of script: 08_pangenome.sh

