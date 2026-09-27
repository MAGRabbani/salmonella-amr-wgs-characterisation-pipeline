#!/bin/bash
# =============================================================================
# Script       : 07_annotation.sh
# Description  : Genome annotation of isolate assemblies with Bakta.
#                Produces GFF3 (for Panaroo pan-genome), GenBank, TSV feature
#                tables and FAA/FFN per isolate. Genus/species set to Salmonella
#                enterica for consistent gene calling.
#                Runs as an SGE array job (one assembly per task).
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 03_assembly.sh first. Bakta DB downloaded (BAKTA_DB).
# Usage        : Run from ${PROJECT_DIR} :  qsub scripts/07_annotation.sh
# Dependencies : conda env apha_annotation (bakta 1.9.4)
# Last updated : 26/09/2026
# =============================================================================

# =============================================================================
# Grid Engine options
#   -t 1-32      : array job, one task per assembly
#   -tc 10       : throttle - max 10 tasks at once
#   -m as        : mail on abort/suspend
# =============================================================================
#$ -N sc_07_annotation
#$ -cwd
#$ -M m.a.g.rabbani@roslin.ed.ac.uk
#$ -m as
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -P roslin_smith_grp
#$ -o logs/
#$ -e logs/
#$ -t 1-32
#$ -tc 10

set -euo pipefail

# ---- Load config ------------------------------------------------------------
source config/config.sh
THREADS=4

# ---- Load environment -------------------------------------------------------
. /etc/profile.d/modules.sh
module load igmm/apps/anaconda/2023.03
source activate apha_annotation

# ---- Output directory -------------------------------------------------------
mkdir -p "${ANNOT_DIR}"

# ---- Select assembly for this array task ------------------------------------
ASM=$(ls "${ASSEMBLY_DIR}"/*.fasta | sed -n "${SGE_TASK_ID}p")
[[ -n "${ASM}" ]] || { echo "ERROR: No assembly at line ${SGE_TASK_ID}"; exit 1; }
ACC=$(basename "${ASM}" .fasta)

echo "Annotation started for ${ACC} as task=${SGE_TASK_ID} on $(date)"

# ---- Bakta annotation -------------------------------------------------------
# --genus/--species: consistent gene calling for Salmonella enterica
# --locus-tag: per-isolate tag for clean downstream IDs
# --prefix: names all output files by accession
# --skip-plot: skips circular genome plots (not needed; saves time/deps)
OUT="${ANNOT_DIR}/${ACC}"
mkdir -p "${OUT}"
bakta \
    --db "${BAKTA_DB}" \
    --genus Salmonella \
    --species enterica \
    --prefix "${ACC}" \
    --locus-tag "${ACC}" \
    --output "${OUT}" \
    --threads "${THREADS}" \
    --skip-plot \
    --force \
    "${ASM}"

echo "Annotation completed for ${ACC} as task=${SGE_TASK_ID} on $(date)"

# End of script: 07_annotation.sh

