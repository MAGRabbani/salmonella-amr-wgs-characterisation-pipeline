#!/bin/bash
# =============================================================================
# Script       : 05_amr_virulence.sh
# Description  : AMR and virulence gene detection on isolate assemblies.
#                AMRFinderPlus (organism-aware: acquired genes + point mutations)
#                + ABRicate (CARD, ResFinder, VFDB).
#                Runs as an SGE array job (one assembly per task).
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 03_assembly.sh and setup_databases.sh first.
# Usage        : Run from ${PROJECT_DIR} :  qsub scripts/05_amr_virulence.sh
# Dependencies : conda env apha_amrPlasfinder (ncbi-amrfinderplus, abricate)
# Last updated : 26/09/2026
# =============================================================================

# =============================================================================
# Grid Engine options
#   -t 1-32      : array job, one task per assembly
#   -tc 10       : throttle - max 10 tasks at once
#   -m as        : mail on abort/suspend
# =============================================================================
#$ -N sc_05_amr_virulence
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
source activate apha_amrPlasfinder

# ---- Output directories -----------------------------------------------------
AMRFINDER_DIR="${AMR_DIR}/01_amrfinderplus"
ABRICATE_DIR="${AMR_DIR}/02_abricate"
mkdir -p "${AMRFINDER_DIR}" "${ABRICATE_DIR}"

# ---- Select assembly for this array task ------------------------------------
# One isolate per task, indexed from the assembly file list.
ASM=$(ls "${ASSEMBLY_DIR}"/*.fasta | sed -n "${SGE_TASK_ID}p")
[[ -n "${ASM}" ]] || { echo "ERROR: No assembly at line ${SGE_TASK_ID}"; exit 1; }
ACC=$(basename "${ASM}" .fasta)

echo "AMR/virulence started for ${ACC} as task=${SGE_TASK_ID} on $(date)"

# ---- 1) AMRFinderPlus (organism-aware: acquired genes + point mutations) -----
# --organism Salmonella enables Salmonella point-mutation detection
# (e.g. gyrA/parC for fluoroquinolone resistance).
amrfinder \
    --nucleotide "${ASM}" \
    --organism Salmonella \
    --threads "${THREADS}" \
    --name "${ACC}" \
    --output "${AMRFINDER_DIR}/${ACC}.amrfinder.tsv"

# ---- 2) ABRicate (acquired genes + virulence, multiple databases) ------------
for DB in card resfinder vfdb; do
    abricate \
        --db "${DB}" \
        --threads "${THREADS}" \
        "${ASM}" > "${ABRICATE_DIR}/${ACC}.${DB}.tab"
done

echo "AMR/virulence completed for ${ACC} as task=${SGE_TASK_ID} on $(date)"

# End of script: 05_amr_virulence.sh

