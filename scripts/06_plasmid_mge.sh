#!/bin/bash
# =============================================================================
# Script       : 06_plasmid_mge.sh
# Description  : Plasmid reconstruction and replicon typing on isolate assemblies.
#                MOB-suite (mob_recon: reconstructs plasmids, assigns replicon/
#                relaxase/mob type, predicts mobility, places genes) + standalone
#                PlasmidFinder (independent replicon call for cross-check).
#                Runs as an SGE array job (one assembly per task).
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 03_assembly.sh first. MOB-suite DB is built-in;
#                PlasmidFinder DB path set in config (PLASMIDFINDER_DB).
# Usage        : Run from ${PROJECT_DIR} :  qsub scripts/06_plasmid_mge.sh
# Dependencies : conda env apha_amrPlasfinder (mob_suite, plasmidfinder)
# Last updated : 26/09/2026
# =============================================================================

# =============================================================================
# Grid Engine options
#   -t 1-32      : array job, one task per assembly
#   -tc 10       : throttle - max 10 tasks at once
#   -m as        : mail on abort/suspend
# =============================================================================
#$ -N sc_06_plasmid_mge
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
MOB_DIR="${PLASMID_DIR}/01_mob_recon"
PF_DIR="${PLASMID_DIR}/02_plasmidfinder"
mkdir -p "${MOB_DIR}" "${PF_DIR}"

# ---- Select assembly for this array task ------------------------------------
# One isolate per task, indexed from the assembly file list.
ASM=$(ls "${ASSEMBLY_DIR}"/*.fasta | sed -n "${SGE_TASK_ID}p")
[[ -n "${ASM}" ]] || { echo "ERROR: No assembly at line ${SGE_TASK_ID}"; exit 1; }
ACC=$(basename "${ASM}" .fasta)

echo "Plasmid/MGE started for ${ACC} as task=${SGE_TASK_ID} on $(date)"

# ---- 1) MOB-suite mob_recon (reconstruction + typing) -----------------------
# Plain defaults, full built-in NCBI plasmid DB, no taxonomy filtering.
# --sample_id used as a label only; --force overwrites a clean re-run.
MOB_OUT="${MOB_DIR}/${ACC}"
rm -rf "${MOB_OUT}"
mob_recon \
    --infile "${ASM}" \
    --outdir "${MOB_OUT}" \
    --sample_id "${ACC}" \
    --num_threads "${THREADS}" \
    --force

# ---- 2) PlasmidFinder (independent replicon call) ---------------------------
# -x extended output (results_tab.tsv). -p points at the DB dir holding config.
PF_OUT="${PF_DIR}/${ACC}"
mkdir -p "${PF_OUT}"
plasmidfinder.py \
    -i "${ASM}" \
    -o "${PF_OUT}" \
    -p "${PLASMIDFINDER_DB}" \
    -x

echo "Plasmid/MGE completed for ${ACC} as task=${SGE_TASK_ID} on $(date)"

# End of script: 06_plasmid_mge.sh

