#!/bin/bash
# =============================================================================
# Script       : 04_serotype_mlst.sh
# Description  : MLST + serotyping (SISTR, SeqSero2) of isolate assemblies.
#                Order: MLST (lineage) -> SISTR (primary serovar) ->
#                SeqSero2 (confirmatory serovar).
#                Runs as an SGE array job (one assembly per task).
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 03_assembly.sh first.
# Usage        : Run from ${PROJECT_DIR} :  qsub scripts/04_serotype_mlst.sh
# Dependencies : conda env apha_typing (mlst, sistr_cmd, seqsero2)
# Last updated : 25/09/2026
# =============================================================================

# =============================================================================
# Grid Engine options
#   -t 1-32      : array job, one task per assembly
#   -tc 10       : throttle - max 10 tasks at once
#   -m as        : mail on abort/suspend
# =============================================================================
#$ -N sc_04_serotype_mlst
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
source activate "${ENV_TYPING}"

# ---- Output directories -----------------------------------------------------
MLST_DIR="${SEROTYPE_DIR}/01_mlst"
SISTR_DIR="${SEROTYPE_DIR}/02_sistr"
SEQSERO_DIR="${SEROTYPE_DIR}/03_seqsero2"
mkdir -p "${MLST_DIR}" "${SISTR_DIR}" "${SEQSERO_DIR}"

# ---- Select assembly for this array task ------------------------------------
# One isolate per task, indexed from the assembly file list.
ASM=$(ls "${ASSEMBLY_DIR}"/*.fasta | sed -n "${SGE_TASK_ID}p")
[[ -n "${ASM}" ]] || { echo "ERROR: No assembly at line ${SGE_TASK_ID}"; exit 1; }
ACC=$(basename "${ASM}" .fasta)

echo "Typing started for ${ACC} as task=${SGE_TASK_ID} on $(date)"

# ---- 1) MLST (7-gene sequence type; lineage) --------------------------------
mlst "${ASM}" > "${MLST_DIR}/${ACC}.mlst.tab"

# ---- 2) SISTR (primary serovar; cgMLST + antigen) ---------------------------
sistr \
    --qc \
    -f tab \
    -o "${SISTR_DIR}/${ACC}.sistr.tab" \
    -t "${THREADS}" \
    "${ASM}"

# ---- 3) SeqSero2 (confirmatory serovar; antigen / Kauffmann-White) ----------
SeqSero2_package.py \
    -m k \
    -t 4 \
    -i "${ASM}" \
    -d "${SEQSERO_DIR}/${ACC}" \
    -p "${THREADS}"

echo "Typing completed for ${ACC} as task=${SGE_TASK_ID} on $(date)"

# End of script: 04_serotype_mlst.sh
