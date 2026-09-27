#!/bin/bash
# =============================================================================
# Script       : 01b_fastqc_trim.sh
# Description  : Quality control of trimmed short-read isolates (FastQC).
#                Runs as an SGE array job (one accession per task).
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 01a_trim.sh first.
# Usage        : Run from ${PROJECT_DIR} :  qsub scripts/01b_fastqc_trim.sh
# Dependencies : conda env apha_wgs (fastqc)
# Last updated : 24/09/2026
# =============================================================================

# =============================================================================
# Grid Engine options
#   -t 1-33      : array job, one task per short-read accession
#   -tc 10       : throttle - max 10 tasks at once
#   -m as        : mail on abort/suspend
# =============================================================================
#$ -N sc_01b_fastqc_trim
#$ -cwd
#$ -M m.a.g.rabbani@roslin.ed.ac.uk
#$ -m as
#$ -l h_vmem=8G
#$ -pe sharedmem 8
#$ -P roslin_smith_grp
#$ -o logs/
#$ -e logs/
#$ -t 1-33
#$ -tc 10

set -euo pipefail

# ---- Load config ------------------------------------------------------------
source config/config.sh

# ---- Load environment -------------------------------------------------------
. /etc/profile.d/modules.sh
module load igmm/apps/anaconda/2023.03
source activate "${ENV_WGS}"

# ---- Output directory -------------------------------------------------------
mkdir -p "${FASTQC_TRIM_DIR}"

# ---- Select accession for this array task -----------------------------------
# ${SGE_TASK_ID} is the task number Grid Engine gives each array task (1..33).
# sed prints that line from short_read_ids.txt, so each task handles one isolate.
IDS="${SAMPLES_DIR}/short_read_ids.txt"
ACC=$(sed -n "${SGE_TASK_ID}p" "${IDS}")
[[ -n "${ACC}" ]] || { echo "ERROR: No accession at line ${SGE_TASK_ID}"; exit 1; }

# ---- Checks -----------------------------------------------------------------
ls "${TRIM_DIR}/${ACC}"_*.trimmed.fastq.gz >/dev/null 2>&1 \
    || { echo "ERROR: No trimmed FASTQ files for ${ACC} in ${TRIM_DIR}"; exit 1; }

# ---- Run FastQC -------------------------------------------------------------
echo "FastQC (trimmed) started for ${ACC} as task=${SGE_TASK_ID} on $(date)"

fastqc "${TRIM_DIR}/${ACC}"_*.trimmed.fastq.gz -o "${FASTQC_TRIM_DIR}" -t "${THREADS}"

echo "FastQC (trimmed) completed for ${ACC} as task=${SGE_TASK_ID} on $(date)"

# End of script: 01b_fastqc_trim.sh
