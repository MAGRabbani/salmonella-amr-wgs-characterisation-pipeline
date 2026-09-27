#!/bin/bash
# =============================================================================
# Script       : 02a_multiqc_trim.sh
# Description  : Aggregate raw + trimmed FastQC and fastp reports into a single
#                MultiQC report (before/after trimming). Single (non-array) job.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 01_fastqc.sh, 01a_trim.sh, 01b_fastqc_trim.sh first.
# Usage        : Run from ${PROJECT_DIR} :  qsub scripts/02a_multiqc_trim.sh
# Dependencies : conda env apha_wgs (multiqc)
# Last updated : 24/09/2026
# =============================================================================

# =============================================================================
# Grid Engine options
# =============================================================================
#$ -N sc_02a_multiqc_trim
#$ -cwd
#$ -M m.a.g.rabbani@roslin.ed.ac.uk
#$ -m beas
#$ -l h_vmem=8G
#$ -pe sharedmem 1
#$ -P roslin_smith_grp
#$ -o logs/
#$ -e logs/

set -euo pipefail

# ---- Load config ------------------------------------------------------------
source config/config.sh

# ---- Load environment -------------------------------------------------------
. /etc/profile.d/modules.sh
module load igmm/apps/anaconda/2023.03
source activate "${ENV_WGS}"

# ---- Output -----------------------------------------------------------------
mkdir -p "${MULTIQC_DIR}"

# ---- Checks -----------------------------------------------------------------
[[ -d "${FASTQC_RAW_DIR}"  ]] || { echo "ERROR: Raw FastQC dir missing: ${FASTQC_RAW_DIR}"; exit 1; }
[[ -d "${FASTQC_TRIM_DIR}" ]] || { echo "ERROR: Trimmed FastQC dir missing: ${FASTQC_TRIM_DIR}"; exit 1; }
[[ -d "${FASTP_REPORTS}"   ]] || { echo "ERROR: fastp reports dir missing: ${FASTP_REPORTS}"; exit 1; }

# ---- Run MultiQC ------------------------------------------------------------
echo "MultiQC (before/after trimming) started on $(date)"

multiqc "${FASTQC_RAW_DIR}" "${FASTQC_TRIM_DIR}" "${FASTP_REPORTS}" \
    -o "${MULTIQC_DIR}" \
    -n "02a_multiqc_raw_vs_trimmed" \
    --title "APHA Salmonella - QC before vs after trimming" \
    --fn_as_s_name \
    --fullnames \
    -f

echo "MultiQC (before/after trimming) completed on $(date)"

# End of script: 02a_multiqc_trim.sh
