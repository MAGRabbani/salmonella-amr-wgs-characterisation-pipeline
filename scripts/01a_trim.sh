#!/bin/bash
# =============================================================================
# Script       : 01a_trim.sh
# Description  : Adapter and quality trimming of raw short-read isolates (fastp).
#                Removes adapters (PE auto-detect), trims polyG tails,
#                light quality trim, min-length filter.
#                Runs as an SGE array job (one accession per task).
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 02_multiqc.sh first (Based on the decision from multiqc report).
# Usage        : Run from ${PROJECT_DIR} :  qsub scripts/01a_trim.sh
# Dependencies : conda env apha_wgs (fastp)
# Last updated : 24/09/2026
# =============================================================================

# =============================================================================
# Grid Engine options
#   -t 1-33      : array job, one task per short-read accession
#   -tc 10       : throttle - max 10 tasks at once
#   -m as        : mail on abort/suspend
# =============================================================================
#$ -N sc_01a_trim
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

# ---- Output directories -----------------------------------------------------
mkdir -p "${TRIM_DIR}" "${FASTP_REPORTS}"

# ---- Select accession for this array task -----------------------------------
# ${SGE_TASK_ID} is the task number Grid Engine gives each array task (1..33).
# sed prints that line from short_read_ids.txt, so each task handles one isolate.
IDS="${SAMPLES_DIR}/short_read_ids.txt"
ACC=$(sed -n "${SGE_TASK_ID}p" "${IDS}")
[[ -n "${ACC}" ]] || { echo "ERROR: No accession at line ${SGE_TASK_ID}"; exit 1; }

R1="${RAW_DIR}/${ACC}_1.fastq.gz"
R2="${RAW_DIR}/${ACC}_2.fastq.gz"

# ---- Checks -----------------------------------------------------------------
[[ -f "${R1}" && -f "${R2}" ]] || { echo "ERROR: Paired files missing for ${ACC}"; exit 1; }

# ---- Run fastp --------------------------------------------------------------
# --detect_adapter_for_pe   : auto-detect and remove adapters (paired-end)
# --trim_poly_g             : remove polyG tails
# --qualified_quality_phred : base quality threshold for the quality filter
# --length_required 50      : discard reads shorter than 50 bp after trimming
# --thread                  : computing threads
# --json / --html           : per-sample reports (JSON feeds MultiQC)
echo "Trimming started for ${ACC} as task=${SGE_TASK_ID} on $(date)"

fastp \
    -i "${R1}" -I "${R2}" \
    -o "${TRIM_DIR}/${ACC}_1.trimmed.fastq.gz" \
    -O "${TRIM_DIR}/${ACC}_2.trimmed.fastq.gz" \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --thread "${THREADS}" \
    --json "${FASTP_REPORTS}/${ACC}.fastp.json" \
    --html "${FASTP_REPORTS}/${ACC}.fastp.html"

echo "Trimming completed for ${ACC} as task=${SGE_TASK_ID} on $(date)"

# End of script: 01a_trim.sh
