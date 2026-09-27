#!/bin/bash
# =============================================================================
# Script       : 02_multiqc.sh
# Description  : Aggregate FastQC reports into a single MultiQC report.
#                Single (non-array) SGE job.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Usage        : Run from ${PROJECT_DIR} :  qsub scripts/02_multiqc.sh
# Dependencies : conda env apha_wgs (multiqc)
# =============================================================================

# =============================================================================
# Grid Engine options
# =============================================================================
#$ -N sc_02_multiqc
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

# ---- Input / output ---------------------------------------------------------
FASTQC_DIR="${QC_DIR}/01_fastqc"
mkdir -p "${MULTIQC_DIR}"

# ---- Checks -----------------------------------------------------------------
[[ -d "${FASTQC_DIR}" ]] || { echo "ERROR: FastQC dir missing: ${FASTQC_DIR}"; exit 1; }

# ---- Run MultiQC ------------------------------------------------------------
echo "MultiQC started on $(date)"

multiqc "${FASTQC_DIR}" \
    -o "${MULTIQC_DIR}" \
    -n "01_multiqc_of_raw_fastqc" \
    --title "APHA Salmonella - Raw Read QC" \
    -f

echo "MultiQC completed on $(date)"

# End of script: 02_multiqc.sh

