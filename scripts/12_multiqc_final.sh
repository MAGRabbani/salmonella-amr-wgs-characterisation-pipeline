#!/bin/bash
# =============================================================================
# Script       : 12_multiqc_final.sh
# Description  : Aggregate all QC across the pipeline into one MultiQC report.
#                Scans FastQC (raw + trimmed), fastp, QUAST and BUSCO outputs
#                from steps 01-03a. CheckM2 (not a MultiQC module) is copied
#                in as its raw quality_report.tsv for reference.
#                Login-node job (no scheduler): MultiQC is lightweight.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Steps 01, 01a, 01b, 03a complete.
# Usage        : Run from ${PROJECT_DIR} :  bash scripts/12_multiqc_final.sh
# Dependencies : conda env apha_wgs (multiqc)
# =============================================================================

set -euo pipefail

# ---- Load config ------------------------------------------------------------
source config/config.sh

# ---- Load environment -------------------------------------------------------
. /etc/profile.d/modules.sh
module load igmm/apps/anaconda/2023.03
source activate "${ENV_WGS}"

# ---- Output directory -------------------------------------------------------
mkdir -p "${MULTIQC_FINAL_DIR}"

# ---- Collect input directories ---------------------------------------------
# Only scan dirs that exist, so a partial run does not abort MultiQC.
SCAN_DIRS=()
for d in "${FASTQC_RAW_DIR}" \
         "${FASTP_REPORTS}" \
         "${FASTQC_TRIM_DIR}" \
         "${QUAST_DIR}" \
         "${BUSCO_DIR}" ; do
    if [[ -d "${d}" ]]; then
        SCAN_DIRS+=("${d}")
        echo "  will scan: ${d}"
    else
        echo "  WARN: missing, skipped: ${d}"
    fi
done

[[ ${#SCAN_DIRS[@]} -gt 0 ]] || { echo "ERROR: no QC dirs found"; exit 1; }

# ---- Run MultiQC ------------------------------------------------------------
echo "MultiQC aggregation started on $(date)"
multiqc \
    "${SCAN_DIRS[@]}" \
    --outdir "${MULTIQC_FINAL_DIR}" \
    --filename multiqc_final_report \
    --title "APHA Salmonella AMR - final QC" \
    --force

# ---- Attach CheckM2 (not a MultiQC module) ---------------------------------
CHECKM2_REPORT="${CHECKM2_DIR}/quality_report.tsv"
if [[ -s "${CHECKM2_REPORT}" ]]; then
    cp "${CHECKM2_REPORT}" "${MULTIQC_FINAL_DIR}/checkm2_quality_report.tsv"
    echo "  copied CheckM2 quality_report.tsv into report dir"
else
    echo "  WARN: CheckM2 report not found at ${CHECKM2_REPORT}"
fi

echo "MultiQC final report: ${MULTIQC_FINAL_DIR}/multiqc_final_report.html"
echo "step 12 finished on $(date)"

# End of script: 12_multiqc_final.sh

