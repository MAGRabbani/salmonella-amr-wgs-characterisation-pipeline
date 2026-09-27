#!/bin/bash
# =============================================================================
# Script       : 00a_download_raw_reads.sh
# Description  : Download sample FASTQ files from ENA using accessions listed
#                in samplesheet.tsv. Handles paired and single-end reads.
#                wget -c resumes partial downloads and skips completed ones.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Usage        : cd "${PROJECT_DIR}"
#                bash scripts/00a_download_raw_reads.sh
# =============================================================================
set -euo pipefail
source config/config.sh

mkdir -p "${RAW_DIR}" "${LOG_DIR}"
LOG="${LOG_DIR}/00_download.log"
echo "Download started: $(date)" | tee -a "${LOG}"

ACC_LIST=$(awk -F'\t' 'NR>1 {print $2}' "${SAMPLE_SHEET}")

for ACC in ${ACC_LIST}; do
    echo "Processing ${ACC}" | tee -a "${LOG}"

    URLS=$(wget -qO- "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=fastq_ftp&format=tsv" \
           | tail -n +2 | cut -f2)

    if [[ -z "${URLS}" ]]; then
        echo "  No FASTQ found for ${ACC}" | tee -a "${LOG}"
        continue
    fi

    for f in ${URLS//;/ }; do
        echo "  Fetching: $(basename "${f}")" | tee -a "${LOG}"
        wget -c --timeout=60 --tries=5 --waitretry=10 -P "${RAW_DIR}" "https://${f}"
    done
done

echo "Download finished: $(date)" | tee -a "${LOG}"
echo "Total files: $(ls "${RAW_DIR}"/*.fastq.gz 2>/dev/null | wc -l)" | tee -a "${LOG}"
# End of script: 00a_download_raw_reads.sh
