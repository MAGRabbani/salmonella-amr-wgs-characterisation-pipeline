#!/bin/bash
# =============================================================================
# Script       : 00b_download_ref_genomes.sh
# Description  : Download S. Typhimurium LT2 reference genome (GCA_000006945.2)
#                from ENA. Used as reference for the phylogeny step.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Usage        : bash scripts/00b_download_ref_genomes.sh
# =============================================================================
set -euo pipefail
source config/config.sh

mkdir -p "${REF_DIR}"

wget -O "${REF_DIR}/reference.fna.gz" \
    "https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_000006945.2?download=true&gzip=true"

gunzip -f "${REF_DIR}/reference.fna.gz"

echo "Reference sequences: $(grep -c '>' "${REF_GENOME}")"
# End of script: 00b_download_ref_genomes.sh
