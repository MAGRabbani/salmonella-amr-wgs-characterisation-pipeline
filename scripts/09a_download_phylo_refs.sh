#!/bin/bash
# =============================================================================
# Script       : 09a_download_phylo_refs.sh
# Description  : Fetch per-serovar reference genomes for within-serovar
#                phylogeny (snippy mapping). Downloads Enteritidis (P125109),
#                Infantis (Z1323CSL0027, pESI+) and Kentucky (PU131) from ENA;
#                links the existing Typhimurium LT2 reference into place.
#                Run once, interactively on the login node (needs internet).
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Usage        : Run from ${PROJECT_DIR} :  bash scripts/09a_download_phylo_refs.sh
# =============================================================================
set -euo pipefail
source config/config.sh

echo "Downloading per-serovar phylogeny references on $(date)"

# ---- accession -> target path ----------------------------------------------
declare -A REFS=(
    ["Enteritidis"]="GCA_000009505.1"
    ["Infantis"]="GCA_037776255.1"
    ["Kentucky"]="GCA_002952975.1"
)

# ENA assembly FASTA download helper (uses the ENA browser API 'fasta' export)
fetch_ena_fasta () {
    local acc="$1" out="$2"
    # ENA serves assembly sequence FASTA via the browser API
    local url="https://www.ebi.ac.uk/ena/browser/api/fasta/${acc}?download=true&gzip=true"
    echo "  fetching ${acc} ..."
    wget -q -c -O "${out}.gz" "${url}"
    gunzip -f "${out}.gz"
}

for sero in "${!REFS[@]}"; do
    acc="${REFS[$sero]}"
    dir="${PHYLO_REF_DIR}/${sero}"
    mkdir -p "${dir}"
    out="${dir}/${acc}.fasta"
    if [[ -s "${out}" ]]; then
        echo "  ${sero}: ${acc} already present, skipping"
    else
        fetch_ena_fasta "${acc}" "${out}"
    fi
    echo "  ${sero}: $(grep -c '>' "${out}") sequence(s), $(grep -v '>' "${out}" | tr -d '\n' | wc -c) bp"
done

# ---- Typhimurium: link the existing LT2 reference --------------------------
mkdir -p "${PHYLO_REF_DIR}/Typhimurium"
if [[ ! -s "${REF_TYPHIMURIUM}" ]]; then
    if [[ -s "${REF_GENOME}" ]]; then
        cp "${REF_GENOME}" "${REF_TYPHIMURIUM}"
        echo "  Typhimurium: copied existing LT2 from ${REF_GENOME}"
    else
        echo "  Typhimurium: fetching LT2 GCA_000006945.2 fresh"
        fetch_ena_fasta "GCA_000006945.2" "${REF_TYPHIMURIUM}"
    fi
fi
echo "  Typhimurium: $(grep -c '>' "${REF_TYPHIMURIUM}") sequence(s)"

echo ""
echo "References ready under ${PHYLO_REF_DIR}"
echo "Download finished on $(date)"

# End of script: 09a_download_phylo_refs.sh

