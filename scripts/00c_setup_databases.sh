#!/bin/bash
# =============================================================================
# Script       : 00c_setup_databases.sh
# Description  : One-time download of reference databases used by the pipeline.
#                Run interactively on the LOGIN NODE (needs internet).
#                Switches conda env per tool.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Usage        : Run from ${PROJECT_DIR} :  bash scripts/00c_setup_databases.sh
# =============================================================================
set -euo pipefail
source config/config.sh

# ---- Initialise conda -------------------------------------------------------
# Source personal conda directly (NOT the anaconda module, which overrides the
# env's Python and breaks Python-import tools like BUSCO). Path from config.
source "${CONDA_BASE}/etc/profile.d/conda.sh"

echo "=== Database setup started: $(date) ==="
mkdir -p "${REF_DIR}"

# ---- BUSCO lineage (env: apha_wgs) -----------------------------------------
# Enterobacterales lineage for Salmonella; downloaded once for offline use in 03a.
# busco writes to ./busco_downloads/lineages/, so cd into REF_DIR first.
echo "--- BUSCO lineage: enterobacterales_odb10 ---"
conda activate "${ENV_WGS}"
cd "${REF_DIR}"
busco --download enterobacterales_odb10
cd "${PROJECT_DIR}"

# ---- CheckM2 database (env: checkm2) ----------------------------------------
echo "--- CheckM2 DB ---"
conda deactivate || true
conda activate "${ENV_CHECKM2}"
checkm2 database --download --path "${CHECKM2_DB}"

# ---- AMRFinderPlus + ABRicate (env: apha_amrPlasfinder) --------------------
echo "--- AMRFinderPlus DB ---"
conda deactivate || true
conda activate "${ENV_PLASMID}"
amrfinder -u

echo "--- ABRicate DBs ---"
abricate --setupdb
abricate --list

# ---- Bakta database (env: apha_annotation) ---------------------------------
# Light DB (~1.3 GB). Use --type full for the complete DB (~30 GB) if needed.
# Env name kept literal here (not in config), matching 07_annotation.sh.
echo "--- Bakta DB (light) ---"
conda deactivate || true
conda activate apha_annotation
mkdir -p "${REF_DIR}/bakta_db"
bakta_db download --output "${REF_DIR}/bakta_db" --type light

echo "=== Database setup finished: $(date) ==="

# End of script: 00c_setup_databases.sh

