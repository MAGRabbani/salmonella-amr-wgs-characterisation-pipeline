#!/bin/bash
# =============================================================================
# Script       : 03a_assembly_qc.sh
# Description  : Assembly QC of all isolate assemblies.
#                QUAST (metrics) + BUSCO (gene completeness, offline) +
#                CheckM2 (completeness/contamination). Single (non-array) job.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 03_assembly.sh and setup_databases.sh first.
# Usage        : Run from ${PROJECT_DIR} :  qsub scripts/03a_assembly_qc.sh
# Dependencies : apha_wgs (quast, busco); apha_checkm2 (checkm2)
# Last updated : 25/09/2026
# =============================================================================

# =============================================================================
# Grid Engine options
# =============================================================================
#$ -N sc_03a_assembly_qc
#$ -cwd
#$ -M m.a.g.rabbani@roslin.ed.ac.uk
#$ -m beas
#$ -l h_vmem=16G
#$ -pe sharedmem 4
#$ -P roslin_smith_grp
#$ -o logs/
#$ -e logs/

set -euo pipefail

# ---- Load config ------------------------------------------------------------
source config/config.sh

# This job runs on 4 slots; override THREADS for tool commands below.
THREADS=4

# ---- Load environment -------------------------------------------------------
. /etc/profile.d/modules.sh
module load igmm/apps/anaconda/2023.03
source activate "${ENV_WGS}"


# ---- Output directories -----------------------------------------------------
mkdir -p "${QUAST_DIR}" "${BUSCO_DIR}" "${CHECKM2_DIR}"

# ---- Checks -----------------------------------------------------------------
N_ASM=$(ls "${ASSEMBLY_DIR}"/*.fasta 2>/dev/null | wc -l)
[[ "${N_ASM}" -gt 0 ]] || { echo "ERROR: No assemblies in ${ASSEMBLY_DIR}"; exit 1; }
echo "Found ${N_ASM} assemblies for QC on $(date)"

# =============================================================================
# 1) QUAST  (env: apha_wgs)  - assembly metrics on all isolates at once
# =============================================================================
conda activate "${ENV_WGS}"
echo "QUAST started on $(date)"

quast.py "${ASSEMBLY_DIR}"/*.fasta \
    -o "${QUAST_DIR}" \
    --threads "${THREADS}" \
    --min-contig 200

echo "QUAST completed on $(date)"

# =============================================================================
# 2) BUSCO  (env: apha_wgs)  - gene completeness vs Enterobacterales (offline)
# =============================================================================
echo "BUSCO started on $(date)"
cd "${BUSCO_DIR}"
for ASM in "${ASSEMBLY_DIR}"/*.fasta; do
    NAME=$(basename "${ASM}" .fasta)
    busco \
        -i "${ASM}" \
        -o "${NAME}" \
        -l "${BUSCO_LINEAGE}" \
        -m genome \
        -c "${THREADS}" \
        --offline \
        -f
done
cd "${PROJECT_DIR}"
echo "BUSCO completed on $(date)"

# =============================================================================
# 3) CheckM2  (env: apha_checkm2)  - completeness & contamination
# =============================================================================
conda deactivate || true
conda activate "${ENV_CHECKM2}"
export CHECKM2DB="${CHECKM2_DB_FILE}"

echo "CheckM2 started on $(date)"

checkm2 predict \
    --input "${ASSEMBLY_DIR}" \
    --extension .fasta \
    --output-directory "${CHECKM2_DIR}" \
    --threads "${THREADS}" \
	--lowmem \
    --force

echo "CheckM2 completed on $(date)"
echo "Assembly QC finished on $(date)"

# End of script: 03a_assembly_qc.sh
