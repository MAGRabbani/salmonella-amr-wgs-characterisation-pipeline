#!/bin/bash
# =============================================================================
# Script       : 10b_longread_assembly_stat.sh
# Description  : QUAST assembly statistics for the 3 Infantis hybrid assemblies from step 10.
#                Reports contig N50, largest contig, total length and GC (reference-free),
#                complementing the plasmid overview in step 10a.
#                Confirms hybrid-assembly contiguity behind the pESI reconstruction
#                (short reads left pESI fragmented at step 06).
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 10_longread_hybrid.sh first (hybrid assemblies present).
# Usage        : Run from ${PROJECT_DIR} :  bash scripts/10b_longread_assembly_stat.sh
# Dependencies : conda env apha_wgs (quast)
# Note         : Login-node script (QUAST is light; no qsub needed).
#                Hybrid assemblies: LONGREAD_DIR/<short_acc>/02_unicycler/assembly.fasta
#                Pairs come from SAMPLES_DIR/longread_pairs.tsv (col 2 = short_acc)
# =============================================================================

set -uo pipefail
source config/config.sh
THREADS=4

. /etc/profile.d/modules.sh
module load igmm/apps/anaconda/2023.03
source activate apha_wgs

# ---- locate the long-read pairs file ----------------------------------------
PAIRS="${SAMPLES_DIR}/longread_pairs.tsv"
[[ -f "${PAIRS}" ]] || { echo "ERROR: ${PAIRS} not found"; exit 1; }

# ---- output directory -------------------------------------------------------
OUT_DIR="${LONGREAD_DIR}/hybrid_assembly_stats"
mkdir -p "${OUT_DIR}"

# ---- collect hybrid assemblies and labels from the pairs file ---------------
ASSEMBLIES=()
LABELS=""
while IFS=$'\t' read -r LONG_ACC SHORT_ACC BIOSAMPLE; do
    # skip blank lines
    [[ -z "${SHORT_ACC}" ]] && continue
    # final Unicycler assembly for this pair (col 2 = short_acc)
    ASM="${LONGREAD_DIR}/${SHORT_ACC}/02_unicycler/assembly.fasta"
    [[ -s "${ASM}" ]] || { echo "ERROR: assembly not found: ${ASM}"; exit 1; }
    ASSEMBLIES+=("${ASM}")
    LABELS="${LABELS:+${LABELS},}${SHORT_ACC}"
done < "${PAIRS}"

# ---- fail if nothing was collected ------------------------------------------
[[ ${#ASSEMBLIES[@]} -gt 0 ]] || { echo "ERROR: no assemblies found"; exit 1; }

# ---- run QUAST (reference-free: contiguity stats only) ----------------------
echo "Running QUAST on ${#ASSEMBLIES[@]} hybrid assemblies -> ${OUT_DIR}"
quast.py \
    "${ASSEMBLIES[@]}" \
    -l "${LABELS}" \
    --threads "${THREADS}" \
    -o "${OUT_DIR}"

echo "Done. Report: ${OUT_DIR}/report.tsv"

# End of script: 10b_longread_assembly_stat.sh

