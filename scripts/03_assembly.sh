#!/bin/bash
# =============================================================================
# Script       : 03_assembly.sh
# Description  : De novo genome assembly of trimmed short-read isolates (Shovill).
#                Runs as an SGE array job (one accession per task).
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 01a_trim.sh first.
# Usage        : Run from ${PROJECT_DIR} :  qsub scripts/03_assembly.sh
# Dependencies : conda env apha_wgs (shovill)
# Last updated : 24/09/2026
# =============================================================================

# =============================================================================
# Grid Engine options
#   -t 1-33      : array job, one task per short-read accession
#   -tc 8        : throttle - max 8 tasks at once
#   -m as        : mail on abort/suspend
# =============================================================================
#$ -N sc_03_assembly
#$ -cwd
#$ -M m.a.g.rabbani@roslin.ed.ac.uk
#$ -m as
#$ -l h_vmem=8G
#$ -pe sharedmem 8
#$ -P roslin_smith_grp
#$ -o logs/
#$ -e logs/
#$ -t 1-33
#$ -tc 8

set -euo pipefail

# ---- Load config ------------------------------------------------------------
source config/config.sh

# ---- Load environment -------------------------------------------------------
. /etc/profile.d/modules.sh
module load igmm/apps/anaconda/2023.03
source activate "${ENV_WGS}"

# ---- Output directory -------------------------------------------------------
mkdir -p "${ASSEMBLY_DIR}"

# ---- Select accession for this array task -----------------------------------
# ${SGE_TASK_ID} is the task number Grid Engine gives each array task (1..33).
# sed prints that line from short_read_ids.txt, so each task handles one isolate.
IDS="${SAMPLES_DIR}/short_read_ids.txt"
ACC=$(sed -n "${SGE_TASK_ID}p" "${IDS}")
[[ -n "${ACC}" ]] || { echo "ERROR: No accession at line ${SGE_TASK_ID}"; exit 1; }

R1="${TRIM_DIR}/${ACC}_1.trimmed.fastq.gz"
R2="${TRIM_DIR}/${ACC}_2.trimmed.fastq.gz"

# ---- Checks -----------------------------------------------------------------
[[ -f "${R1}" && -f "${R2}" ]] || { echo "ERROR: Trimmed pair missing for ${ACC}"; exit 1; }

# ---- Run Shovill ------------------------------------------------------------
# Memory: total = h_vmem x sharedmem (8G x 8 = 64G); --cpus = sharedmem (8);
#         --ram (24) should not be ~80% of total memory - 64G to leave heap headroom ((avoids exit 137 OOM / JVM native alloc failure).).
# --outdir     : per-sample output directory
# --R1 / --R2  : trimmed paired reads
# --gsize      : expected genome size (speeds up read subsampling)
# --minlen 200 : drop contigs shorter than 200 bp (assembly noise)
# --mincov 5   : drop contigs below 5x coverage (junk/contamination)
# --cpus 8     : threads (matches -pe sharedmem 8)
# --ram 24     : memory cap in GB (well below ~80% of total 64G)
# --force      : overwrite existing outdir (safe rerun)


OUT="${ASSEMBLY_DIR}/${ACC}"

echo "Assembly started for ${ACC} as task=${SGE_TASK_ID} on $(date)"

shovill \
    --outdir "${OUT}" \
    --R1 "${R1}" \
    --R2 "${R2}" \
    --gsize 4.8M \
    --minlen 200 \
    --mincov 5 \
    --cpus 8 \
    --ram 24 \
    --force

# ---- Rename final assembly with accession -----------------------------------
# Shovill writes contigs.fa; copy to an accession-named file for downstream steps.
cp "${OUT}/contigs.fa" "${ASSEMBLY_DIR}/${ACC}.fasta"

echo "Assembly completed for ${ACC} as task=${SGE_TASK_ID} on $(date)"

# End of script: 03_assembly.sh
