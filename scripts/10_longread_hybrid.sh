#!/bin/bash
# =============================================================================
# Script       : 10_longread_hybrid.sh
# Description  : Hybrid assembly of the 3 Infantis long-read pairs.
#                Per pair: lightweight ONT read-length summary, Filtlong
#                length/quality filtering, then Unicycler hybrid assembly
#                (trimmed Illumina + filtered ONT). Resolves the pESI
#                megaplasmid that short reads alone left fragmented (step 06).
#                NanoPlot removed (numpy/OpenBLAS vmem explosion on large nodes);
#                read metrics done with awk instead.
#                SGE array job: one task per pair (-t 1-3).
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 00_download.sh (ONT reads), 01a_trim.sh (Illumina) first.
# Usage        : Run from ${PROJECT_DIR} :  qsub scripts/10_longread_hybrid.sh
# Dependencies : conda env apha_longread (filtlong, unicycler, flye)
# Note         : longread_pairs.tsv columns: 1 long_acc, 2 short_acc, 3 biosample
#                ONT reads: RAW_DIR/<long_acc>_1.fastq.gz (single-end)
#                Illumina : TRIM_DIR/<short_acc>_{1,2}.trimmed.fastq.gz
# =============================================================================
#$ -N sc_10_longread_hybrid
#$ -cwd
#$ -M m.a.g.rabbani@roslin.ed.ac.uk
#$ -m as
#$ -l h_vmem=32G
#$ -pe sharedmem 2
#$ -l h_rt=48:00:00
#$ -P roslin_smith_grp
#$ -o logs/
#$ -e logs/
#$ -t 1-3

set -euo pipefail
source config/config.sh
THREADS=2

. /etc/profile.d/modules.sh
module load igmm/apps/anaconda/2023.03
source activate apha_longread

# ---- keep virtual memory sane on large HPC nodes ----------------------------
export MALLOC_ARENA_MAX=2
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

PAIRS="${SAMPLES_DIR}/longread_pairs.tsv"
[[ -f "${PAIRS}" ]] || { echo "ERROR: ${PAIRS} not found"; exit 1; }

# ---- select this task's pair ------------------------------------------------
LINE=$(sed -n "${SGE_TASK_ID}p" "${PAIRS}")
LONG_ACC=$(echo "${LINE}"  | cut -f1)
SHORT_ACC=$(echo "${LINE}" | cut -f2)
BIOSAMPLE=$(echo "${LINE}" | cut -f3)

ONT="${RAW_DIR}/${LONG_ACC}_1.fastq.gz"
R1="${TRIM_DIR}/${SHORT_ACC}_1.trimmed.fastq.gz"
R2="${TRIM_DIR}/${SHORT_ACC}_2.trimmed.fastq.gz"

for f in "${ONT}" "${R1}" "${R2}"; do
    [[ -s "${f}" ]] || { echo "ERROR: missing input ${f}"; exit 1; }
done

WORK="${LONGREAD_DIR}/${SHORT_ACC}"
FILT="${WORK}/${LONG_ACC}.filtlong.fastq.gz"
UNI_DIR="${WORK}/02_unicycler"
QC="${WORK}/01_ont_readstats.txt"
mkdir -p "${WORK}"

echo "Hybrid assembly: ${SHORT_ACC} (long=${LONG_ACC}, biosample=${BIOSAMPLE})"
echo "  task ${SGE_TASK_ID} started on $(date)"

# ---- 1) lightweight ONT read metrics (awk, no numpy) -----------------------
echo "  [1/3] ONT read stats ..."
zcat "${ONT}" | awk '
    NR%4==2 { n++; L=length($0); tot+=L; len[n]=L; if(L>max)max=L }
    END{
        # N50
        asort(len)
        half=tot/2; run=0
        for(i=n;i>=1;i--){ run+=len[i]; if(run>=half){ n50=len[i]; break } }
        printf "raw ONT reads : %d\n", n
        printf "raw ONT bases : %d (%.1f Mb)\n", tot, tot/1e6
        printf "raw ONT N50   : %d\n", n50
        printf "raw ONT max   : %d\n", max
    }' > "${QC}"
cat "${QC}" | sed 's/^/    /'

# ---- 2) Filtlong: keep the best long reads ---------------------------------
echo "  [2/3] Filtlong ..."
filtlong \
    --min_length 1000 \
    --keep_percent 90 \
    --target_bases 500000000 \
    "${ONT}" | gzip > "${FILT}"
echo "    filtered ONT: $(zcat "${FILT}" | awk 'NR%4==2{n++; b+=length($0)} END{print n" reads, "b" bp"}')"

# ---- 3) Unicycler: hybrid assembly -----------------------------------------
echo "  [3/3] Unicycler ..."
unicycler \
    -1 "${R1}" \
    -2 "${R2}" \
    -l "${FILT}" \
    -o "${UNI_DIR}" \
    --threads "${THREADS}" \
    --keep 1

# ---- quick report -----------------------------------------------------------
ASM="${UNI_DIR}/assembly.fasta"
if [[ -s "${ASM}" ]]; then
    echo "  assembly: $(grep -c '>' "${ASM}") contigs"
    echo "  circular contigs:"
    grep '>' "${ASM}" | grep -i "circular=true" | sed 's/^/    /' || echo "    (none flagged circular)"
    cp "${ASM}" "${WORK}/${SHORT_ACC}.hybrid.fasta"
fi

echo "Hybrid assembly: ${SHORT_ACC} completed on $(date)"

# End of script: 10_longread_hybrid.sh

