#!/bin/bash
# =============================================================================
# Script       : 09b_phylo_snippy.sh
# Description  : Within-serovar SNP phylogeny (APHA outbreak-style relatedness).
#                Per serovar: snippy maps each isolate's trimmed reads to a
#                serovar-specific reference, snippy-core builds a core-SNP
#                alignment, Gubbins removes recombination, IQ-TREE builds the
#                recombination-free ML tree.
#                SGE array job: one task per serovar (-t 1-4).
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 01a_trim.sh, 04a typing, 09a_download_phylo_refs.sh first.
# Usage        : Run from ${PROJECT_DIR} :  qsub scripts/09b_phylo_snippy.sh
# Dependencies : conda env apha_phylo (snippy, snp-sites, gubbins, iqtree)
# Note         : Reads named <ACC>_1.trimmed.fastq.gz / _2 in TRIM_DIR.
#                Serovar lists derived from typing_summary.tsv (col 5).
# =============================================================================
#$ -N sc_09b_phylo_snippy
#$ -cwd
#$ -M m.a.g.rabbani@roslin.ed.ac.uk
#$ -m as
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -P roslin_smith_grp
#$ -o logs/
#$ -e logs/
#$ -t 1-4

set -euo pipefail
source config/config.sh
THREADS=4

. /etc/profile.d/modules.sh
module load igmm/apps/anaconda/2023.03
source activate apha_phylo

TYPING_TSV="${SEROTYPE_DIR}/typing_summary.tsv"

# ---- map array task -> serovar + reference ---------------------------------
SEROVARS=(Typhimurium Enteritidis Infantis Kentucky)
REFS=("${REF_TYPHIMURIUM}" "${REF_ENTERITIDIS}" "${REF_INFANTIS}" "${REF_KENTUCKY}")

IDX=$(( SGE_TASK_ID - 1 ))
SEROVAR="${SEROVARS[$IDX]}"
REF="${REFS[$IDX]}"

[[ -s "${REF}" ]] || { echo "ERROR: reference missing for ${SEROVAR}: ${REF}"; exit 1; }

WORK="${PHYLO_DIR}/02_snippy/${SEROVAR}"
SNIPPY_OUT="${WORK}/snippy"
mkdir -p "${SNIPPY_OUT}"

echo "Snippy phylogeny: ${SEROVAR} (task ${SGE_TASK_ID}) started on $(date)"

# ---- isolate list for this serovar (from typing col 5) ---------------------
mapfile -t ISOLATES < <(awk -F'\t' -v s="${SEROVAR}" 'NR>1 && $5==s {print $1}' "${TYPING_TSV}")
echo "  ${SEROVAR}: ${#ISOLATES[@]} isolates"
[[ "${#ISOLATES[@]}" -ge 3 ]] || { echo "ERROR: <3 isolates for ${SEROVAR}, cannot build tree"; exit 1; }

# ---- 1) snippy per isolate --------------------------------------------------
for ACC in "${ISOLATES[@]}"; do
    R1="${TRIM_DIR}/${ACC}_1.trimmed.fastq.gz"
    R2="${TRIM_DIR}/${ACC}_2.trimmed.fastq.gz"
    [[ -s "${R1}" && -s "${R2}" ]] || { echo "  WARN: reads missing for ${ACC}, skipping"; continue; }

    if [[ -f "${SNIPPY_OUT}/${ACC}/snps.aligned.fa" ]]; then
        echo "  snippy ${ACC}: already done"
        continue
    fi
    echo "  snippy ${ACC} ..."
    snippy \
        --outdir "${SNIPPY_OUT}/${ACC}" \
        --ref "${REF}" \
        --R1 "${R1}" \
        --R2 "${R2}" \
        --cpus "${THREADS}" \
        --force
done

# ---- 2) snippy-core: core-SNP alignment ------------------------------------
cd "${SNIPPY_OUT}"
snippy-core \
    --ref "${REF}" \
    --prefix "${WORK}/core" \
    "${SNIPPY_OUT}"/*/
cd "${PROJECT_DIR}"

# snippy-core writes core.full.aln (whole-genome) and core.aln (SNP-only)
FULL_ALN="${WORK}/core.full.aln"
[[ -s "${FULL_ALN}" ]] || { echo "ERROR: snippy-core produced no alignment for ${SEROVAR}"; exit 1; }

# ---- 3) clean alignment for Gubbins (snippy-clean-full-aln) -----------------
# Gubbins dislikes non-ACGTN chars; snippy ships a cleaner.
CLEAN_ALN="${WORK}/core.full.clean.aln"
snippy-clean_full_aln "${FULL_ALN}" > "${CLEAN_ALN}"

# ---- 4) Gubbins: remove recombination --------------------------------------
GUB_DIR="${WORK}/gubbins"
mkdir -p "${GUB_DIR}"
cd "${GUB_DIR}"
run_gubbins.py \
    --prefix "${SEROVAR}" \
    --threads "${THREADS}" \
    "${CLEAN_ALN}" || {
        echo "  WARN: Gubbins failed for ${SEROVAR} (often too few SNPs / isolates);"
        echo "        falling back to SNP tree without recombination removal."
    }
cd "${PROJECT_DIR}"

# ---- 5) IQ-TREE on the recombination-filtered SNP alignment -----------------
# Prefer Gubbins filtered polymorphic sites; fall back to snippy core SNPs.
GUB_SNP="${GUB_DIR}/${SEROVAR}.filtered_polymorphic_sites.fasta"
TREE_OUT="${WORK}/iqtree"
mkdir -p "${TREE_OUT}"

if [[ -s "${GUB_SNP}" ]]; then
    echo "  building tree from Gubbins-filtered SNPs"
    TREE_IN="${GUB_SNP}"
    EXTRA=""
else
    echo "  building tree from snippy core SNPs (no recombination removal)"
    TREE_IN="${WORK}/core.aln"
    EXTRA="+ASC"   # ascertainment bias correction for SNP-only alignment
fi

iqtree2 \
    -s "${TREE_IN}" \
    -m MFP \
    -B 1000 \
    -T "${THREADS}" \
    --prefix "${TREE_OUT}/${SEROVAR}_snp" \
    -redo

echo "Snippy phylogeny: ${SEROVAR} completed on $(date)"
echo "  core alignment : ${FULL_ALN}"
echo "  SNP tree       : ${TREE_OUT}/${SEROVAR}_snp.treefile"

# End of script: 09b_phylo_snippy.sh

