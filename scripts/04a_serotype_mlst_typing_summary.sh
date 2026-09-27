#!/bin/bash
# =============================================================================
# Script       : 04a_serotype_mlst_typing_summary.sh
# Description  : Join MLST + SISTR + SeqSero2 into one per-isolate typing table.
#                Captures classical 7-gene ST and cgMLST ST (high-resolution,
#                as used by UKHSA/EnteroBase), subspecies, dual-method serovar
#                with agreement, antigenic profile and SISTR QC.
#                Light text parsing; can run on the login node with bash.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 04_serotype_mlst.sh first.
# Usage        : Run from ${PROJECT_DIR} :  bash scripts/04a_serotype_mlst_typing_summary.sh
# =============================================================================
set -euo pipefail
source config/config.sh

MLST_DIR="${SEROTYPE_DIR}/01_mlst"
SISTR_DIR="${SEROTYPE_DIR}/02_sistr"
SEQSERO_DIR="${SEROTYPE_DIR}/03_seqsero2"
OUT_TSV="${SEROTYPE_DIR}/typing_summary.tsv"

echo "Building typing summary on $(date)"

# ---- Helper: pull a SISTR column value by header name -----------------------
sistr_val() {
    # $1 = sistr file, $2 = column header
    awk -F'\t' -v col="$2" '
        NR==1{for(i=1;i<=NF;i++) if($i==col) c=i}
        NR==2{print $c}' "$1"
}

# ---- Header -----------------------------------------------------------------
printf "isolate\tMLST_ST\tcgMLST_ST\tsubspecies\tSISTR_serovar\tSISTR_cgmlst\tSISTR_antigen\tSISTR_QC\tSeqSero2_serovar\tantigenic_profile\tagreement\n" > "${OUT_TSV}"

# ---- Loop over assemblies ---------------------------------------------------
for ASM in "${ASSEMBLY_DIR}"/*.fasta; do
    ACC=$(basename "${ASM}" .fasta)

    # --- MLST: ST is field 3 ---
    MLST_FILE="${MLST_DIR}/${ACC}.mlst.tab"
    if [[ -f "${MLST_FILE}" ]]; then
        ST=$(awk -F'\t' '{print $3}' "${MLST_FILE}")
    else
        ST="NA"
    fi

    # --- SISTR: cgMLST ST, subspecies, serovar, cgmlst serovar, antigen, QC ---
    SISTR_FILE="${SISTR_DIR}/${ACC}.sistr.tab"
    if [[ -f "${SISTR_FILE}" ]]; then
        CG_ST=$(sistr_val "${SISTR_FILE}" "cgmlst_ST")
        SUBSP=$(sistr_val "${SISTR_FILE}" "cgmlst_subspecies")
        S_SEROVAR=$(sistr_val "${SISTR_FILE}" "serovar")
        S_CGMLST=$(sistr_val "${SISTR_FILE}" "serovar_cgmlst")
        S_ANTIGEN=$(sistr_val "${SISTR_FILE}" "serovar_antigen")
        S_QC=$(sistr_val "${SISTR_FILE}" "qc_status")
    else
        CG_ST="NA"; SUBSP="NA"; S_SEROVAR="NA"; S_CGMLST="NA"; S_ANTIGEN="NA"; S_QC="NA"
    fi

    # --- SeqSero2: predicted serotype + antigenic profile ---
    SS_FILE="${SEQSERO_DIR}/${ACC}/SeqSero_result.txt"
    if [[ -f "${SS_FILE}" ]]; then
        SS_SEROVAR=$(grep "Predicted serotype:" "${SS_FILE}" | sed 's/Predicted serotype:[[:space:]]*//')
        SS_PROFILE=$(grep "Predicted antigenic profile:" "${SS_FILE}" | sed 's/Predicted antigenic profile:[[:space:]]*//')
    else
        SS_SEROVAR="NA"; SS_PROFILE="NA"
    fi

    # --- Agreement (SISTR serovar vs SeqSero2 serovar) ---
    if [[ "${S_SEROVAR}" == "${SS_SEROVAR}" && "${S_SEROVAR}" != "NA" ]]; then
        AGREE="YES"
    else
        AGREE="CHECK"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${ACC}" "${ST}" "${CG_ST}" "${SUBSP}" "${S_SEROVAR}" "${S_CGMLST}" \
        "${S_ANTIGEN}" "${S_QC}" "${SS_SEROVAR}" "${SS_PROFILE}" "${AGREE}" >> "${OUT_TSV}"
done

# ---- Short summary ----------------------------------------------------------
# Column positions: 1 isolate, 2 MLST_ST, 3 cgMLST_ST, 4 subspecies,
# 5 SISTR_serovar, 8 SISTR_QC, 11 agreement
TOTAL=$(($(wc -l < "${OUT_TSV}") - 1))
AGREE_N=$(awk -F'\t' 'NR>1 && $11=="YES"' "${OUT_TSV}" | wc -l)
QC_PASS=$(awk -F'\t' 'NR>1 && $8=="PASS"' "${OUT_TSV}" | wc -l)

echo ""
echo "=============================================="
echo " Typing summary"
echo "=============================================="
echo " Table               : ${OUT_TSV}"
echo " Total isolates      : ${TOTAL}"
echo " SISTR/SeqSero2 agree: ${AGREE_N}/${TOTAL}"
echo " SISTR QC PASS       : ${QC_PASS}/${TOTAL}"
echo ""
echo " Serovar counts (SISTR):"
awk -F'\t' 'NR>1{print $5}' "${OUT_TSV}" | sort | uniq -c | sort -rn | sed 's/^/   /'
echo ""
echo " 7-gene ST counts (MLST):"
awk -F'\t' 'NR>1{print $2}' "${OUT_TSV}" | sort | uniq -c | sort -rn | sed 's/^/   /'
echo ""
echo " Subspecies:"
awk -F'\t' 'NR>1{print $4}' "${OUT_TSV}" | sort | uniq -c | sort -rn | sed 's/^/   /'
echo "=============================================="

if command -v column >/dev/null 2>&1; then
    echo ""
    column -t -s $'\t' "${OUT_TSV}"
fi

echo ""
echo "Typing summary finished on $(date)"

# End of script: 04a_serotype_mlst_typing_summary.sh

