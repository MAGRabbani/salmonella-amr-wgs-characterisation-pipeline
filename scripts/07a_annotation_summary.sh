#!/bin/bash
# =============================================================================
# Script       : 07a_annotation_summary.sh
# Description  : Summarise Bakta annotation across all isolates.
#                Per-isolate feature table (CDS, tRNA, rRNA, ncRNA, hypothetical,
#                coding density, genome size) parsed from each <ACC>.txt, serovar
#                annotated and sorted. Flags isolates outside expected CDS range
#                as a pre-pangenome QC check.
#                Light text parsing; run on the login node with bash.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 07_annotation.sh and 04a_serotype_mlst_typing_summary.sh first.
# Usage        : Run from ${PROJECT_DIR} :  bash scripts/07a_annotation_summary.sh
# Note         : typing_summary.tsv col 5 = SISTR_serovar.
# =============================================================================
set -euo pipefail
source config/config.sh

TYPING_TSV="${SEROTYPE_DIR}/typing_summary.tsv"
OUT_TABLE="${ANNOT_DIR}/01_annotation_summary.tsv"

echo "Building annotation summary on $(date)"

# ---- serovar lookup ---------------------------------------------------------
declare -A SERO
while IFS=$'\t' read -r iso c2 c3 c4 sero rest; do
    [[ "${iso}" == "isolate" ]] && continue
    SERO["${iso}"]="${sero}"
done < "${TYPING_TSV}"

# ---- per-isolate feature table ---------------------------------------------
printf "isolate\tserovar\tsize_Mb\tGC\tcoding_density\tCDS\ttRNA\trRNA\tncRNA\thypothetical\tflag\n" > "${OUT_TABLE}"

for D in "${ANNOT_DIR}"/*/; do
    ACC=$(basename "${D}")
    TXT="${D}${ACC}.txt"
    [[ -f "${TXT}" ]] || continue
    s="${SERO[${ACC}]:-NA}"

    len=$(awk -F': ' '/^Length:/{print $2}' "${TXT}")
    gc=$(awk -F': '  '/^GC:/{print $2}'     "${TXT}")
    cod=$(awk -F': ' '/^coding density:/{print $2}' "${TXT}")
    cds=$(awk -F': ' '/^CDSs:/{print $2}'   "${TXT}")
    trna=$(awk -F': ' '/^tRNAs:/{print $2}' "${TXT}")
    rrna=$(awk -F': ' '/^rRNAs:/{print $2}' "${TXT}")
    ncrna=$(awk -F': ' '/^ncRNAs:/{print $2}' "${TXT}")
    hyp=$(awk -F': ' '/^hypotheticals:/{print $2}' "${TXT}")

    mb=$(awk -v l="${len}" 'BEGIN{printf "%.2f", l/1000000}')

    # QC flag: expected Salmonella CDS ~4400-5000
    flag="OK"
    if [[ "${cds}" -lt 4200 || "${cds}" -gt 5000 ]]; then flag="CHECK_CDS"; fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${ACC}" "${s}" "${mb}" "${gc}" "${cod}" "${cds}" "${trna}" "${rrna}" "${ncrna}" "${hyp}" "${flag}" \
        >> "${OUT_TABLE}"
done

# ---- sort by serovar then isolate ------------------------------------------
(head -1 "${OUT_TABLE}"; tail -n +2 "${OUT_TABLE}" | sort -k2,2 -k1,1) > "${OUT_TABLE}.tmp" \
    && mv "${OUT_TABLE}.tmp" "${OUT_TABLE}"

# ---- screen summary ---------------------------------------------------------
N=$(($(wc -l < "${OUT_TABLE}") - 1))
N_FLAG=$(awk -F'\t' 'NR>1 && $11!="OK"' "${OUT_TABLE}" | wc -l)

echo ""
echo "=============================================="
echo " Annotation summary (Bakta)"
echo "=============================================="
echo " Isolates annotated : ${N}"
echo " Flagged (CDS)      : ${N_FLAG}"
echo ""
echo " Mean CDS per serovar:"
awk -F'\t' 'NR>1{sum[$2]+=$6; n[$2]++} END{for(s in n) printf "   %-14s %.0f  (n=%d)\n", s, sum[s]/n[s], n[s]}' "${OUT_TABLE}" | sort
echo ""
echo " Full table: ${OUT_TABLE}"
echo "=============================================="
echo ""
echo "Annotation summary finished on $(date)"

# End of script: 07a_annotation_summary.sh

