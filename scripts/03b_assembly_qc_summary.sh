#!/bin/bash
# =============================================================================
# Script       : 03b_qc_summary.sh
# Description  : Join QUAST + BUSCO + CheckM2 (+coverage) into one per-isolate
#                QC table, with pass/fail verdict against thresholds.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 03a_assembly_qc.sh first.
# Usage        : Run from ${PROJECT_DIR} :  bash scripts/03b_qc_summary.sh
# =============================================================================
set -euo pipefail
source config/config.sh

OUT_TSV="${ASSEMBLY_QC_DIR}/assembly_qc_summary.tsv"
QUAST_REPORT="${QUAST_DIR}/transposed_report.tsv"
CHECKM2_REPORT="${CHECKM2_DIR}/quality_report.tsv"
COVERAGE_FILE="${PROJECT_DIR}/coverage_check.tsv"

# ---- Checks -----------------------------------------------------------------
[[ -f "${QUAST_REPORT}"   ]] || { echo "ERROR: missing ${QUAST_REPORT}"; exit 1; }
[[ -f "${CHECKM2_REPORT}" ]] || { echo "ERROR: missing ${CHECKM2_REPORT}"; exit 1; }

echo "Building QC summary on $(date)"

# ---- Thresholds -------------------------------------------------------------
MIN_COMPLETENESS=95
MAX_CONTAMINATION=5
MIN_N50=50000
MIN_SIZE=4500000
MAX_SIZE=5200000
MAX_CONTIGS=500
MIN_COVERAGE=30

# ---- Helper: pull a QUAST column value by header name for a given isolate ---
quast_val() {
    # $1 = isolate, $2 = exact column header
    awk -F'\t' -v s="$1" -v col="$2" '
        NR==1{for(i=1;i<=NF;i++) if($i==col) c=i}
        $1==s{print $c}' "${QUAST_REPORT}"
}

# ---- Header -----------------------------------------------------------------
printf "isolate\tcontigs\tN50\tL50\tlargest\tsize_bp\tGC\tNper100kb\tcoverage\tBUSCO_C\tCheckM2_Comp\tCheckM2_Contam\tverdict\n" > "${OUT_TSV}"

# ---- Loop over assemblies ---------------------------------------------------
for ASM in "${ASSEMBLY_DIR}"/*.fasta; do
    ACC=$(basename "${ASM}" .fasta)

    CONTIGS=$(quast_val "${ACC}" "# contigs")
    N50=$(quast_val "${ACC}" "N50")
    L50=$(quast_val "${ACC}" "L50")
    LARGEST=$(quast_val "${ACC}" "Largest contig")
    SIZE=$(quast_val "${ACC}" "Total length")
    GC=$(quast_val "${ACC}" "GC (%)")
    NPER=$(quast_val "${ACC}" "# N's per 100 kbp")
    [[ -z "${NPER}" ]] && NPER="0.00"

    # BUSCO %C
    BUSCO_FILE=$(ls "${BUSCO_DIR}/${ACC}"/short_summary.*."${ACC}".txt 2>/dev/null | head -1)
    if [[ -n "${BUSCO_FILE}" ]]; then
        BUSCO_C=$(grep -oP 'C:\K[0-9.]+' "${BUSCO_FILE}" | head -1)
    else
        BUSCO_C="NA"
    fi

    # CheckM2
    COMP=$(awk -F'\t' -v s="${ACC}" '$1==s{print $2}' "${CHECKM2_REPORT}")
    CONTAM=$(awk -F'\t' -v s="${ACC}" '$1==s{print $3}' "${CHECKM2_REPORT}")

    # Coverage (optional)
    if [[ -f "${COVERAGE_FILE}" ]]; then
        COV=$(awk -F'\t' -v s="${ACC}" '$1==s{print $2}' "${COVERAGE_FILE}")
        [[ -z "${COV}" ]] && COV="NA"
    else
        COV="NA"
    fi

    # ---- Verdict ----
    VERDICT="PASS"
    awk -v c="${COMP}"    -v t="${MIN_COMPLETENESS}"  'BEGIN{exit !(c+0>=t)}' || VERDICT="FAIL"
    awk -v c="${CONTAM}"  -v t="${MAX_CONTAMINATION}" 'BEGIN{exit !(c+0<=t)}' || VERDICT="FAIL"
    awk -v c="${N50}"     -v t="${MIN_N50}"           'BEGIN{exit !(c+0>=t)}' || VERDICT="FAIL"
    awk -v c="${SIZE}"    -v t="${MIN_SIZE}"          'BEGIN{exit !(c+0>=t)}' || VERDICT="FAIL"
    awk -v c="${SIZE}"    -v t="${MAX_SIZE}"          'BEGIN{exit !(c+0<=t)}' || VERDICT="FAIL"
    awk -v c="${CONTIGS}" -v t="${MAX_CONTIGS}"       'BEGIN{exit !(c+0<=t)}' || VERDICT="FAIL"
    if [[ "${COV}" != "NA" ]]; then
        awk -v c="${COV}" -v t="${MIN_COVERAGE}"      'BEGIN{exit !(c+0>=t)}' || VERDICT="FAIL"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${ACC}" "${CONTIGS}" "${N50}" "${L50}" "${LARGEST}" "${SIZE}" \
        "${GC}" "${NPER}" "${COV}" "${BUSCO_C}" "${COMP}" "${CONTAM}" "${VERDICT}" >> "${OUT_TSV}"
done

# ---- Short summary ----------------------------------------------------------
TOTAL=$(($(wc -l < "${OUT_TSV}") - 1))
PASS=$(awk -F'\t' 'NR>1 && $13=="PASS"' "${OUT_TSV}" | wc -l)
FAIL=$((TOTAL - PASS))

echo ""
echo "=============================================="
echo " Assembly QC summary"
echo "=============================================="
echo " Table          : ${OUT_TSV}"
echo " Total isolates : ${TOTAL}"
echo " PASS           : ${PASS}"
echo " FAIL           : ${FAIL}"
echo ""
echo " Thresholds: Comp>=${MIN_COMPLETENESS}%, Contam<=${MAX_CONTAMINATION}%, N50>=${MIN_N50},"
echo "             Size ${MIN_SIZE}-${MAX_SIZE} bp, Contigs<=${MAX_CONTIGS}, Coverage>=${MIN_COVERAGE}x"
echo "=============================================="

if command -v column >/dev/null 2>&1; then
    echo ""
    column -t -s $'\t' "${OUT_TSV}"
fi

echo ""
echo "QC summary finished on $(date)"

# End of script: 03b_qc_summary.sh

