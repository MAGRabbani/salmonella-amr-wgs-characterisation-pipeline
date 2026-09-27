#!/bin/bash
# =============================================================================
# Script       : 10a_longread_summary.sh
# Description  : Summarise hybrid assemblies (any number of isolate pairs).
#                Plasmid calls come from MOB-suite, which groups fragmented
#                contigs into their true plasmid, so no contig-size heuristics
#                are used and the summary is reproducible for any sample / any
#                plasmid. Two outputs (Option A: overview + detail):
#                  (1) 01_assembly_overview.tsv : one row per isolate (aggregates)
#                  (2) 02_plasmids_detail.tsv   : one row per plasmid
#                All bp columns are plain integers (sortable / machine-readable);
#                the remark column carries human-readable kb.
#                Runs in the plasmid env (apha_amrPlasfinder); self-activating.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 10_longread_hybrid.sh first.
# Usage        : Run from ${PROJECT_DIR} :  bash scripts/10a_longread_summary.sh
# Columns (overview):
#   isolate           short-read (Illumina) accession = sample name
#   long_acc          matched Nanopore accession used in the hybrid
#   n_contigs         number of contigs (assembly fragmentation; fewer = better)
#   assembly_size_bp  whole assembly size (chromosome + all plasmids)
#   n_plasmids        number of plasmids reconstructed by MOB-suite
#   total_plasmid_bp  combined bp of ALL plasmids (per-plasmid sizes in detail)
#   largest_replicon  replicon type of the largest plasmid (e.g. IncFIB)
#   remark            plain-language interpretation, auto-generated from data
# =============================================================================
set -uo pipefail
source config/config.sh

# ---- activate the plasmid env (MOB-suite) ----------------------------------
. /etc/profile.d/modules.sh 2>/dev/null || true
module load igmm/apps/anaconda/2023.03 2>/dev/null || true
source activate "${ENV_PLASMID}" 2>/dev/null || conda activate "${ENV_PLASMID}"

PAIRS="${SAMPLES_DIR}/longread_pairs.tsv"
OUT_ISO="${LONGREAD_DIR}/01_assembly_overview.tsv"
OUT_PLASMID="${LONGREAD_DIR}/02_plasmids_detail.tsv"

echo "Building hybrid-assembly summary on $(date)"

printf "isolate\tlong_acc\tn_contigs\tassembly_size_bp\tn_plasmids\ttotal_plasmid_bp\tlargest_replicon\tremark\n" > "${OUT_ISO}"
printf "isolate\tplasmid_cluster\tsize_bp\trep_type\tmobility\n" > "${OUT_PLASMID}"

while IFS=$'\t' read -r LONG_ACC SHORT_ACC BIOSAMPLE; do
    [[ -z "${SHORT_ACC}" ]] && continue
    ASM="${LONGREAD_DIR}/${SHORT_ACC}/${SHORT_ACC}.hybrid.fasta"
    if [[ ! -s "${ASM}" ]]; then
        echo "  WARN: no assembly for ${SHORT_ACC}"
        continue
    fi

    # ---- assembly contiguity from headers (length=NNN) ----------------------
    LENS=$(grep '>' "${ASM}" | sed -n 's/.*length=\([0-9]*\).*/\1/p')
    n_contigs=$(echo "${LENS}" | grep -c .)
    assembly_size=$(echo "${LENS}" | awk '{s+=$1} END{print s+0}')

    # ---- MOB-suite: authoritative plasmid calls -----------------------------
    MOB_OUT="${LONGREAD_DIR}/${SHORT_ACC}/03_mob_recon"
    rm -rf "${MOB_OUT}"
    n_plasmids=0
    total_plasmid=0
    largest_rep="-"
    largest_size=0
    largest_mob="-"

    if mob_recon --infile "${ASM}" --outdir "${MOB_OUT}" \
            --sample_id "${SHORT_ACC}" --num_threads 4 --force 2>/dev/null; then
        MT="${MOB_OUT}/mobtyper_results.txt"
        if [[ -f "${MT}" ]]; then
            while IFS=$'\t' read -r cl sz rep mob; do
                printf "%s\t%s\t%s\t%s\t%s\n" "${SHORT_ACC}" "${cl}" "${sz}" "${rep}" "${mob}" >> "${OUT_PLASMID}"
                n_plasmids=$(( n_plasmids + 1 ))
                total_plasmid=$(( total_plasmid + sz ))
                if [[ "${sz}" -gt "${largest_size}" ]]; then
                    largest_size="${sz}"; largest_rep="${rep}"; largest_mob="${mob}"
                fi
            done < <(awk -F'\t' '
                NR==1{ for(i=1;i<=NF;i++) h[$i]=i; next }
                {
                    sid=$h["sample_id"]; cl=sid; sub(/^[^:]*:/,"",cl)
                    sz  = ("size" in h)? $h["size"]+0 : 0
                    rep = ("rep_type(s)" in h)? $h["rep_type(s)"] : "NA"
                    mob = ("predicted_mobility" in h)? $h["predicted_mobility"] : "NA"
                    print cl"\t"sz"\t"rep"\t"mob
                }' "${MT}")
        fi
    else
        echo "  WARN: mob_recon failed for ${SHORT_ACC}"
    fi

    if [[ "${n_plasmids}" -eq 0 ]]; then
        printf "%s\t(none)\t0\t-\t-\n" "${SHORT_ACC}" >> "${OUT_PLASMID}"
    fi

    # ---- auto-generated remark (data-driven; kb for readability) ------------
    if [[ "${n_plasmids}" -eq 0 ]]; then
        remark="no plasmid detected"
    else
        lkb=$(( (largest_size + 500) / 1000 ))
        if [[ "${n_plasmids}" -eq 1 ]]; then plural="plasmid"; else plural="plasmids"; fi
        remark="${n_plasmids} ${plural}; largest ~${lkb} kb ${largest_rep} (${largest_mob})"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${SHORT_ACC}" "${LONG_ACC}" "${n_contigs}" "${assembly_size}" \
        "${n_plasmids}" "${total_plasmid}" "${largest_rep}" "${remark}" >> "${OUT_ISO}"
done < "${PAIRS}"

# ---- screen summary ---------------------------------------------------------
echo ""
echo "=============================================="
echo " Hybrid assembly summary (Unicycler + MOB-suite)"
echo "=============================================="
echo ""
echo " Per-isolate overview:"
column -t -s $'\t' "${OUT_ISO}" | sed 's/^/   /'
echo ""
echo " Plasmid detail (one row per plasmid; MOB-suite):"
column -t -s $'\t' "${OUT_PLASMID}" | sed 's/^/   /'
echo ""
echo " Column meanings (overview):"
echo "   n_contigs        : assembly fragmentation (fewer = better)"
echo "   assembly_size_bp : whole assembly size (chromosome + all plasmids)"
echo "   n_plasmids       : plasmids reconstructed by MOB-suite"
echo "   total_plasmid_bp : combined bp of ALL plasmids (per-plasmid sizes in detail table)"
echo "   largest_replicon : replicon type of the largest plasmid"
echo "   remark           : plain-language interpretation (auto-generated)"
echo ""
echo " Outputs:"
echo "   overview : ${OUT_ISO}"
echo "   plasmids : ${OUT_PLASMID}"
echo "=============================================="
echo ""
echo "Hybrid-assembly summary finished on $(date)"

# End of script: 10a_longread_summary.sh

