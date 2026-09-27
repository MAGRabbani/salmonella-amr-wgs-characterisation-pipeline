#!/bin/bash
# =============================================================================
# Script       : 06a_plasmid_summary.sh
# Description  : Summarise plasmid / MGE results.
#                From MOB-suite: per-plasmid table (cluster, replicon, relaxase,
#                mob type, predicted mobility, size); AMR-gene location
#                (plasmid vs chromosome) by joining AMRFinderPlus contig calls
#                to mob_recon contig_report molecule types.
#                From PlasmidFinder: independent replicon calls.
#                Cross-check: MOB-suite vs PlasmidFinder replicon concordance.
#                Light text parsing; run on the login node with bash.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 06_plasmid_mge.sh, 05_amr_virulence.sh, 04a typing first.
# Usage        : Run from ${PROJECT_DIR} :  bash scripts/06a_plasmid_summary.sh
# Note         : typing_summary.tsv col 5 = SISTR_serovar.
#                mob_recon: contig_report.txt (per contig) +
#                mobtyper_results.txt (per plasmid; sample_id = ACC:CLUSTER).
#                contig_id in contig_report carries trailing text -> trim.
# =============================================================================
set -euo pipefail
source config/config.sh

MOB_DIR="${PLASMID_DIR}/01_mob_recon"
PF_DIR="${PLASMID_DIR}/02_plasmidfinder"
TYPING_TSV="${SEROTYPE_DIR}/typing_summary.tsv"
AMRFINDER_DIR="${AMR_DIR}/01_amrfinderplus"

OUT_PLASMIDS="${PLASMID_DIR}/01_plasmids_per_isolate.tsv"
OUT_REPLICONS="${PLASMID_DIR}/02_replicons_by_serovar.tsv"
OUT_AMRLOC="${PLASMID_DIR}/03_amr_gene_location.tsv"
OUT_CONC="${PLASMID_DIR}/04_replicon_concordance.tsv"

echo "Building plasmid/MGE summary on $(date)"

# ---- serovar lookup ---------------------------------------------------------
declare -A SERO
while IFS=$'\t' read -r iso c2 c3 c4 sero rest; do
    [[ "${iso}" == "isolate" ]] && continue
    SERO["${iso}"]="${sero}"
done < "${TYPING_TSV}"

# =============================================================================
# 1) Per-plasmid table from MOB-suite mobtyper_results.txt
#    One row per plasmid. Cluster taken from sample_id (ACC:CLUSTER).
#    Columns pulled by header name (build-safe).
# =============================================================================
printf "isolate\tserovar\tplasmid_cluster\tsize\trep_types\trelaxase_types\tmobility\n" > "${OUT_PLASMIDS}"
for D in "${MOB_DIR}"/*/; do
    ACC=$(basename "${D}")
    MT="${D}mobtyper_results.txt"
    [[ -f "${MT}" ]] || continue
    s="${SERO[${ACC}]:-NA}"
    awk -F'\t' -v acc="${ACC}" -v sero="${s}" '
        NR==1{ for(i=1;i<=NF;i++) h[$i]=i; next }
        {
            sid = $h["sample_id"]; cl=sid; sub(/^[^:]*:/,"",cl)
            sz  = ("size" in h)               ? $h["size"]               : "NA"
            rep = ("rep_type(s)" in h)        ? $h["rep_type(s)"]        : "NA"
            rlx = ("relaxase_type(s)" in h)   ? $h["relaxase_type(s)"]   : "NA"
            mob = ("predicted_mobility" in h) ? $h["predicted_mobility"] : "NA"
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", acc, sero, cl, sz, rep, rlx, mob
        }
    ' "${MT}" >> "${OUT_PLASMIDS}"
done

# =============================================================================
# 2) Replicon types collapsed by serovar (distinct isolates per serovar/rep)
# =============================================================================
printf "serovar\trep_type\tn_isolates\n" > "${OUT_REPLICONS}"
awk -F'\t' 'NR>1 && $5!="NA" && $5!="-" {
    n=split($5,a,",")
    for(i=1;i<=n;i++){
        gsub(/^ +| +$/,"",a[i])
        if(!(a[i] SUBSEP $2 SUBSEP $1 in seen)){ seen[a[i] SUBSEP $2 SUBSEP $1]=1; cnt[$2 SUBSEP a[i]]++ }
    }
}
END{ for(k in cnt){ split(k,p,SUBSEP); print p[1]"\t"p[2]"\t"cnt[k] } }' "${OUT_PLASMIDS}" \
    | sort -k1,1 -k2,2 >> "${OUT_REPLICONS}"

# =============================================================================
# 3) AMR gene location: plasmid vs chromosome
#    Join AMRFinderPlus AMR hits (contig in col 3) to mob_recon contig_report
#    molecule_type. contig_report contig_id trimmed to first token.
#    AMRFinder cols: 3=Contig id, 7=Gene symbol, 10=Element type, 11=Subtype
# =============================================================================
printf "isolate\tserovar\tgene\tsubtype\tcontig\tmolecule_type\tcluster\n" > "${OUT_AMRLOC}"
for D in "${MOB_DIR}"/*/; do
    ACC=$(basename "${D}")
    CR="${D}contig_report.txt"
    AF="${AMRFINDER_DIR}/${ACC}.amrfinder.tsv"
    [[ -f "${CR}" && -f "${AF}" ]] || continue
    s="${SERO[${ACC}]:-NA}"
    awk -F'\t' -v acc="${ACC}" -v sero="${s}" '
        FNR==NR{
            if(FNR==1){ for(i=1;i<=NF;i++) hc[$i]=i; next }
            c=$hc["contig_id"]; sub(/ .*/,"",c)
            mt[c]=$hc["molecule_type"]; cl[c]=$hc["primary_cluster_id"]; next
        }
        FNR==1{ next }
        $10=="AMR"{
            g=$7; subt=$11; c=$3
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", acc, sero, g, subt, c, \
                   (c in mt?mt[c]:"NA"), (c in cl?cl[c]:"-")
        }
    ' "${CR}" "${AF}" >> "${OUT_AMRLOC}"
done

# =============================================================================
# 4) Replicon concordance: MOB-suite vs PlasmidFinder (per isolate)
# =============================================================================
printf "isolate\tserovar\tn_mob\tn_pf\tshared\tmob_only\tpf_only\n" > "${OUT_CONC}"
for D in "${MOB_DIR}"/*/; do
    ACC=$(basename "${D}")
    s="${SERO[${ACC}]:-NA}"

    mob=$(awk -F'\t' -v a="${ACC}" '$1==a && $5!="NA" && $5!="-"{
            n=split($5,x,","); for(i=1;i<=n;i++){gsub(/^ +| +$/,"",x[i]); print tolower(x[i])}
          }' "${OUT_PLASMIDS}" | sort -u)

    PF="${PF_DIR}/${ACC}/results_tab.tsv"
    pf=""
    [[ -f "${PF}" ]] && pf=$(awk -F'\t' 'NR>1{print tolower($2)}' "${PF}" | sort -u)

    nM=$(echo "${mob}" | grep -c . || true)
    nP=$(echo "${pf}"  | grep -c . || true)
    shared=$(comm -12 <(echo "${mob}") <(echo "${pf}") | grep -c . || true)
    monly=$(comm -23 <(echo "${mob}") <(echo "${pf}") | grep -c . || true)
    ponly=$(comm -13 <(echo "${mob}") <(echo "${pf}") | grep -c . || true)

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${ACC}" "${s}" "${nM}" "${nP}" "${shared}" "${monly}" "${ponly}" >> "${OUT_CONC}"
done

# =============================================================================
# Short summary to screen
# =============================================================================
N_ISO=$(ls -d "${MOB_DIR}"/*/ 2>/dev/null | wc -l)
N_PLAS=$(($(wc -l < "${OUT_PLASMIDS}") - 1))
N_AMR_PLAS=$(awk -F'\t' 'NR>1 && $6=="plasmid"' "${OUT_AMRLOC}" | wc -l)
N_AMR_CHR=$(awk -F'\t' 'NR>1 && $6=="chromosome"' "${OUT_AMRLOC}" | wc -l)

echo ""
echo "=============================================="
echo " Plasmid / MGE summary (MOB-suite primary)"
echo "=============================================="
echo " Isolates processed     : ${N_ISO}"
echo " Plasmids detected      : ${N_PLAS}"
echo " AMR hits on plasmid    : ${N_AMR_PLAS}"
echo " AMR hits on chromosome : ${N_AMR_CHR}"
echo ""
echo " Outputs:"
echo "   plasmids/isolate  : ${OUT_PLASMIDS}"
echo "   replicons/serovar : ${OUT_REPLICONS}"
echo "   AMR gene location : ${OUT_AMRLOC}"
echo "   concordance       : ${OUT_CONC}"
echo ""
echo " Replicon types by serovar:"
awk -F'\t' 'NR>1{printf "   %-14s %-20s %s\n", $1, $2, $3}' "${OUT_REPLICONS}"
echo ""
echo " AMR gene location by serovar (plasmid vs chromosome):"
awk -F'\t' 'NR>1{k=$2"\t"$6; c[k]++} END{for(x in c) print "   "x"\t"c[x]}' "${OUT_AMRLOC}" | sort
echo "=============================================="
echo ""
echo "Plasmid/MGE summary finished on $(date)"

# End of script: 06a_plasmid_summary.sh

