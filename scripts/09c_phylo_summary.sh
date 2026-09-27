#!/bin/bash
# =============================================================================
# Script       : 09c_phylo_summary.sh
# Description  : Summarise phylogeny results (core-gene tree + per-serovar SNP
#                trees). Per serovar: raw and QC-filtered pairwise SNP-distance
#                matrices (snp-dists). QC drops any isolate whose median pairwise
#                SNP distance exceeds SNP_QC_MAX (cross-lineage / failed mapping).
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 09_phylo_core.sh, 09b_phylo_snippy.sh; snp-dists installed.
# Usage        : conda activate apha_phylo; bash scripts/09c_phylo_summary.sh
# =============================================================================
set -uo pipefail
source config/config.sh

SNIPPY_BASE="${PHYLO_DIR}/02_snippy"
OUT_OVERVIEW="${PHYLO_DIR}/03_phylo_overview.tsv"
DIST_DIR="${PHYLO_DIR}/03_snp_distances"
QC_LOG="${PHYLO_DIR}/03_excluded_isolates.tsv"
mkdir -p "${DIST_DIR}"

SEROVARS=(Typhimurium Enteritidis Infantis Kentucky)
SNP_QC_MAX=10000

echo "Building phylogeny summary on $(date)"

printf "serovar\tn_isolates\tn_excluded\tn_snp_sites\tgubbins\tmin_snp\tmedian_snp\tmax_snp\n" > "${OUT_OVERVIEW}"
printf "serovar\texcluded_isolate\tmedian_snp_to_others\treason\n" > "${QC_LOG}"

for S in "${SEROVARS[@]}"; do
    WORK="${SNIPPY_BASE}/${S}"
    if [[ ! -d "${WORK}" ]]; then
        echo "  ${S}: no output, skipping"
        continue
    fi

    GUB="${WORK}/gubbins/${S}.filtered_polymorphic_sites.fasta"
    CORE_SNP="${WORK}/core.aln"
    if [[ -s "${GUB}" ]]; then
        ALN="${GUB}"; GUB_FLAG="yes"
    elif [[ -s "${CORE_SNP}" ]]; then
        ALN="${CORE_SNP}"; GUB_FLAG="fallback"
    else
        echo "  ${S}: no alignment found, skipping"
        continue
    fi

    DIST_RAW="${DIST_DIR}/${S}_snp_dist.tsv"
    if ! snp-dists -q "${ALN}" > "${DIST_RAW}" 2>/dev/null; then
        echo "  ${S}: snp-dists failed"
        continue
    fi

    # identify isolates to exclude (median dist to others > SNP_QC_MAX)
    EXC_FILE="${DIST_DIR}/${S}.exclude.txt"
    awk -F'\t' -v S="${S}" -v qc="${SNP_QC_MAX}" -v qclog="${QC_LOG}" '
        NR==1{ for(i=2;i<=NF;i++) name[i]=$i; ncol=NF; next }
        {
            row=$1
            delete d; k=0
            for(i=2;i<=ncol;i++){
                if(name[i]==row) continue
                if(name[i]=="Reference") continue
                d[++k]=$i+0
            }
            if(k==0) next
            asort(d)
            med = (k%2)? d[int(k/2)+1] : (d[k/2]+d[k/2+1])/2
            if(row!="Reference" && med>qc){
                print row
                printf "%s\t%s\t%d\tmedian SNP to others > %d (cross-lineage/failed)\n", S, row, med, qc >> qclog
            }
        }
    ' "${DIST_RAW}" > "${EXC_FILE}"

    n_excl=$(grep -c . "${EXC_FILE}" 2>/dev/null)
    n_excl=$(echo "${n_excl:-0}" | tr -d '[:space:]')

    # build QC-filtered alignment (drop Reference + excluded isolates)
    CLEAN_ALN="${DIST_DIR}/${S}.clean.aln"
    awk -v exf="${EXC_FILE}" '
        BEGIN{ while((getline l < exf)>0) if(l!="") drop[l]=1 }
        /^>/{
            name=substr($0,2); sub(/[ \\t].*/,"",name)
            keep=(name!="Reference" && !(name in drop))
        }
        keep
    ' "${ALN}" > "${CLEAN_ALN}"

    DIST_CLEAN="${DIST_DIR}/${S}_snp_dist_filtered.tsv"
    snp-dists -q "${CLEAN_ALN}" > "${DIST_CLEAN}" 2>/dev/null || true

    STATS=$(awk -F'\\t' '
        NR==1{ ncol=NF; next }
        { for(i=2;i<=ncol;i++){ if(i-1==NR-1) continue; v[++k]=$i+0 } }
        END{
            m=0; for(i=1;i<=k;i++) if(v[i]>0) w[++m]=v[i]
            if(m==0){ print "0 0 0"; exit }
            asort(w)
            printf "%d %d %d", w[1], w[int((m+1)/2)], w[m]
        }
    ' "${DIST_CLEAN}" 2>/dev/null)
    [[ -z "${STATS}" ]] && STATS="0 0 0"
    mn=$(echo "${STATS}" | awk '{print $1}')
    md=$(echo "${STATS}" | awk '{print $2}')
    mx=$(echo "${STATS}" | awk '{print $3}')

    n_snp_sites=$(awk '/^>/{next}{print length($0); exit}' "${ALN}")
    n_snp_sites=$(echo "${n_snp_sites:-0}" | tr -d '[:space:]')
    n_total=$(grep -c '^>' "${ALN}")
    n_total=$(echo "${n_total:-0}" | tr -d '[:space:]')
    has_ref=$(grep -c '^>Reference' "${ALN}")
    has_ref=$(echo "${has_ref:-0}" | tr -d '[:space:]')
    n_isolates=$(( n_total - has_ref - n_excl ))

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${S}" "${n_isolates}" "${n_excl}" "${n_snp_sites}" "${GUB_FLAG}" "${mn}" "${md}" "${mx}" >> "${OUT_OVERVIEW}"
done

echo ""
echo "=============================================="
echo " Phylogeny summary"
echo "=============================================="
echo ""
echo " Core-gene tree (all 32): ${PHYLO_DIR}/01_core_gene_tree/core_tree.treefile"
echo ""
echo " Per-serovar SNP relatedness (reference + QC-failed isolates excluded):"
column -t "${OUT_OVERVIEW}" | sed 's/^/   /'
echo ""
if [[ $(wc -l < "${QC_LOG}") -gt 1 ]]; then
    echo " Isolates excluded by relatedness QC (median SNP > ${SNP_QC_MAX}):"
    tail -n +2 "${QC_LOG}" | column -t | sed 's/^/   /'
    echo ""
fi
echo " Distance matrices:"
echo "   raw      : ${DIST_DIR}/<serovar>_snp_dist.tsv"
echo "   filtered : ${DIST_DIR}/<serovar>_snp_dist_filtered.tsv"
echo "=============================================="
echo ""
echo "Phylogeny summary finished on $(date)"

# End of script: 09c_phylo_summary.sh

