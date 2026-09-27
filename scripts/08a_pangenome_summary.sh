#!/bin/bash
# =============================================================================
# Script       : 08a_pangenome_summary.sh
# Description  : Summarise Panaroo pan-genome results.
#                Reports core/soft-core/shell/cloud counts, total genes, and a
#                per-serovar accessory-gene profile (mean genes per isolate by
#                serovar) from gene_presence_absence.Rtab. Also flags serovar-
#                specific genes (>=90% in one serovar, <=10% elsewhere) as a
#                light comparative view.
#                Light text parsing; run on the login node with bash.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 08_pangenome.sh and 04a typing summary first.
# Usage        : Run from ${PROJECT_DIR} :  bash scripts/08a_pangenome_summary.sh
# Note         : typing_summary.tsv col 5 = SISTR_serovar.
#                Rtab: row 1 = header (Gene <tab> isolate1 isolate2 ...),
#                cols 2..N are 0/1 presence per isolate.
# =============================================================================
set -euo pipefail
source config/config.sh

TYPING_TSV="${SEROTYPE_DIR}/typing_summary.tsv"
RTAB="${PANGENOME_DIR}/gene_presence_absence.Rtab"
PSTATS="${PANGENOME_DIR}/summary_statistics.txt"

OUT_SEROV="${PANGENOME_DIR}/01_accessory_by_serovar.tsv"
OUT_SPEC="${PANGENOME_DIR}/02_serovar_specific_genes.tsv"

[[ -f "${RTAB}" ]] || { echo "ERROR: ${RTAB} not found; run 08_pangenome.sh first"; exit 1; }

echo "Building pan-genome summary on $(date)"

# ---- 1) Headline stats from Panaroo ----------------------------------------
echo ""
echo "=============================================="
echo " Pan-genome summary (Panaroo)"
echo "=============================================="
if [[ -f "${PSTATS}" ]]; then
    sed 's/^/   /' "${PSTATS}"
fi
TOTAL_GENES=$(($(wc -l < "${RTAB}") - 1))
N_ISO=$(($(head -1 "${RTAB}" | awk '{print NF}') - 1))
echo ""
echo "   Total gene clusters : ${TOTAL_GENES}"
echo "   Isolates            : ${N_ISO}"

# ---- 2) Per-serovar mean genes per isolate ---------------------------------
printf "serovar\tn_isolates\tmean_genes_per_isolate\n" > "${OUT_SEROV}"
awk -F'\t' -v typ="${TYPING_TSV}" '
    BEGIN{
        while((getline l < typ)>0){ n=split(l,a,"\t"); if(a[1]!="isolate") sero[a[1]]=a[5] }
    }
    NR==1{
        for(i=2;i<=NF;i++){ grp[i]=(sero[$i]?sero[$i]:"NA") }
        ncol=NF
        next
    }
    {
        for(i=2;i<=ncol;i++) if($i==1) cnt[i]++
    }
    END{
        for(i=2;i<=ncol;i++){ sum[grp[i]]+=cnt[i]; ng[grp[i]]++ }
        for(g in ng) printf "%s\t%d\t%.0f\n", g, ng[g], sum[g]/ng[g]
    }
' "${RTAB}" | sort -k1,1 >> "${OUT_SEROV}"

# ---- 3) Serovar-specific genes ---------------------------------------------
# Gene present in >=90% of one serovar AND <=10% in every other serovar.
printf "serovar\tn_specific_genes\n" > "${OUT_SPEC}"
awk -F'\t' -v typ="${TYPING_TSV}" '
    BEGIN{
        while((getline l < typ)>0){ n=split(l,a,"\t"); if(a[1]!="isolate") sero[a[1]]=a[5] }
    }
    NR==1{
        for(i=2;i<=NF;i++){ grp[i]=(sero[$i]?sero[$i]:"NA"); tot[grp[i]]++ }
        for(g in tot) serovars[g]=1
        ncol=NF
        next
    }
    {
        delete pres
        for(i=2;i<=ncol;i++) if($i==1) pres[grp[i]]++
        for(sv in serovars){
            in_frac = (tot[sv]>0)? pres[sv]/tot[sv] : 0
            ok=1
            for(other in serovars){
                if(other==sv) continue
                out_frac = (tot[other]>0)? pres[other]/tot[other] : 0
                if(out_frac>0.10){ ok=0; break }
            }
            if(in_frac>=0.90 && ok) spec[sv]++
        }
    }
    END{ for(sv in serovars) printf "%s\t%d\n", sv, (spec[sv]?spec[sv]:0) }
' "${RTAB}" | sort -k1,1 >> "${OUT_SPEC}"

# ---- screen output ----------------------------------------------------------
echo ""
echo " Mean genes per isolate by serovar:"
awk -F'\t' 'NR>1{printf "   %-14s n=%-3s %s\n", $1, $2, $3}' "${OUT_SEROV}"
echo ""
echo " Serovar-specific genes (>=90% in-group, <=10% elsewhere):"
awk -F'\t' 'NR>1{printf "   %-14s %s\n", $1, $2}' "${OUT_SPEC}"
echo ""
echo " Outputs:"
echo "   accessory/serovar     : ${OUT_SEROV}"
echo "   serovar-specific genes: ${OUT_SPEC}"
echo "=============================================="
echo ""
echo "Pan-genome summary finished on $(date)"

# End of script: 08a_pangenome_summary.sh

