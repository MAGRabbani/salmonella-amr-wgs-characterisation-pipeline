#!/bin/bash
# =============================================================================
# Script       : 05a_amr_summary.sh
# Description  : Summarise AMR and virulence results.
#                From AMRFinderPlus: gene long-table (acquired vs point), a
#                presence/absence matrix (serovar/ST annotated), per-isolate
#                drug-class phenotype and MDR flag.
#                Cross-check: AMRFinderPlus vs ResFinder (acquired genes).
#                Virulence: light VFDB gene count per isolate.
#                Light text parsing; run on the login node with bash.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 05_amr_virulence.sh and 04a typing summary first.
# Usage        : Run from ${PROJECT_DIR} :  bash scripts/05a_amr_summary.sh
# Note         : typing_summary.tsv columns: 1 isolate, 2 MLST_ST, 3 cgMLST_ST,
#                4 subspecies, 5 SISTR_serovar (serovar read from column 5).
# =============================================================================
set -euo pipefail
source config/config.sh

AMRFINDER_DIR="${AMR_DIR}/01_amrfinderplus"
ABRICATE_DIR="${AMR_DIR}/02_abricate"
TYPING_TSV="${SEROTYPE_DIR}/typing_summary.tsv"

OUT_LONG="${AMR_DIR}/01_amr_genes_long.tsv"
OUT_MATRIX="${AMR_DIR}/02_amr_gene_matrix.tsv"
OUT_ANNOT="${AMR_DIR}/03_amr_matrix_annotated.tsv"
OUT_PHENO="${AMR_DIR}/04_amr_phenotype_mdr.tsv"
OUT_CONC="${AMR_DIR}/05_amr_concordance.tsv"
OUT_VF="${AMR_DIR}/06_virulence_counts.tsv"

echo "Building AMR/virulence summary on $(date)"

# =============================================================================
# 1) AMRFinderPlus long table (AMR class only), tag acquired vs point
#    Cols: 1=Name 7=Gene 10=ElementType 11=ElementSubtype 12=Class
# =============================================================================
printf "isolate\tgene\tclass\tsubtype\n" > "${OUT_LONG}"
for F in "${AMRFINDER_DIR}"/*.amrfinder.tsv; do
    awk -F'\t' 'NR>1 && $10=="AMR" {print $1"\t"$7"\t"$12"\t"$11}' "${F}" >> "${OUT_LONG}"
done

# =============================================================================
# 2) Presence/absence matrix (isolate x gene)
# =============================================================================
awk -F'\t' '
    NR>1 { iso[$1]=1; gene[$2]=1; seen[$1"\t"$2]=1 }
    END {
        printf "isolate"
        n=0; for (g in gene) glist[++n]=g
        for (i=1;i<=n;i++) for (j=i+1;j<=n;j++) if (glist[j]<glist[i]){t=glist[i];glist[i]=glist[j];glist[j]=t}
        for (i=1;i<=n;i++) printf "\t%s", glist[i]
        printf "\n"
        m=0; for (s in iso) ilist[++m]=s
        for (i=1;i<=m;i++) for (j=i+1;j<=m;j++) if (ilist[j]<ilist[i]){t=ilist[i];ilist[i]=ilist[j];ilist[j]=t}
        for (a=1;a<=m;a++){
            printf "%s", ilist[a]
            for (i=1;i<=n;i++) printf "\t%s", (seen[ilist[a]"\t"glist[i]]?1:0)
            printf "\n"
        }
    }' "${OUT_LONG}" > "${OUT_MATRIX}"

# =============================================================================
# 3) Annotated matrix: prepend serovar + ST, sort by serovar
#    serovar from typing col 5, ST from typing col 2
# =============================================================================
awk -F'\t' -v typ="${TYPING_TSV}" '
    BEGIN{ while((getline l < typ)>0){ n=split(l,a,"\t"); if(a[1]!="isolate"){st[a[1]]=a[2]; sero[a[1]]=a[5]} } }
    NR==1 { print "serovar\tST\t"$0; next }
    { s=(sero[$1]?sero[$1]:"NA"); t=(st[$1]?st[$1]:"NA"); print s"\t"t"\t"$0 }
' "${OUT_MATRIX}" | (read -r h; echo "$h"; sort -k1,1 -k2,2) > "${OUT_ANNOT}"

# =============================================================================
# 4) Per-isolate drug-class phenotype + MDR flag (>=3 classes)
#    serovar from typing col 5, ST from typing col 2
# =============================================================================
printf "isolate\tserovar\tST\tn_classes\tclasses\tMDR\n" > "${OUT_PHENO}"
awk -F'\t' -v typ="${TYPING_TSV}" '
    BEGIN{ while((getline l < typ)>0){ n=split(l,a,"\t"); if(a[1]!="isolate"){st[a[1]]=a[2]; sero[a[1]]=a[5]} } }
    NR>1 { cls[$1"\t"$3]=1; iso[$1]=1 }
    END{
        for (i in iso) ilist[++m]=i
        for (x=1;x<=m;x++) for (y=x+1;y<=m;y++) if (ilist[y]<ilist[x]){t=ilist[x];ilist[x]=ilist[y];ilist[y]=t}
        for (x=1;x<=m;x++){
            s=ilist[x]; delete seen; nc=0; list=""
            for (k in cls){ split(k,p,"\t"); if(p[1]==s && !(p[2] in seen)){seen[p[2]]=1; nc++; list=(list==""?p[2]:list";"p[2])} }
            mdr=(nc>=3?"MDR":"-")
            printf "%s\t%s\t%s\t%s\t%s\t%s\n", s, (sero[s]?sero[s]:"NA"), (st[s]?st[s]:"NA"), nc, list, mdr
        }
    }' "${OUT_LONG}" | sort -k2,2 -k3,3 >> "${OUT_PHENO}"


# =============================================================================
# 5) Concordance: AMRFinderPlus vs ResFinder (acquired genes only)
#    AMRFinder acquired = subtype AMR (drop POINT). ResFinder col6=GENE.
# =============================================================================
printf "isolate\tn_amrfinder\tn_resfinder\tshared\tamrfinder_only\tresfinder_only\n" > "${OUT_CONC}"
for ASM in "${ASSEMBLY_DIR}"/*.fasta; do
    ACC=$(basename "${ASM}" .fasta)
    AF="${AMRFINDER_DIR}/${ACC}.amrfinder.tsv"
    RF="${ABRICATE_DIR}/${ACC}.resfinder.tab"
    [[ -f "${AF}" && -f "${RF}" ]] || continue

    afg=$(awk -F'\t' 'NR>1 && $10=="AMR" && $11=="AMR"{print tolower($7)}' "${AF}" \
          | sed 's/[-_].*//' | sort -u)
    rfg=$(awk -F'\t' 'NR>1{print tolower($6)}' "${RF}" \
          | sed 's/[-_].*//' | sort -u)

    nA=$(echo "${afg}" | grep -c . || true)
    nR=$(echo "${rfg}" | grep -c . || true)
    shared=$(comm -12 <(echo "${afg}") <(echo "${rfg}") | grep -c . || true)
    aonly=$(comm -23 <(echo "${afg}") <(echo "${rfg}") | grep -c . || true)
    ronly=$(comm -13 <(echo "${afg}") <(echo "${rfg}") | grep -c . || true)

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${ACC}" "${nA}" "${nR}" "${shared}" "${aonly}" "${ronly}" >> "${OUT_CONC}"
done

# =============================================================================
# 6) Virulence: VFDB gene count per isolate (light summary)
#    ABRicate VFDB: col6=GENE.
# =============================================================================
printf "isolate\tn_virulence_genes\n" > "${OUT_VF}"
for ASM in "${ASSEMBLY_DIR}"/*.fasta; do
    ACC=$(basename "${ASM}" .fasta)
    VF="${ABRICATE_DIR}/${ACC}.vfdb.tab"
    [[ -f "${VF}" ]] || continue
    n=$(awk -F'\t' 'NR>1' "${VF}" | wc -l)
    printf "%s\t%s\n" "${ACC}" "${n}" >> "${OUT_VF}"
done

# =============================================================================
# Short summary to screen
# =============================================================================
N_ISO=$(($(wc -l < "${OUT_MATRIX}") - 1))
N_GENE=$(($(head -1 "${OUT_MATRIX}" | awk -F'\t' '{print NF-1}')))
N_MDR=$(awk -F'\t' 'NR>1 && $6=="MDR"' "${OUT_PHENO}" | wc -l)

echo ""
echo "=============================================="
echo " AMR / virulence summary (AMRFinderPlus primary)"
echo "=============================================="
echo " Isolates            : ${N_ISO}"
echo " Distinct AMR genes  : ${N_GENE}"
echo " MDR isolates (>=3)  : ${N_MDR}/${N_ISO}"
echo ""
echo " Outputs:"
echo "   long table        : ${OUT_LONG}"
echo "   gene matrix       : ${OUT_MATRIX}"
echo "   annotated matrix  : ${OUT_ANNOT}"
echo "   phenotype + MDR   : ${OUT_PHENO}"
echo "   concordance       : ${OUT_CONC}"
echo "   virulence counts  : ${OUT_VF}"
echo ""
echo " Gene frequency across panel (AMR class):"
awk -F'\t' 'NR>1{print $2}' "${OUT_LONG}" | sort | uniq -c | sort -rn | sed 's/^/   /'
echo ""
echo " Distinct AMR genes per serovar:"
awk -F'\t' -v typ="${TYPING_TSV}" '
    BEGIN{ while((getline l < typ)>0){ n=split(l,a,"\t"); if(a[1]!="isolate") sero[a[1]]=a[5] } }
    NR>1{ key=sero[$1]"\t"$2; if(!(key in s)){s[key]=1; c[sero[$1]]++} }
    END{ for(x in c) printf "   %-14s %s\n", x, c[x] }
' "${OUT_LONG}"
echo "=============================================="
echo ""
echo "AMR/virulence summary finished on $(date)"

# End of script: 05a_amr_summary.sh

