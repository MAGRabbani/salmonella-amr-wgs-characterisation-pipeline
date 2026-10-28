#!/usr/bin/env python3
# =============================================================================
# Script       : 13_report_summary.py
# Description  : Consolidates the pipeline's per-step summary tables into one
#                document-ready report (Markdown + console): dataset overview
#                (typing), acquired AMR and mutations, multidrug resistance,
#                resistance by antimicrobial class, a critically-important ESBL
#                callout, and resistance-determinant location (chromosome vs
#                plasmid). Missing inputs are skipped with a WARN. Fully dynamic:
#                counts, classes and ranges are derived from the files.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Steps 04a, 05a and 06a complete (summary tables present).
# Usage        : Run from ${PROJECT_DIR}:
#                  source config/config.sh
#                  python3 scripts/13_report_summary.py
# Dependencies : Python 3 standard library only (no pandas, no conda env).
# Note         : Reads paths from environment variables exported by config.sh.
#                Writes to ${RESULTS_DIR}/13_summary_report/ (created if absent).
# =============================================================================

import os
import re
import sys
import csv
from collections import defaultdict, Counter

# ---- read config paths from environment (config.sh must be sourced) ---------
def need(var):
    val = os.environ.get(var)
    if not val:
        sys.exit(f"ERROR: {var} not set. Run 'source config/config.sh' first.")
    return val

RESULTS_DIR  = need("RESULTS_DIR")
SEROTYPE_DIR = need("SEROTYPE_DIR")
AMR_DIR      = need("AMR_DIR")
PLASMID_DIR  = need("PLASMID_DIR")

TYPING    = os.path.join(SEROTYPE_DIR, "typing_summary.tsv")
AMR_MAT   = os.path.join(AMR_DIR, "03_amr_matrix_annotated.tsv")
PHENO_MDR = os.path.join(AMR_DIR, "04_amr_phenotype_mdr.tsv")
GENE_LOC  = os.path.join(PLASMID_DIR, "03_amr_gene_location.tsv")

OUT_DIR = os.path.join(RESULTS_DIR, "13_summary_report")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_MD  = os.path.join(OUT_DIR, "summary_report.md")

# ---- antimicrobial class from gene name (pattern-based, order matters) -------
def amr_class(gene):
    g = gene.lower()
    if g.startswith("gyra") or g.startswith("parc") or g.startswith("gyrb") \
       or g.startswith("pare") or g.startswith("qnr"):
        return "Fluoroquinolones"
    if g.startswith("bla"):
        return "Beta-lactams"
    if g.startswith(("aac", "aad", "aph", "arm", "rmt", "ant", "str")):
        return "Aminoglycosides"
    if g.startswith("sul"):
        return "Sulphonamides"
    if g.startswith("tet"):
        return "Tetracyclines"
    if g.startswith("dfr"):
        return "Trimethoprim"
    if g.startswith(("flo", "cat", "cml", "cmx")):
        return "Phenicols"
    if g.startswith("sat"):
        return "Streptothricin"
    if g.startswith("mcr"):
        return "Colistin"
    if g.startswith(("erm", "mef", "mph")):
        return "Macrolides"
    return "Other"

def is_esbl(gene):
    # CTX-M enzymes are the ESBLs in this dataset
    return "ctx-m" in gene.lower()

# ---- helpers ----------------------------------------------------------------
def read_tsv(path):
    if not os.path.isfile(path):
        print(f"WARN: not found, section skipped: {path}")
        return None, None
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return reader.fieldnames, list(reader)

def find_col(fieldnames, *cands):
    if not fieldnames:
        return None
    low = {f.lower(): f for f in fieldnames}
    for c in cands:
        if c.lower() in low:
            return low[c.lower()]
    for f in fieldnames:
        for c in cands:
            if c.lower() in f.lower():
                return f
    return None

def is_true(v):
    return str(v).strip().lower() in ("1", "1.0", "true", "yes", "y", "mdr")

report = []
def line(s=""):
    report.append(s)

line("# Salmonella AMR WGS pipeline — summary report\n")

# ---- 1) dataset overview (from typing_summary) ------------------------------
total_isolates = None
fn, rows = read_tsv(TYPING)
if rows:
    c_iso  = find_col(fn, "isolate", "sample", "run_accession", "accession", "id")
    c_sero = find_col(fn, "serovar", "serotype")
    c_st   = find_col(fn, "st", "sequence_type", "mlst")
    isolates = set()
    serovar_counts = Counter()
    st_by_serovar = defaultdict(set)
    for r in rows:
        iso  = r.get(c_iso) if c_iso else None
        sero = (r.get(c_sero) if c_sero else "Unknown") or "Unknown"
        if iso:
            isolates.add(iso)
            serovar_counts[sero] += 1
            if c_st and r.get(c_st):
                st_by_serovar[sero].add(r[c_st])
    total_isolates = len(isolates)
    line("## Dataset overview\n")
    line(f"- Total isolates analysed: **{total_isolates}**\n")
    line("| Serovar | n | Sequence type(s) |")
    line("|---|---|---|")
    for sero in sorted(serovar_counts):
        sts = ", ".join(sorted(st_by_serovar.get(sero, []))) or "-"
        line(f"| {sero} | {serovar_counts[sero]} | {sts} |")
    line("")

# ---- 2) acquired AMR and mutations (from annotated matrix) ------------------
amr_rows = None
amr_cols = None
c_iso = c_sero = None
fn, rows = read_tsv(AMR_MAT)
if rows:
    amr_rows = rows
    c_iso  = find_col(fn, "isolate")
    c_sero = find_col(fn, "serovar")
    c_st   = find_col(fn, "st", "sequence_type")
    meta = {c_iso, c_sero, c_st}
    amr_cols = [f for f in fn if f not in meta]
    n_with_amr = len(rows)
    line("## Acquired AMR and resistance mutations\n")
    if total_isolates:
        clean = total_isolates - n_with_amr
        line(f"- Isolates with at least one acquired gene or resistance "
             f"mutation: **{n_with_amr}** of {total_isolates} "
             f"({clean} with none).\n")
    else:
        line(f"- Isolates with at least one acquired gene or resistance "
             f"mutation: **{n_with_amr}**.\n")
    line("| Isolate | Serovar | Genes/mutations detected |")
    line("|---|---|---|")
    for r in rows:
        n = sum(1 for g in amr_cols if is_true(r.get(g, "0")))
        line(f"| {r.get(c_iso)} | {r.get(c_sero)} | {n} |")
    line("")

# ---- 3) multidrug resistance (from phenotype/MDR table) --------------------
fn, rows = read_tsv(PHENO_MDR)
if rows:
    c_mdr = find_col(fn, "mdr", "is_mdr", "mdr_flag", "multidrug")
    if c_mdr:
        n_mdr = sum(1 for r in rows if is_true(r.get(c_mdr, "")))
        line(f"- Multidrug-resistant (resistance to ≥3 classes): "
             f"**{n_mdr}** isolates.\n")
    else:
        print("WARN: MDR column not found in phenotype table; MDR line skipped.")

# ---- 4) resistance by antimicrobial class (from the matrix) -----------------
if amr_rows and amr_cols:
    # per serovar: set of classes with >=1 positive isolate, and isolate counts
    classes_by_serovar = defaultdict(lambda: Counter())  # serovar -> class -> n iso
    for r in amr_rows:
        sero = r.get(c_sero, "Unknown")
        hit_classes = set()
        for g in amr_cols:
            if is_true(r.get(g, "0")):
                hit_classes.add(amr_class(g))
        for cl in hit_classes:
            classes_by_serovar[sero][cl] += 1
    line("## Resistance by antimicrobial class\n")
    line("| Serovar | Classes detected (isolates affected) |")
    line("|---|---|")
    for sero in sorted(classes_by_serovar):
        cc = classes_by_serovar[sero]
        bits = ", ".join(f"{cl} ({cc[cl]})" for cl in sorted(cc))
        line(f"| {sero} | {bits} |")
    line("")

# ---- 5) critically-important ESBL callout -----------------------------------
# genes from the matrix; location/cluster from the gene-location file if present
if amr_rows and amr_cols:
    esbl_cols = [g for g in amr_cols if is_esbl(g)]
    # build location lookup: (isolate, gene) -> cluster (from gene_loc)
    loc_lookup = {}
    fnL, rowsL = read_tsv(GENE_LOC)
    if rowsL:
        cL_iso = find_col(fnL, "isolate")
        cL_gene = find_col(fnL, "gene")
        cL_mol = find_col(fnL, "molecule_type", "molecule", "location")
        cL_clu = find_col(fnL, "cluster")
        for r in rowsL:
            key = (r.get(cL_iso), r.get(cL_gene))
            mol = (r.get(cL_mol) or "").strip()
            clu = r.get(cL_clu) if cL_clu else "-"
            if mol.lower() == "plasmid" and clu and clu != "-":
                loc_lookup[key] = f"plasmid {clu}"
            elif mol:
                loc_lookup[key] = mol
    if esbl_cols:
        line("## Critically important resistance (ESBLs)\n")
        line("Extended-spectrum beta-lactamases (CTX-M) confer resistance to "
             "third-generation cephalosporins, a critically important class in "
             "human medicine.\n")
        line("| Isolate | Serovar | ESBL gene | Location |")
        line("|---|---|---|---|")
        any_esbl = False
        for r in amr_rows:
            iso = r.get(c_iso)
            sero = r.get(c_sero)
            for g in esbl_cols:
                if is_true(r.get(g, "0")):
                    any_esbl = True
                    loc = loc_lookup.get((iso, g), "-")
                    line(f"| {iso} | {sero} | {g} | {loc} |")
        if not any_esbl:
            line("| - | - | none detected | - |")
        line("")

# ---- 6) resistance determinant location (chromosome vs plasmid) ------------
fn, rows = read_tsv(GENE_LOC)
if rows:
    c_sL   = find_col(fn, "serovar")
    c_mol  = find_col(fn, "molecule_type", "molecule", "location")
    c_clu  = find_col(fn, "cluster")
    chrom = Counter()
    plas  = Counter()
    clusters = defaultdict(set)
    for r in rows:
        sero = (r.get(c_sL) if c_sL else "Unknown") or "Unknown"
        mol  = (r.get(c_mol) or "").strip().lower() if c_mol else ""
        if mol == "chromosome":
            chrom[sero] += 1
        elif mol == "plasmid":
            plas[sero] += 1
            cl = r.get(c_clu) if c_clu else None
            if cl and cl != "-":
                clusters[sero].add(cl)
    line("## Resistance determinant location (chromosome vs plasmid)\n")
    line("| Serovar | Chromosomal | Plasmid-borne | Plasmid clusters |")
    line("|---|---|---|---|")
    for sero in sorted(set(list(chrom) + list(plas))):
        cl = ", ".join(sorted(clusters.get(sero, []))) or "-"
        line(f"| {sero} | {chrom.get(sero, 0)} | {plas.get(sero, 0)} | {cl} |")
    line("")

# ---- write report -----------------------------------------------------------
with open(OUT_MD, "w") as fh:
    fh.write("\n".join(report) + "\n")

print(f"Wrote {OUT_MD}\n")
print("\n".join(report))

# End of script: 13_report_summary.py

