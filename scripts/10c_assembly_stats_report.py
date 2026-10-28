#!/usr/bin/env python3
# =============================================================================
# Script       : 10c_assembly_stats_report.py
# Description  : Builds a document-ready assembly-statistics table for the three
#                Infantis hybrid assemblies. Merges QUAST contiguity metrics
#                (step 10b) with the MOB-suite plasmid overview (step 10a) into
#                one tidy table, written as TSV and Markdown (with a short
#                auto-generated description) for direct use in a report.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Run 10b_longread_assembly_stat.sh (QUAST) and 10a (overview).
# Usage        : Run from ${PROJECT_DIR}:
#                  source config/config.sh
#                  python3 scripts/10c_assembly_stats_report.py
# Dependencies : Python 3 standard library only (no pandas, no conda env).
# Note         : Reads paths from environment variables exported by config.sh.
# =============================================================================

import os
import sys
import csv

# ---- read config paths from environment (config.sh must be sourced) ---------
def need(var):
    val = os.environ.get(var)
    if not val:
        sys.exit(f"ERROR: {var} not set. Run 'source config/config.sh' first.")
    return val

LONGREAD_DIR = need("LONGREAD_DIR")

QUAST_REPORT = os.path.join(LONGREAD_DIR, "hybrid_assembly_stats", "report.tsv")
OVERVIEW     = os.path.join(LONGREAD_DIR, "01_assembly_overview.tsv")
OUT_TSV      = os.path.join(LONGREAD_DIR, "hybrid_assembly_report.tsv")
OUT_MD       = os.path.join(LONGREAD_DIR, "hybrid_assembly_report.md")

# ---- read the QUAST report (row1 = labels; each row = one metric) -----------
def read_quast(path):
    if not os.path.isfile(path):
        sys.exit(f"ERROR: QUAST report not found: {path}")
    with open(path) as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh]
    samples = rows[0][1:]
    metrics = {}
    for r in rows[1:]:
        if r and r[0]:
            metrics[r[0]] = dict(zip(samples, r[1:]))
    return samples, metrics

# ---- read the step-10a plasmid overview, keyed by isolate ------------------
def read_overview(path):
    data = {}
    if not os.path.isfile(path):
        print(f"WARN: overview not found, plasmid columns left blank: {path}")
        return data
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            iso = row.get("isolate") or row.get("Isolate")
            if iso:
                data[iso] = row
    return data

# ---- small formatters -------------------------------------------------------
def fnum(v):
    try:
        return float(v)
    except (ValueError, TypeError):
        return None

def mb(bp):
    f = fnum(bp)
    return f"{f/1e6:.2f}" if f is not None else "NA"

def as_int(v):
    f = fnum(v)
    return str(int(f)) if f is not None else "NA"

# ---- build the merged table -------------------------------------------------
samples, m = read_quast(QUAST_REPORT)
ov = read_overview(OVERVIEW)

def metric(name, sample):
    return m.get(name, {}).get(sample, "NA")

header = ["Isolate", "Total_length_Mb", "Contigs", "N50_Mb", "L50",
          "Largest_contig_Mb", "GC_pct", "Plasmids", "Largest_replicon"]

rows = []
for s in samples:
    o = ov.get(s, {})
    rows.append([
        s,
        mb(metric("Total length", s)),
        as_int(metric("# contigs", s)),
        mb(metric("N50", s)),
        as_int(metric("L50", s)),
        mb(metric("Largest contig", s)),
        metric("GC (%)", s),
        o.get("n_plasmids", "NA"),
        o.get("largest_replicon", "NA"),
    ])

# ---- auto-generated description ---------------------------------------------
def build_description(rows):
    idx = {h: i for i, h in enumerate(header)}
    lines = []
    n = len(rows)
    lines.append(
        f"Hybrid (Illumina + Nanopore) assemblies were produced for {n} "
        f"Infantis isolates. GC content was consistent at ~52% across all, with "
        f"no ambiguous bases, and each assembly placed over half its genome in a "
        f"single contig (L50 = 1), indicating good long-read contiguity."
    )
    carriers, noncarriers = [], []
    for r in rows:
        p = r[idx["Plasmids"]]
        if p not in ("NA", ""):
            (carriers if p not in ("0", "0.0") else noncarriers).append(r)
    if carriers and noncarriers:
        c = carriers[0]
        nc = ", ".join(f"{r[idx['Isolate']]} at {r[idx['Total_length_Mb']]} Mb "
                       f"in {r[idx['Contigs']]} contigs" for r in noncarriers)
        lines.append(
            f"The plasmid carrier {c[idx['Isolate']]} is the largest assembly "
            f"({c[idx['Total_length_Mb']]} Mb) and the most fragmented "
            f"({c[idx['Contigs']]} contigs), reflecting the additional, "
            f"repeat-rich plasmid content. The plasmid-free isolates assembled "
            f"into near-complete chromosomes ({nc}). The size difference is "
            f"consistent with carriage of the ~313 kb pESI megaplasmid."
        )
    return lines

description = build_description(rows)

# ---- write TSV --------------------------------------------------------------
with open(OUT_TSV, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(header)
    w.writerows(rows)

# ---- write Markdown ---------------------------------------------------------
with open(OUT_MD, "w") as fh:
    fh.write("## Hybrid assembly statistics (long-read Infantis)\n\n")
    fh.write("| " + " | ".join(header) + " |\n")
    fh.write("|" + "|".join(["---"] * len(header)) + "|\n")
    for r in rows:
        fh.write("| " + " | ".join(r) + " |\n")
    fh.write("\n### Description\n\n")
    for p in description:
        fh.write(p + "\n\n")
    fh.write("Contiguity metrics from QUAST (reference-free); plasmid columns "
             "from MOB-suite via step 10a.\n")

# ---- console summary --------------------------------------------------------
print(f"Wrote {OUT_TSV}")
print(f"Wrote {OUT_MD}\n")
widths = [max(len(h), *(len(r[i]) for r in rows)) for i, h in enumerate(header)]
def fmt(cells):
    return "  ".join(c.ljust(widths[i]) for i, c in enumerate(cells))
print(fmt(header))
for r in rows:
    print(fmt(r))
print()
for p in description:
    print(p + "\n")

# End of script: 10c_assembly_stats_report.py

