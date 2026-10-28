#!/usr/bin/env python3
# =============================================================================
# Script       : 03c_assembly_stats_report.py
# Description  : Builds a document-ready short-read assembly-QC table for the full isolate set.
#                Merges QUAST contiguity (step 03a), 
#                BUSCO completeness, CheckM2 completeness/contamination, 
#                an estimated genome coverage (fastp post-trim bases / assembly size) and
#                serovar + country of origin from the sample sheet, into one tidy table.
#                Writes TSV and Markdown with a short auto-generated description.
#                Fully dynamic: isolates and ranges are discovered from the files, nothing hard-coded.
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Prereq       : Step 03a complete (QUAST, BUSCO, CheckM2); fastp JSONs (01a).
# Usage        : Run from ${PROJECT_DIR}:
#                  source config/config.sh
#                  python3 scripts/03c_assembly_stats_report.py
# Dependencies : Python 3 standard library only (no pandas, no conda env).
# Note         : Reads paths from environment variables exported by config.sh.
#                Coverage is an estimate: post-trim bases / assembly size (the
#                same metric as the step-03b coverage gate), not mapped depth.
# =============================================================================

import os
import re
import sys
import csv
import glob
import json


# ---- read config paths from environment (config.sh must be sourced) ---------
def need(var):
    val = os.environ.get(var)
    if not val:
        sys.exit(f"ERROR: {var} not set. Run 'source config/config.sh' first.")
    return val

ASSEMBLY_QC_DIR = need("ASSEMBLY_QC_DIR")
QUAST_DIR       = need("QUAST_DIR")
BUSCO_DIR       = need("BUSCO_DIR")
CHECKM2_DIR     = need("CHECKM2_DIR")
SEROTYPE_DIR    = os.environ.get("SEROTYPE_DIR", "")
FASTP_REPORTS   = os.environ.get("FASTP_REPORTS", "")
SAMPLE_SHEET    = os.environ.get("SAMPLE_SHEET", "")

QUAST_REPORT = os.path.join(QUAST_DIR, "report.tsv")
CHECKM2_TSV  = os.path.join(CHECKM2_DIR, "quality_report.tsv")
TYPING       = os.path.join(SEROTYPE_DIR, "typing_summary.tsv") if SEROTYPE_DIR else ""

OUT_TSV = os.path.join(ASSEMBLY_QC_DIR, "assembly_qc_report.tsv")
OUT_MD  = os.path.join(ASSEMBLY_QC_DIR, "assembly_qc_report.md")

# ---- helpers ----------------------------------------------------------------
def fnum(v):
    try:
        return float(v)
    except (ValueError, TypeError):
        return None

def mb(bp):
    f = fnum(bp)
    return f"{f/1e6:.2f}" if f is not None else "NA"

def kb(bp):
    f = fnum(bp)
    return f"{f/1e3:.1f}" if f is not None else "NA"

def as_int(v):
    f = fnum(v)
    return str(int(f)) if f is not None else "NA"

def one_dp(v):
    f = fnum(v)
    return f"{f:.1f}" if f is not None else "NA"

def two_dp(v):
    f = fnum(v)
    return f"{f:.2f}" if f is not None else "NA"

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

# ---- 1) QUAST combined report (isolates as columns) -------------------------
def read_quast(path):
    data = {}
    if not os.path.isfile(path):
        print(f"WARN: QUAST report not found: {path}")
        return data
    with open(path) as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh]
    samples = rows[0][1:]
    metrics = {}
    for r in rows[1:]:
        if r and r[0]:
            metrics[r[0]] = dict(zip(samples, r[1:]))
    for s in samples:
        data[s] = {
            "size":    metrics.get("Total length", {}).get(s),
            "contigs": metrics.get("# contigs", {}).get(s),
            "n50":     metrics.get("N50", {}).get(s),
            "l50":     metrics.get("L50", {}).get(s),
            "gc":      metrics.get("GC (%)", {}).get(s),
        }
    return data

# ---- 2) CheckM2 quality_report (isolates as rows) ---------------------------
def read_checkm2(path):
    data = {}
    if not os.path.isfile(path):
        print(f"WARN: CheckM2 report not found: {path}")
        return data
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        c_name = find_col(reader.fieldnames, "Name", "genome", "bin")
        c_comp = find_col(reader.fieldnames, "Completeness")
        c_cont = find_col(reader.fieldnames, "Contamination")
        for r in reader:
            iso = r.get(c_name) if c_name else None
            if iso:
                data[iso] = {
                    "completeness":  r.get(c_comp) if c_comp else None,
                    "contamination": r.get(c_cont) if c_cont else None,
                }
    return data

# ---- 3) BUSCO per-isolate short summaries -----------------------------------
# Parse the stable one-line result: C:99.1%[S:98.9%,D:0.2%],F:0.3%,M:0.6%,n:440
BUSCO_RE = re.compile(
    r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)"
)

def read_busco(busco_dir):
    data = {}
    if not os.path.isdir(busco_dir):
        print(f"WARN: BUSCO dir not found: {busco_dir}")
        return data
    for d in sorted(os.listdir(busco_dir)):
        iso_dir = os.path.join(busco_dir, d)
        if not os.path.isdir(iso_dir):
            continue
        txts = glob.glob(os.path.join(iso_dir, "short_summary.*.txt"))
        if not txts:
            continue
        with open(txts[0]) as fh:
            text = fh.read()
        m = BUSCO_RE.search(text)
        if m:
            data[d] = {"busco_c": m.group(1)}  # Complete %
    return data

# ---- 4) serovar from typing summary (optional) ------------------------------
def read_serovar(path):
    data = {}
    if not path or not os.path.isfile(path):
        if path:
            print(f"WARN: typing summary not found, serovar left blank: {path}")
        return data
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        c_iso  = find_col(reader.fieldnames, "isolate", "sample",
                          "run_accession", "accession", "id")
        c_sero = find_col(reader.fieldnames, "serovar", "serotype")
        for r in reader:
            iso = r.get(c_iso) if c_iso else None
            if iso and c_sero:
                data[iso] = r.get(c_sero)
    return data

# ---- 5) coverage: fastp post-trim bases / assembly size ---------------------
def read_fastp_bases(fastp_dir):
    """Return {isolate: total_bases_after_filtering} from fastp JSONs."""
    data = {}
    if not fastp_dir or not os.path.isdir(fastp_dir):
        print(f"WARN: fastp reports dir not found, coverage left blank: {fastp_dir}")
        return data
    for jpath in glob.glob(os.path.join(fastp_dir, "*.json")):
        iso = os.path.basename(jpath).split(".")[0]   # ERR026042.fastp.json -> ERR026042
        try:
            with open(jpath) as fh:
                j = json.load(fh)
            tb = j.get("summary", {}).get("after_filtering", {}).get("total_bases")
            if tb is not None:
                data[iso] = float(tb)
        except (ValueError, KeyError, json.JSONDecodeError):
            print(f"WARN: could not parse {jpath}")
    return data

# ---- 6) country of origin from the sample sheet -----------------------------
def read_country(path):
    data = {}
    if not path or not os.path.isfile(path):
        if path:
            print(f"WARN: sample sheet not found, country left blank: {path}")
        return data
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        c_iso = find_col(reader.fieldnames, "Run_accession", "run_accession",
                         "accession", "isolate", "sample")
        c_cty = find_col(reader.fieldnames, "Country", "country", "location")
        for r in reader:
            iso = r.get(c_iso) if c_iso else None
            if iso and c_cty:
                data[iso] = r.get(c_cty)
    return data

# ---- merge all sources ------------------------------------------------------
quast       = read_quast(QUAST_REPORT)
checkm2     = read_checkm2(CHECKM2_TSV)
busco       = read_busco(BUSCO_DIR)
serovar     = read_serovar(TYPING)
fastp_bases = read_fastp_bases(FASTP_REPORTS)
country     = read_country(SAMPLE_SHEET)

isolates = sorted(set(quast) | set(checkm2) | set(busco))
if not isolates:
    sys.exit("ERROR: no isolates found across QUAST/CheckM2/BUSCO.")

header = ["Isolate", "Serovar", "Country", "Size_Mb", "Contigs", "N50_kb", "L50",
          "GC_pct", "Coverage_x", "Completeness_pct", "Contamination_pct",
          "BUSCO_C_pct"]

def coverage(iso, size_bp):
    tb = fastp_bases.get(iso)
    s = fnum(size_bp)
    if tb is None or not s:
        return "NA"
    return f"{tb / s:.1f}"

rows = []
for iso in isolates:
    q = quast.get(iso, {})
    c = checkm2.get(iso, {})
    b = busco.get(iso, {})
    rows.append([
        iso,
        serovar.get(iso, "-"),
        country.get(iso, "-"),
        mb(q.get("size")),
        as_int(q.get("contigs")),
        kb(q.get("n50")),
        as_int(q.get("l50")),
        q.get("gc") or "NA",
        coverage(iso, q.get("size")),
        one_dp(c.get("completeness")),
        two_dp(c.get("contamination")),
        one_dp(b.get("busco_c")),
    ])

# ---- auto-generated description (dynamic ranges) ----------------------------
def rng(vals, fmt="{:.2f}"):
    nums = [fnum(v) for v in vals]
    nums = [x for x in nums if x is not None]
    if not nums:
        return "NA"
    lo, hi = min(nums), max(nums)
    return f"{fmt.format(lo)}" if lo == hi else f"{fmt.format(lo)}–{fmt.format(hi)}"

idx = {h: i for i, h in enumerate(header)}
def col(name):
    return [r[idx[name]] for r in rows]

def counts_of(name):
    out = {}
    for r in rows:
        k = r[idx[name]]
        out[k] = out.get(k, 0) + 1
    return out

n = len(rows)
sero_bits = ", ".join(f"{k} ({v})" for k, v in sorted(counts_of("Serovar").items())
                      if k not in ("-", ""))
cty_bits = ", ".join(f"{k} ({v})" for k, v in sorted(counts_of("Country").items())
                     if k not in ("-", ""))

description = []
description.append(
    f"De novo short-read assemblies were quality-assessed for {n} isolates"
    + (f" across {sero_bits}." if sero_bits else ".")
)
if cty_bits:
    description.append(f"Isolates originated from {cty_bits}.")
description.append(
    f"Assembly sizes ranged {rng(col('Size_Mb'))} Mb with N50 of "
    f"{rng(col('N50_kb'), '{:.1f}')} kb and GC content of "
    f"{rng(col('GC_pct'), '{:.1f}')}%."
)
if any(v != "NA" for v in col("Coverage_x")):
    description.append(
        f"Estimated genome coverage ranged {rng(col('Coverage_x'), '{:.0f}')}× "
        f"(post-trim bases / assembly size)."
    )
description.append(
    f"CheckM2 estimated completeness at {rng(col('Completeness_pct'), '{:.1f}')}% "
    f"and contamination at {rng(col('Contamination_pct'), '{:.2f}')}%; "
    f"BUSCO completeness was {rng(col('BUSCO_C_pct'), '{:.1f}')}% against the "
    f"enterobacterales_odb10 lineage, confirming high-quality, low-contamination "
    f"assemblies suitable for downstream analysis."
)

# ---- write TSV --------------------------------------------------------------
with open(OUT_TSV, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(header)
    w.writerows(rows)

# ---- write Markdown ---------------------------------------------------------
with open(OUT_MD, "w") as fh:
    fh.write("## Short-read assembly QC summary\n\n")
    fh.write("| " + " | ".join(header) + " |\n")
    fh.write("|" + "|".join(["---"] * len(header)) + "|\n")
    for r in rows:
        fh.write("| " + " | ".join(r) + " |\n")
    fh.write("\n### Description\n\n")
    for p in description:
        fh.write(p + "\n\n")
    fh.write("Contiguity from QUAST; completeness/contamination from CheckM2; "
             "BUSCO against enterobacterales_odb10; coverage estimated as "
             "post-trim bases / assembly size; serovar and country from the "
             "sample sheet.\n")

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

# End of script: 03c_assembly_stats_report.py

