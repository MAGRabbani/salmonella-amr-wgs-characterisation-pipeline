#!/bin/bash
# =============================================================================
# Script       : 00_make_sample_lists.sh
# Description  : Build sample lists from samplesheet.tsv:
#                short_read_ids.txt (steps 01-09) and longread_pairs.tsv (step 10).
# Author       : Md Ataul Goni Rabbani <m.a.g.rabbani@roslin.ed.ac.uk>
# Institution  : The Roslin Institute, University of Edinburgh
# Usage        : bash scripts/00_make_sample_lists.sh
# =============================================================================
set -euo pipefail
source config/config.sh

# Short-read accessions (column 2 where Read_type contains short_read)
awk -F'\t' 'NR>1 && $13 ~ /short_read/ {print $2}' "${SAMPLE_SHEET}" \
    > "${SAMPLES_DIR}/short_read_ids.txt"

# Long-read pairs: nanopore  illumina  sample (matched by Sample column)
awk -F'\t' 'NR>1 && $13=="long_read_matched"  {long[$11]=$2}
            NR>1 && $13=="short_read_matched" {short[$11]=$2}
            END {for (s in long) print long[s]"\t"short[s]"\t"s}' \
    "${SAMPLE_SHEET}" > "${SAMPLES_DIR}/longread_pairs.tsv"

echo "Short-read: $(wc -l < "${SAMPLES_DIR}/short_read_ids.txt")  Pairs: $(wc -l < "${SAMPLES_DIR}/longread_pairs.tsv")"
# End of script: 00_make_sample_lists.sh
