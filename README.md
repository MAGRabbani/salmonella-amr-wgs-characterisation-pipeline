# Salmonella AMR & Mobile Genetic Elements — WGS Surveillance Pipeline

**Where does the resistance sit — on the chromosome, or on a mobile element that
can pass to the next isolate?** That answer decides whether an AMR problem stays
put or spreads through a flock. This pipeline is built to resolve it.

An end-to-end, reproducible workflow for whole-genome characterisation of
*Salmonella enterica* from poultry — the routine a surveillance lab runs on every
isolate: **serotype and sequence type**, **acquired AMR and virulence genes**,
their **genetic context and link to mobile genetic elements**, and **comparative,
pan-genome and phylogenetic** placement. It handles **short- and long-read** data,
runs in UNIX from bash and R, and is config-driven to scale from a few isolates to
a large collection. Every call is cross-checked against a second tool, so the
output is defensible, not one tool's best guess.

Proof on a four-serovar set: the chromosomal SGI1 block in Typhimurium,
fluoroquinolone resistance in Kentucky, and — confirmed five independent ways —
the ESBL-carrying **pESI megaplasmid** in Infantis. The approach is informed by
APHA's WGS surveillance of *Salmonella* and AMR in livestock, and built on the
author's peer-reviewed WGS work
([Rabbani *et al.*, 2024, *Scientific Data*](https://www.nature.com/articles/s41597-024-04291-z)),
hands-on bioinformatics training, and published tool documentation and tutorials.

> 🧬 *Salmonella* Enteritidis · Typhimurium · Infantis · Kentucky — serotyping,
> MLST, AMR & virulence, plasmids/MGEs, pan-genome, phylogeny, hybrid assembly,
> all from one config.



## Pipeline overview

Whole-genome sequencing pipeline for four poultry-associated *Salmonella*
serovars (Enteritidis, Typhimurium, Infantis, Kentucky), following a public-
health surveillance logic — identify the isolate, find its resistance
determinants, work out **where** those determinants sit (mobile or not), and
place each isolate against its serovar background.

**Flow:** build sample lists → download → read QC → trim → assembly →
assembly QC → serotyping/MLST → AMR/virulence → plasmid/MGE → annotation →
pan-genome → phylogeny → long-read hybrid → figures → final MultiQC.

- All scripts source a single central config: `config/config.sh`
  (edit paths there; nothing else needs changing).
- Numbered scripts (`00`–`12`) are the pipeline steps; `a`/`b`/`c` suffixes are
  summary or comparison helpers for the preceding step.
- Most steps run on the SGE cluster via `qsub`; several are array jobs
  (one task per isolate). Run all scripts from the project root.
- Four serovars, 32 short-read isolates + 3 matched Nanopore runs.

## What each step does (at a glance)

| Step | Script                                | Does                                                      | Key output                                    |
| ---- | ------------------------------------- | --------------------------------------------------------- | --------------------------------------------- |
| 00   | `00_make_sample_lists.sh`             | Build sample sheet + read-set lists (run first)           | `samples/*.tsv`, `samples/*.txt`              |
| 00a  | `00a_download_raw_reads.sh`           | Fetch raw FASTQ from ENA (Illumina + ONT)                 | `results/00_raw_reads/`                       |
| 00b  | `00b_download_ref_genomes.sh`         | Download the main reference genome (LT2)                  | `ref/reference.fna`                           |
| 00c  | `00c_setup_databases.sh`              | Stage AMR/plasmid/annotation databases                    | `ref/` databases                              |
| 01   | `01_fastqc.sh`                        | Raw-read quality control (FastQC), array job              | `results/01_qc/01_fastqc/`                    |
| 01a  | `01a_trim.sh`                         | Adapter/quality trimming (fastp), array job               | `results/01_qc/01a_trim/`                     |
| 01b  | `01b_fastqc_trim.sh`                  | FastQC on trimmed reads, array job                        | `results/01_qc/01b_fastqc_trim/`              |
| 02   | `02_multiqc.sh`                       | Aggregate raw-read QC                                     | `results/02_multiqc/`                         |
| 02a  | `02a_multiqc_trim.sh`                 | Aggregate trimmed-read QC                                 | `results/02_multiqc/`                         |
| 03   | `03_assembly.sh`                      | De novo assembly (Shovill/SPAdes), array job              | `results/03_assembly/<iso>/`                  |
| 03a  | `03a_assembly_qc.sh`                  | QUAST + BUSCO + CheckM2, array job                        | `results/03a_assembly_qc/`                    |
| 03b  | `03b_assembly_qc_summary.sh`          | Assembly QC summary + coverage gate                       | `results/03a_assembly_qc/`                    |
| 04   | `04_serotype_mlst.sh`                 | MLST → SISTR → SeqSero2, array job                        | `results/04_serotype_mlst/`                   |
| 04a  | `04a_serotype_mlst_typing_summary.sh` | Serovar/ST concordance table                              | `results/04_serotype_mlst/typing_summary.tsv` |
| 05   | `05_amr_virulence.sh`                 | AMRFinderPlus + ABRicate, array job                       | `results/05_amr_virulence/`                   |
| 05a  | `05a_amr_summary.sh`                  | Gene matrix, phenotype/MDR, virulence counts              | `results/05_amr_virulence/`                   |
| 06   | `06_plasmid_mge.sh`                   | MOB-suite + PlasmidFinder, array job                      | `results/06_plasmid_mge/`                     |
| 06a  | `06a_plasmid_summary.sh`              | AMR gene location (chromosome vs plasmid)                 | `results/06_plasmid_mge/`                     |
| 07   | `07_annotation.sh`                    | Genome annotation (Bakta), array job                      | `results/07_annotation/<iso>/`                |
| 07a  | `07a_annotation_summary.sh`           | Per-isolate CDS/tRNA/rRNA summary                         | `results/07_annotation/`                      |
| 08   | `08_pangenome.sh`                     | Pan-genome (Panaroo)                                      | `results/08_pangenome/`                       |
| 08a  | `08a_pangenome_summary.sh`            | Core/accessory + serovar-specific genes                   | `results/08_pangenome/`                       |
| 09   | `09_phylo_core.sh`                    | Whole-set core-gene tree (IQ-TREE)                        | `results/09_phylogeny/01_core_gene_tree/`     |
| 09a  | `09a_download_phylo_refs.sh`          | Download per-serovar phylo refs (Ent/Inf/Kent) + link LT2 | `ref/phylo_refs/`                             |
| 09b  | `09b_phylo_snippy.sh`                 | Per-serovar SNP trees (Snippy→Gubbins→IQ-TREE), array     | `results/09_phylogeny/02_snippy/`             |
| 09c  | `09c_phylo_summary.sh`                | SNP-distance matrices + outlier exclusion                 | `results/09_phylogeny/03_snp_distances/`      |
| 10   | `10_longread_hybrid.sh`               | Hybrid assembly of 3 Infantis pairs, array job            | `results/10_longread_hybrid/`                 |
| 10a  | `10a_longread_summary.sh`             | Assembly + plasmid overview                               | `results/10_longread_hybrid/`                 |
| 11   | `11_figures.sh` (+ `11a`–`11e`.R)     | Presentation figures (R: ggtree, pheatmap)                | `results/11_figures/`                         |
| 12   | `12_multiqc_final.sh`                 | Final aggregated QC across the run                        | `results/12_multiqc_final/`                   |

> Directory paths above match the variables defined in `config/config.sh`.
> Logs are written to `logs/`; scratch to `tmp/`.



## Repository layout

```
apha_salmonella_amr/
├── config/config.sh     # central paths & variables (edit before running)
├── scripts/             # numbered pipeline scripts + R figure scripts
├── envs/                # conda env files (one .yml per environment)
├── samples/             # accession list, sample sheet, long-read pairs
├── results/             # per-step outputs (00_raw_reads ... 12_multiqc_final)
├── ref/                 # references, databases, phylo refs (not in repo)
├── logs/                # job logs
└── README.md
```



## Requirements

- **conda** (Miniconda or Anaconda) — the pipeline builds a per-stage environment
  from each `envs/*.yml`. Set `CONDA_BASE` in `config.sh` to your conda install.
- **Reference genome and databases** — the LT2 reference, per-serovar phylogeny
  references, Bakta DB, CheckM2 DB, BUSCO lineage and PlasmidFinder DB. These are
  **large** and are **not** stored in this repository; the setup steps
  (00b, 00c) download them, and their paths are defined in `config.sh`.
- **Internet access** on the machine running the setup and download steps
  (00a–00c fetch reads, references and databases from ENA/NCBI).
- **A job scheduler** for the per-isolate steps — these are written as **SGE /
  `qsub`** array jobs. On a different scheduler (SLURM, PBS) or a single machine,
  the analysis commands inside each script are unchanged; only the job-submission
  header and the array indexing need adapting. The summary and figure steps are
  plain `bash` and run anywhere.

> **Memory note:** the per-isolate job headers request memory per task. Two tools
> were tuned to stay within a bounded allocation — NanoPlot was dropped (its
> numpy/OpenBLAS backend reserved very large amounts of virtual memory) and
> Shovill's Java `--ram` is held below the total requested. Adjust the resource
> requests in the script headers to suit your system.



## Configuration

All scripts source **`config/config.sh`** — edit it once, nothing else needs
changing. Key variables:

| Variable                           | Meaning                                   |
| ---------------------------------- | ----------------------------------------- |
| `PROJECT_DIR`                      | Project root                              |
| `SAMPLES_DIR` / `SAMPLE_SHEET`     | Samples directory and master sample sheet |
| `CONDA_BASE`                       | Personal conda base (holds conda.sh)      |
| `REF_DIR` / `BAKTA_DB`             | Reference + annotation database paths     |
| `CHECKM2_DB` / `BUSCO_LINEAGE`     | Assembly-QC database paths                |
| `PLASMIDFINDER_DB`                 | PlasmidFinder database                    |
| `REF_TYPHIMURIUM` … `REF_KENTUCKY` | Per-serovar phylogeny references          |
| `ENV_*`                            | Conda environment names (one per stage)   |
| `RESULTS_DIR` (+ per-step dirs)    | Output locations for each step            |
| `THREADS`                          | Compute default (8; overridden per step)  |



## Environment strategy

Each stage has its own conda environment — SPAdes, CheckM2, Bakta and the R
stack do not share a Python happily, so they are kept apart. The `.yml` files in
`envs/` recreate them.

Create every environment once, before running the pipeline. Initialise your conda
first (personal install as used here, or the Eddie anaconda module), then create
each env:

```bash
# initialise conda (personal install, as used here)
source "${CONDA_BASE}/etc/profile.d/conda.sh"

conda env create -f envs/apha_wgs.yml
conda env create -f envs/checkm2.yml
conda env create -f envs/apha_typing.yml
conda env create -f envs/apha_amrPlasfinder.yml
conda env create -f envs/apha_annotation.yml
conda env create -f envs/apha_pangenome.yml
conda env create -f envs/apha_phylo.yml
conda env create -f envs/apha_longread.yml
conda env create -f envs/apha_figures.yml
```

| Env                  | Used for                                             |
| -------------------- | ---------------------------------------------------- |
| `apha_wgs`           | FastQC, fastp, MultiQC, Shovill/SPAdes, QUAST, BUSCO |
| `checkm2`            | CheckM2 (kept separate — Python conflict)            |
| `apha_typing`        | mlst, SISTR, SeqSero2                                |
| `apha_amrPlasfinder` | AMRFinderPlus, ABRicate, PlasmidFinder, MOB-suite    |
| `apha_annotation`    | Bakta, Prokka                                        |
| `apha_pangenome`     | Panaroo, MAFFT                                       |
| `apha_phylo`         | Snippy, Gubbins, IQ-TREE, snp-dists                  |
| `apha_longread`      | Filtlong, Unicycler, Flye                            |
| `apha_figures`       | R: ggtree, tidyverse, pheatmap, patchwork            |

> **Conda + module quirk:** in batch scripts, `source activate <env>` works
> after `module load igmm/apps/anaconda/2023.03`, whereas `conda activate` does
> not (no `conda init` under the module). Login-node summary scripts source a
> personal conda and use `conda activate`. `config.sh` itself does **not**
> initialise conda.
>
> **Login-node steps need conda initialised.** The summary and figure steps call
> `conda activate <env>` directly, so conda must be initialised in your shell.
> Run `conda init` once, or `source "${CONDA_BASE}/etc/profile.d/conda.sh"` at
> the start of your session. Set `CONDA_BASE` in `config.sh` to your conda
> install; the setup script (00c) sources it for you.



## Dataset

36 sequencing runs in total: **32 short-read isolates** that assembled and
passed QC, plus **3 matched Nanopore (GridION)** runs for the Infantis long-read
pairs.

| Serovar     | n    | Dominant ST       | Notable                                  |
| ----------- | ---- | ----------------- | ---------------------------------------- |
| Enteritidis | 10   | ST11 (+ 1 ST616)  | clean AMR background                     |
| Infantis    | 9    | ST32              | pESI+ cluster + 2 pESI-negative isolates |
| Typhimurium | 7    | ST19              | SGI1 chromosomal MDR                     |
| Kentucky    | 6    | ST198 (+ 1 ST314) | FQ-resistant; ST314 is a distant outlier |

One isolate (**ERR026016**) was dropped at the coverage gate — 17× coverage and
a 771 kb assembly, well below the 30× threshold. QC reports are kept for all 33
for audit; analysis runs on the 32 that passed (enforced via
`samples/short_read_ids.txt`).

## How samples are chosen

The pipeline starts from a curated list of ENA/SRA run accessions. Samples were
selected as publicly available Illumina paired-end runs of *Salmonella enterica*
from **poultry**, across the four target serovars (Enteritidis, Typhimurium,
Infantis, Kentucky), plus three Infantis isolates with matched Nanopore
(GridION) runs for the hybrid-assembly branch.

The run accessions are listed in `samples/accessions.txt` (one per line),
**provided in this repository**. Step-00 reads this file to build the sample
sheet and read-set lists. To run on your own data, edit `accessions.txt` and
re-run Step-00 — nothing else needs changing.



## How samples are chosen

Isolate metadata was extracted from **Enterobase** — publicly available
*Salmonella enterica* genomes from **poultry** across the four target serovars
(Enteritidis, Typhimurium, Infantis, Kentucky), spanning several countries,
years and submitting labs (including APHA and its VLA Weybridge predecessor).
Three South Korean Infantis isolates were chosen with matched Nanopore (GridION)
runs to demonstrate the hybrid-assembly branch and pESI reconstruction.

The full metadata is committed as `samples/samplesheet.tsv` — one row per
sequencing run, carrying accession, serovar, source, country, year, submitting
lab, BioProject, platform, read type and known AMR profile. Step-00 reads this
sheet to build the read-set lists the pipeline runs on. To run on your own data,
replace `samplesheet.tsv` and re-run Step-00.

---



## Step-00: Build sample lists

Reads the master sample sheet and writes the two read-set lists the rest of the
pipeline iterates over. This is the **first** script to run.

### Prerequisite

- No conda env needed — core shell utilities only.
- Sample sheet present at `${SAMPLE_SHEET}` (`samples/samplesheet.tsv`).

### What it does

- Extracts every short-read accession (`Read_type` matching `short_read`) into
  `short_read_ids.txt` — this covers plain `short_read` and the matched
  `short_read_matched` Infantis isolates.
- Pairs the matched long/short Infantis runs by their `Sample` (BioSample) into
  `longread_pairs.tsv` (nanopore accession, illumina accession, sample).

### Sample sheet columns

`samples/samplesheet.tsv` is tab-separated, one sequencing run per row. Columns
used by the pipeline are shown in bold; the rest are metadata for provenance.

| #      | Column              | Notes                                                    |
| ------ | ------------------- | -------------------------------------------------------- |
| 1      | Sl_No               | serial number                                            |
| **2**  | **Run_accession**   | ENA/SRA accession — the download/lookup key              |
| 3      | Enterobase_name     | Enterobase strain name                                   |
| 4      | Barcode             | original lab barcode                                     |
| **5**  | **Serovar**         | expected serovar (cross-checked against Step-04)         |
| 6      | Source              | isolation source (poultry)                               |
| 7      | Country             | country of origin                                        |
| 8      | Year                | year of isolation                                        |
| 9      | Submitting_lab      | submitting laboratory                                    |
| 10     | BioProject          | ENA/SRA BioProject                                       |
| **11** | **Sample**          | BioSample — used to pair matched long/short runs         |
| 12     | Instrument_platform | sequencing platform                                      |
| **13** | **Read_type**       | `short_read`, `short_read_matched`, `long_read_matched`  |
| 14     | AMR_profile         | expected AMR (from Enterobase; cross-checked in Step-05) |
| 15     | Notes               | free-text (coverage, resistance context, flags)          |

The expected serovar and AMR profile here are only cross-checks — the pipeline's
serovar calls come from SISTR/SeqSero2 (Step-04a) and its AMR calls from
AMRFinderPlus (Step-05), so a mismatch against the sheet is itself informative.

### Outputs (`samples/`)

- `short_read_ids.txt` — short-read accessions, one per line (33 runs; iterated
  by Steps 01–09).
- `longread_pairs.tsv` — the 3 matched Infantis pairs
  (nanopore, illumina, sample; used by Step-10).

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script (core shell, no env needed)
bash scripts/00_make_sample_lists.sh
```

### Notes

- Array jobs pull the isolate for task `N` with
  `sed -n "${SGE_TASK_ID}p" short_read_ids.txt`, so **`-t 1-N` must match the
  line count** of the list the step reads.
- `short_read_ids.txt` starts at 33; ERR026016 is removed at the Step-03b
  coverage gate, leaving the 32-isolate analysis set for Steps 03 onward. QC
  (Steps 01–02) still runs on all 33 for audit.

---



## Step-00a: Download raw reads

Downloads the FASTQ for every run in the sample sheet from ENA, resolving each
accession to its FASTQ URLs via the ENA API. A **login-node bash script**
(`wget -c` resumes partial downloads and skips completed ones, so it is safe to
re-run). Run after the sample sheet is in place.

### Prerequisite

- No conda env needed — uses system `wget` and `awk`.
- Sample sheet present at `${SAMPLE_SHEET}` (`samples/samplesheet.tsv`).

### What it does

- Reads every accession from column 2 of the sample sheet (all 36 runs —
  33 Illumina + 3 ONT).
- Queries the ENA `filereport` API for each accession's `fastq_ftp` URLs.
- Downloads each FASTQ into `${RAW_DIR}` with `wget -c` (resume/skip), retrying
  on timeout.
- Logs progress to `logs/00_download.log` and reports the final file count.

### Outputs (`results/00_raw_reads/`)

- Paired Illumina: `<accession>_1.fastq.gz`, `<accession>_2.fastq.gz`.
- ONT: `<accession>.fastq.gz` (single-end).

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script (login-node; system wget, no env needed)
bash scripts/00a_download_raw_reads.sh
```

### Notes

- Safe to re-run — `wget -c` resumes incomplete files and skips completed ones.
- Both Illumina and ONT reads are fetched here; the ONT runs feed only the
  Step-10 hybrid branch, the Illumina reads feed everything else.
- If an accession returns no FASTQ URL it is logged and skipped, not fatal.

---



## Step-00b: Download reference genome

Downloads the *S.* Typhimurium LT2 reference (GCA_000006945.2) from ENA into
`${REF_DIR}`. This is the pipeline's main reference genome; the per-serovar
phylogeny references are fetched separately at Step-09a. A **login-node bash
script** — run once.

### Ref genomes details

| Serovar     | Accession       | Strain / note                     |
| ----------- | --------------- | --------------------------------- |
| Typhimurium | GCA_000006945.2 | LT2                               |
| Enteritidis | GCA_000009505.1 | P125109                           |
| Infantis    | GCA_037776255.1 | Z1323CSL0027, pESI+ (4 replicons) |
| Kentucky    | GCA_002952975.1 | PU131, complete                   |

### Prerequisite

- No conda env needed — uses system `wget`.
- `config.sh` paths `REF_DIR` and `REF_GENOME` set.

### What it does

- Downloads the LT2 assembly FASTA (gzipped) via the ENA browser API.
- Unzips it to `${REF_GENOME}` (`ref/reference.fna`) and reports the sequence
  count.

### Outputs (`ref/`)

- `reference.fna` — *S.* Typhimurium LT2 (GCA_000006945.2).

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script (login-node; system wget, no env needed)
bash scripts/00b_download_ref_genomes.sh
```

### Notes

- This provides the single LT2 reference only. The Enteritidis, Infantis and
  Kentucky phylogeny references are downloaded at Step-09a, which also copies this
  LT2 into `phylo_refs/Typhimurium/`.

---

## Step-00c: Set up reference databases

One-time download of the databases the pipeline needs — the BUSCO lineage,
CheckM2 DB, AMRFinderPlus DB, ABRicate DBs and the Bakta light DB. A
**login-node bash script** that **switches conda env per tool internally**, so
you run it directly without activating anything first. Run once, interactively
(needs internet).

### Prerequisite

- The relevant conda envs already created: `apha_wgs`, `checkm2`,
  `apha_amrPlasfinder`, `apha_annotation` (the script activates each in turn).
- Personal conda initialised (the script sources it directly — see note).
- Sufficient disk under `${REF_DIR}` (Bakta light ~1.3 GB, CheckM2 DB ~3 GB).

### What it does

- **BUSCO** (env `apha_wgs`) — downloads `enterobacterales_odb10` into
  `${REF_DIR}/busco_downloads/`.
- **CheckM2** (env `checkm2`) — downloads the DIAMOND DB to `${CHECKM2_DB}`.
- **AMRFinderPlus** (env `apha_amrPlasfinder`) — `amrfinder -u` updates the DB.
- **ABRicate** (same env) — `--setupdb` builds the bundled DBs (card, resfinder,
  vfdb, plasmidfinder).
- **Bakta** (env `apha_annotation`) — downloads the light DB to
  `${REF_DIR}/bakta_db`.

### Outputs (`ref/`)

- `busco_downloads/lineages/enterobacterales_odb10/`
- `checkm2_db/`
- `bakta_db/db-light/`
- AMRFinderPlus and ABRicate DBs (installed within their env prefixes).

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script (login-node; activates each env internally, no pre-activation needed)
bash scripts/00c_setup_databases.sh
```

### Notes

- **Personal conda, not the module:** the script sources the personal conda
  (`.../anaconda/etc/profile.d/conda.sh`) directly rather than
  `module load anaconda`, because the module overrides the env's Python and
  breaks import-based tools like BUSCO. Edit that path if your conda lives
  elsewhere.
- **PlasmidFinder DB** comes bundled via ABRicate's `--setupdb` here; the
  standalone PlasmidFinder DB used in Step-06 (`${PLASMIDFINDER_DB}`) ships with
  the `apha_amrPlasfinder` env at `$CONDA_PREFIX/db/plasmidfinder/src`.
- Bakta light DB is used by default; for the full DB (~30 GB) change
  `--type light` to `--type full`.

---



## Step-01: Read QC (FastQC)

Quality-control of the raw paired-end Illumina reads. Runs as an **SGE array
job** — one task per isolate.

### Prerequisite

- Conda env (apha_wgs):
- Step-00a complete (raw reads present); Step-00 sample list built.

### What it does

- Task `N` processes the isolate on line `N` (`SGE_TASK_ID` → `sed`).
- Runs FastQC on both read files per isolate.
- Writes reports to `results/01_qc/01_fastqc/`.

### Outputs (`results/01_qc/01_fastqc/`)

- `<accession>_{1,2}_fastqc.{html,zip}` per isolate.

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
qsub scripts/01_fastqc.sh
```

### Notes

- **Conda Env: Activating conda env (apha_wgs) need once upto step-?**
- **Array size must match isolate count** — adjust `-t 1-N`. Test `-t 1-2` first.
- `-tc 10` throttles concurrent tasks (courteous on the shared cluster).
- Output is consumed by Step-02 (MultiQC).

---

## Step-01a: Adapter/quality trimming (fastp)

Trims adapters and low-quality bases from the raw reads. Runs as an **SGE array
job** (one task per isolate). fastp writes a JSON (for MultiQC) and an HTML
report per isolate.

### Prerequisite

- Conda env (apha_wgs):
- Step-00a complete (raw reads present).

### What it does

- Locates paired reads (`*_1.fastq.gz` / `*_2.fastq.gz`) for each isolate.
- Runs fastp with adapter detection + quality trimming.
- Writes trimmed reads to `${TRIM_DIR}` and reports to `fastp_reports/`.

### Outputs (`results/01_qc/01a_trim/`)

- `<acc>_{1,2}.trimmed.fastq.gz` — trimmed reads (feed every downstream step).
- `fastp_reports/<acc>.{json,html}` — per-isolate reports.

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
qsub scripts/01a_trim.sh
```

### Notes

- **Array size must match isolate count** — adjust `-t 1-N`.
- Trimmed reads are the canonical read set from here on — assembly (Step-03) and
  the hybrid branch (Step-10) both read from `${TRIM_DIR}`.

---

## Step-01b: QC of trimmed reads (FastQC)

Quality-controls the fastp-trimmed reads to confirm trimming did what it should.
Runs as an **SGE array job** (one task per isolate).

### Prerequisite

- Conda env (apha_wgs):
- Step-01a complete (trimmed reads present).

### What it does

- Runs FastQC on the trimmed read pair per isolate.
- Writes reports to `results/01_qc/01b_fastqc_trim/`.

### Outputs (`results/01_qc/01b_fastqc_trim/`)

- `<accession>_{1,2}.trimmed_fastqc.{html,zip}` per isolate.

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
qsub scripts/01b_fastqc_trim.sh
```

### Notes

- **Array size must match isolate count** — adjust `-t 1-N`.
- Output is aggregated alongside the raw FastQC in Step-02a.

---

## Step-02: Aggregate raw-read QC (MultiQC)

Combines the raw-read FastQC reports into a single interactive MultiQC report. A
**single job** — run after Step-01 finishes for all isolates.

### Prerequisite

- Conda env (apha_wgs):
- Step-01 complete.

### What it does

- Scans `results/01_qc/01_fastqc/` for FastQC outputs.
- Produces one consolidated report at `results/02_multiqc/`.

### Outputs (`results/02_multiqc/`)

- `multiqc_report.html` — aggregated raw-read QC.

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
qsub scripts/02_multiqc.sh
```

---

## Step-02a: Aggregate trimmed-read QC (MultiQC)

Combines the trimmed-read FastQC reports plus the fastp JSONs into a single
report, documenting the effect of trimming. A **single job** — run after
Step-01b finishes for all isolates.

### Prerequisite

- Conda env (apha_wgs):
- Step-01a and Step-01b complete.

### What it does

- Scans `results/01_qc/01b_fastqc_trim/` and the fastp JSONs.
- Produces one consolidated report at `results/02_multiqc/`.

### Outputs (`results/02_multiqc/`)

- Trimmed-read MultiQC report (raw vs trimmed comparison).

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
qsub scripts/02a_multiqc_trim.sh
```

---



## Step-03: De novo assembly (Shovill)

Assembles each isolate's trimmed reads into a draft genome. Runs as an **SGE
array job** (one task per isolate).

### Prerequisite

- Conda env (apha_wgs):
- Step-01a complete (trimmed reads present).

### What it does

- Runs Shovill (SPAdes under the hood) on each trimmed read pair.
- Parameters: `--minlen 200 --mincov 5`, estimated genome size `4.8M`.
- Writes one assembly directory per isolate.

### Outputs (`results/03_assembly/<iso>/`)

- `contigs.fa` — draft assembly.
- Shovill logs and intermediate SPAdes output.

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
qsub scripts/03_assembly.sh
```

### Notes

- **Array size must match isolate count** — adjust `-t 1-N`.
- **Memory:** Shovill's `--ram` must be set **well below** the total requested,
  or the SPAdes Java step hits an out-of-memory error. On Eddie, remember total
  memory is `h_vmem × sharedmem` — set `--ram` under that.

---

## Step-03a: Assembly QC (QUAST + BUSCO + CheckM2)

Runs three complementary assembly-quality tools on every draft: QUAST for
contiguity, BUSCO for gene-completeness, CheckM2 for completeness and
contamination. Runs as an **SGE array job**.

### Prerequisite

- Conda envs (apha_wgs, checkm2):
- Databases staged (Step-00c): BUSCO lineage (enterobacterales_odb10),
  CheckM2 DB.
- Step-03 assemblies present.

### What it does

- **QUAST** — contig count, N50, total length, GC.
- **BUSCO** — single-copy ortholog completeness (enterobacterales_odb10).
- **CheckM2** — completeness and contamination estimates.

### Outputs (`results/03a_assembly_qc/`)

- `01_quast/` — QUAST report (single combined `report.tsv` across all isolates).
- `02_busco/` — per-isolate BUSCO short summaries.
- `03_checkm2/quality_report.tsv` — CheckM2 completeness/contamination table.

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
qsub scripts/03a_assembly_qc.sh
```

### Notes

- **Module vs conda Python clash:** the anaconda module overrides Python and
  breaks import-based tools (BUSCO). This step uses `conda activate` *after* the
  module load for QUAST/BUSCO so the correct Python is picked up.
- CheckM2 runs in its own env because its Python pin conflicts with `apha_wgs`.

---

## Step-03b: Assembly QC summary & coverage gate

Consolidates the three QC tools into one table and applies the coverage/quality
gate that defines the analysis set. A **login-node summary** script.

### Prerequisite

- Conda env (apha_wgs):
- Step-03a complete.

### What it does

- Merges QUAST, BUSCO and CheckM2 metrics per isolate into one summary table.
- Applies the **coverage gate** (≥30×) and quality thresholds; flags failures.

### Outputs (`results/03a_assembly_qc/`)

- Combined QC summary table (per-isolate metrics + pass/fail).

### Result on this dataset

All 32 retained isolates **PASS**: CheckM2 100% complete, 0.04–0.64%
contamination; BUSCO 99.1–99.5%; size 4.63–5.03 Mb; N50 92.5–547 kb;
GC 51.9–52.3%; zero ambiguous bases.

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
bash scripts/03b_assembly_qc_summary.sh
```

### Notes

- **Coverage gate:** ERR026016 was excluded here — 17× coverage, 771 kb
  assembly, below the 30× threshold. QC reports are retained for all 33 for
  audit; the analysis set is the 32 that pass, and `short_read_ids.txt` is where
  the 33→32 reduction is enforced for every step from 03 onward.
- Summary scripts use `set -uo pipefail` (not `-e`) to avoid mid-loop aborts.

---



## Step-04: Serotyping & MLST

Identifies each isolate to serovar and sequence type using three orthogonal
tools, so calls are cross-checked rather than trusted from one method. Runs as
an **SGE array job**.

### Prerequisite

- Conda env (apha_typing):
- Step-03 assemblies present.

### What it does

- **mlst** — 7-gene MLST sequence type.
- **SISTR** — serovar prediction (cgMLST + antigen).
- **SeqSero2** — independent serovar prediction from the antigenic formula.

### Outputs (`results/04_serotype_mlst/`)

- Per-tool results per isolate (mlst, SISTR, SeqSero2).

### Run

```bash
# Create conda env
module load igmm/apps/anaconda/2023.03 		# or use a personal conda
conda env create -f envs/apha_typing.yml
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
qsub scripts/04_serotype_mlst.sh
```

### Notes

- **Array size must match isolate count** — adjust `-t 1-N`.
- Running SISTR and SeqSero2 side by side gives an independent concordance check
  — important in surveillance, where a mis-called serovar changes the
  public-health interpretation.

---

## Step-04a: Typing summary & concordance

Builds the master typing table used by every downstream step for serovar
grouping, and checks SISTR–SeqSero2 agreement. A **login-node summary** script.

### Prerequisite

- Conda env (apha_typing):
- Step-04 complete.

### What it does

- Collates MLST, SISTR and SeqSero2 calls per isolate.
- Reports serovar/ST concordance.

### Outputs (`results/04_serotype_mlst/`)

- `typing_summary.tsv` — the master table (serovar in column 5); consumed by
  Steps 05a, 08, 09 and 11.

### Result on this dataset

**32/32 SISTR–SeqSero2 agree.** Antigenic profiles: Typhimurium 4:i:1,2;
Enteritidis 9:g,m:-; Kentucky 8:i:z6; Infantis 7:r:1,5. Sequence types:
Enteritidis ST11 (+1 ST616), Typhimurium ST19, Infantis ST32, Kentucky ST198
(+1 ST314 outlier).

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script (interactively)
bash scripts/04a_serotype_mlst_typing_summary.sh
```

### Notes

- `typing_summary.tsv` is the single source of serovar labels downstream — the
  figures and phylogeny group by this computed call, not by any expected-serovar
  column in the sample sheet.

---



## Step-05: AMR & virulence detection

Detects acquired AMR genes, resistance point mutations and virulence factors
across all 32 isolates — the core of the project. Runs as an **SGE array job**.

### Prerequisite

- Conda env (apha_amrPlasfinder):
- AMRFinderPlus DB present (Step-00c); Step-03 assemblies.

### What it does

- **AMRFinderPlus** (`--organism Salmonella`) — primary caller, with
  Salmonella-specific point-mutation detection.
- **ABRicate** against **card**, **resfinder** and **vfdb** — cross-check for
  acquired genes and virulence factors.

### Outputs (`results/05_amr_virulence/`)

- Per-isolate AMRFinderPlus and ABRicate results (per database).

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
qsub scripts/05_amr_virulence.sh
```

### Notes

- **Array size must match isolate count** — adjust `-t 1-N`.
- AMRFinderPlus with the organism flag catches the point mutations (gyrA, parC)
  that a plain gene-presence tool misses — the fluoroquinolone resistance in
  Kentucky and Infantis is point-mutation-driven, so this matters.
- Intrinsic *aac(6')-Iaa* appears in the resfinder track only — AMRFinderPlus
  correctly filters it as intrinsic, so it is **not** counted as acquired
  resistance.

---

## Step-05a: AMR summary, phenotype & MDR

Turns the raw AMR calls into the tables and matrices that drive the figures and
the interpretation. A **login-node summary** script.

### Prerequisite

- Conda env (apha_wgs):
- Step-05 complete; `typing_summary.tsv` (Step-04a).

### What it does

- Builds a long table, a gene presence/absence matrix, and a serovar-annotated
  matrix.
- Infers phenotype and flags MDR (resistance to ≥3 antimicrobial classes).
- Reports AMRFinderPlus vs ABRicate concordance and virulence-gene counts.

### Outputs (`results/05_amr_virulence/`)

- `01_amr_genes_long.tsv` — one row per gene call.
- `02_amr_gene_matrix.tsv` — isolate × gene presence/absence.
- `03_amr_matrix_annotated.tsv` — matrix with serovar + ST (drives fig3).
- `04_amr_phenotype_mdr.tsv` — inferred phenotype and MDR flag.
- `05_amr_concordance.tsv` — AMRFinderPlus vs ABRicate agreement.
- `06_virulence_counts.tsv` — virulence-gene counts per isolate.

### Result on this dataset

- **Typhimurium** — SGI1 penta-block (*aadA2, blaCARB-2, floR, sul1, tet(G)*),
  all chromosomal; *gyrA* point mutations.
- **Infantis** — *bla*CTX-M-1/65 (plasmid-borne); *gyrA*_D87G (chromosomal).
- **Kentucky** — *gyrA*_S83F + *gyrA*_D87 + *parC*_S80I; *bla*CMY-2,
  *bla*CTX-M-65, *bla*TEM-1.
- **Enteritidis** — no acquired AMR.

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
bash scripts/05a_amr_summary.sh
```

### Notes

- The chromosome-vs-plasmid location of each gene is not decided here — it comes
  from Step-06a, which cross-references these calls against the MOB-suite plasmid
  assignments. Step-05a records *what* resistance is present; Step-06 records
  *where* it sits.

---

## Step-06: Plasmid & mobile-element typing

Determines **where** the resistance sits — chromosome or plasmid — the question
that matters most for how resistance spreads. Runs as an **SGE array job**.

### Prerequisite

- Conda env (apha_amrPlasfinder):
- PlasmidFinder DB (`${PLASMIDFINDER_DB}`, Step-00c); Step-03 assemblies.

### What it does

- **MOB-suite** — reconstructs and types plasmids (replicon, mobility, cluster).
- **PlasmidFinder** — independent replicon detection.
- Cross-references AMR gene contigs against plasmid vs chromosome calls.

### Outputs (`results/06_plasmid_mge/`)

- MOB-suite `mobtyper_results.txt` and contig reports per isolate.
- PlasmidFinder results per isolate.

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
qsub scripts/06_plasmid_mge.sh
```

### Notes

- **Array size must match isolate count** — adjust `-t 1-N`.
- MOB-suite runs on plain defaults with `--sample_id` labelling and the full DB
  (no taxonomy filter).
- **Parsing gotchas:** MOB-suite `contig_report` contig IDs carry trailing text
  (trim to the first token); mobtyper `sample_id` is written as `ACC:CLUSTER`.

---

## Step-06a: Plasmid summary & AMR gene location

Produces the chromosome-vs-plasmid location call for every AMR gene — the input
to the location strip in the AMR figure, and the payoff of the whole "what vs
where" logic set up in Step-05. A **login-node summary** script.

### Prerequisite

- Conda env (apha_wgs):
- Steps 05 and 06 complete.

### What it does

- Assigns each AMR gene to chromosome or plasmid, per isolate.
- Summarises replicons by serovar and checks replicon concordance
  (MOB-suite vs PlasmidFinder).

### Outputs (`results/06_plasmid_mge/`)

- `01_plasmids_per_isolate.tsv` — plasmids and replicon types per isolate.
- `02_replicons_by_serovar.tsv` — replicon distribution by serovar.
- `03_amr_gene_location.tsv` — chromosome/plasmid call per AMR gene (drives the
  fig3 location strip).
- `04_replicon_concordance.tsv` — MOB-suite vs PlasmidFinder agreement.

### Result on this dataset

- **Typhimurium** — 35 AMR calls chromosomal / 1 plasmid; the SGI1 block is
  chromosomal. pSLT (AB460, IncFIB+IncFII) present in all Typhimurium.
- **Infantis** — 7 chromosomal / 34 plasmid; the pESI cluster (AC358, IncFIB,
  ~273 kb) in 7/9 Infantis (absent from the two short-read pESI-negative
  isolates, SRR30149746 and SRR30149758).
- **Kentucky** — 29 chromosomal / 7 plasmid.
- One *aadA2*-on-pSLT in ERR232523 — an integron-mobility case (resistance
  normally chromosomal in Typhimurium, here caught on a plasmid).

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
bash scripts/06a_plasmid_summary.sh
```

### Notes

- Location calls come from MOB-suite contig assignment, not contig-size guesses —
  a gene on a plasmid-assigned contig is called plasmid-borne regardless of
  contig length.
- This is where the surveillance-relevant distinction is made: chromosomal SGI1
  is stably inherited, whereas the pESI/pSLT plasmid-borne resistance is
  horizontally mobile.

---



## Step-07: Genome annotation (Bakta)

Annotates each assembly, providing the GFF3 input for the pan-genome and the
CDS/feature counts for the annotation summary. Runs as an **SGE array job**.

### Prerequisite

- Conda env (apha_annotation):
- Bakta DB (`${BAKTA_DB}`, db-light, schema v5.1) staged (Step-00c).
- Step-03 assemblies present.

### What it does

- Runs Bakta on each assembly, producing GFF3, GBFF, FAA and feature tables.

### Outputs (`results/07_annotation/<iso>/`)

- Bakta output per isolate (GFF3 used downstream by the pan-genome).

### Run

```bash
# Load conda env
module load igmm/apps/anaconda/2023.03  		# or use a personal conda
conda env create -f envs/apha_annotation.yml
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
qsub scripts/07_annotation.sh
```

### Notes

- **Array size must match isolate count** — adjust `-t 1-N`.
- **pyhmmer pin:** the `apha_annotation` env pins `pyhmmer=0.10.15` —
  version 0.12.3 breaks Bakta with `'str' object has no attribute 'decode'`. The
  provided `envs/apha_annotation.yml` already carries the pin.
- **Bakta GFF3 for Panaroo:** Bakta writes single-`#` comment lines that break
  Panaroo's GFF parser — these are stripped in Step-08.
- The db-light schema includes the amrfinderplus DB, so annotations carry AMR
  context.

---

## Step-07a: Annotation summary

Summarises annotation features per isolate and flags any with an unusually low
CDS count. A **login-node summary** script.

### Prerequisite

- Conda env (apha_wgs):
- Step-07 complete.

### What it does

- Tabulates size, GC, coding density, CDS, tRNA, rRNA, ncRNA and hypothetical
  counts per isolate.

### Outputs (`results/07_annotation/`)

- `01_annotation_summary.tsv` — per-isolate feature counts (feeds the pESI
  evidence panel in fig4).

### Result on this dataset

CDS counts by serovar: Typhimurium ~4715, Kentucky ~4565, Infantis ~4534,
Enteritidis ~4405. All 32 annotated cleanly. Within Infantis, the pESI-carrying
SRR30149748 sits at ~4662 CDS versus ~4265–4342 for the two pESI-negative
isolates — a difference that reappears as one line of the five-way pESI evidence
(Step-10).

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
bash scripts/07a_annotation_summary.sh
```

### Notes

- **CDS flag threshold** relaxed from 4400 to 4200 — Enteritidis and the
  reduced-pESI Infantis isolates are naturally leaner, so the original threshold
  false-flagged them.

---



## Step-08: Pan-genome (Panaroo)

Builds the pan-genome across all 32 isolates — core vs accessory gene structure,
and the core-gene alignment that feeds the whole-set tree. A **single job**.

### Prerequisite

- Conda env (apha_pangenome):
- Step-07 annotations (GFF3) present.

### What it does

- Cleans Bakta GFF3 (strips single-`#` comment lines) into
  `${ANNOT_DIR}/gff_clean/`.
- Runs Panaroo with `--remove-invalid-genes` (drops frame-broken genes from
  fragmented contigs).
- Produces the core-gene alignment (MAFFT).

### Outputs (`results/08_pangenome/`)

- `gene_presence_absence.Rtab` — gene presence/absence matrix (drives fig5).
- Core-gene alignment (feeds Step-09 whole-set tree).
- Panaroo summary and graph files.

### Run

```bash
# Load conda env
module load igmm/apps/anaconda/2023.03  		# or use a personal conda
conda env create -f envs/apha_pangenome.yml
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
qsub scripts/08_pangenome.sh
```

### Result on this dataset

Core 3845, soft-core 38, shell 1809, cloud 571 — total 6263 gene clusters.
Core alignment 116 MB, 32 sequences.

### Notes

- Both fixes are required: raw Bakta GFF3 otherwise breaks Panaroo's parser, and
  frame-broken genes from fragmented contigs otherwise corrupt the alignment.
- The cleaned GFF3 files in `gff_clean/` are the actual Panaroo input, not the
  raw Bakta output.

---

## Step-08a: Pan-genome summary

Summarises accessory structure and identifies serovar-specific genes. A
**login-node summary** script.

### Prerequisite

- Conda env (apha_wgs):
- Step-08 complete.

### What it does

- Summarises core/soft-core/shell/cloud structure.
- Counts genes specific to each serovar.

### Outputs (`results/08_pangenome/`)

- `01_accessory_by_serovar.tsv` — mean genes per isolate, by serovar.
- `02_serovar_specific_genes.tsv` — serovar-specific gene counts (drives the
  fig5 right panel).

### Result on this dataset

Serovar-specific genes: Typhimurium 322, Enteritidis 183, Kentucky 99,
Infantis 71.

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
bash scripts/08a_pangenome_summary.sh
```

### Notes

- The Infantis figure (71) is **artificially low** — the pESI genes fall below
  the 90% presence threshold because they split across the pESI+ and
  pESI-negative Infantis isolates. This is a threshold artefact, not a biological
  floor, and is worth stating explicitly whenever the number is shown.

---



## Step-09: Whole-set core-gene phylogeny

Builds the genome-wide tree across all 32 isolates from the Panaroo core-gene
alignment. A **single job**.

### Prerequisite

- Conda env (apha_phylo):
- Step-08 core-gene alignment present.

### What it does

- Runs IQ-TREE (`-m MFP -B 1000`) on the Panaroo core-gene alignment.

### Outputs (`results/09_phylogeny/01_core_gene_tree/`)

- `core_tree.treefile` — the whole-set tree (rendered as fig1 in Step-11).

### Result on this dataset

74,400 parsimony-informative sites. Four clean serovar clades; Typhimurium
basal and near-clonal, Infantis and Kentucky as sister groups.

### Run

```bash
# Load conda env
module load igmm/apps/anaconda/2023.03   		# or use a personal conda
conda env create -f envs/apha_phylo.yml
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
qsub scripts/09_phylo_core.sh
```

### Notes

- This tree shows serovar separation across the whole set; within-serovar
  resolution comes from the reference-based SNP trees (Step-09b).

---

## Step-09a: Download per-serovar phylogeny references

Fetches the four ENA reference genomes used to anchor the reference-based SNP
trees, and prepares them for Snippy. A **single job** — run once, before
Step-09b.

### Prerequisite

- Conda env (apha_wgs):
- `config.sh` paths for the four `REF_<SEROVAR>` accessions.

### What it does

- Downloads each per-serovar reference into `ref/phylo_refs/<serovar>/` (same
  references as Step-00b), ready for the Snippy mapping step.

### Outputs (`ref/phylo_refs/<serovar>/`)

- One reference FASTA per serovar (see the reference table in Step-00b).

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
qsub scripts/09a_download_phylo_refs.sh
```

### Notes

- Downloads Enteritidis (P125109), Infantis (Z1323CSL0027, pESI+) and Kentucky
  (PU131). For Typhimurium it reuses the LT2 reference from Step-00b (copied into
  `phylo_refs/Typhimurium/`), fetching it fresh only if `${REF_GENOME}` is
  absent.
- Skips any reference already present (`-s` check), so it is safe to re-run.

---

## Step-09b: Per-serovar SNP phylogeny (Snippy → Gubbins → IQ-TREE)

Builds a high-resolution, reference-based SNP tree **within** each serovar,
where the core-gene tree lacks resolution. Runs as an **SGE array job**
(`-t 1-4`, one serovar per task).

### Prerequisite

- Conda env (apha_phylo):
- Step-09a references; trimmed reads (Step-01a); typing table (Step-04a).

### What it does

Per serovar:

1. **Snippy** — map each isolate's reads to the serovar reference; call variants.
2. **snippy-core** + **snippy-clean_full_aln** — build the core SNP alignment
   (includes `Reference` as an implicit outgroup).
3. **Gubbins** — remove recombinant regions.
4. **IQ-TREE** — build the recombination-filtered SNP tree.

### Outputs (`results/09_phylogeny/02_snippy/<Serovar>/`)

- `iqtree/<Serovar>_snp.treefile` — the per-serovar tree (rendered as fig2).
- Snippy, snippy-core and Gubbins intermediates.

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
qsub scripts/09b_phylo_snippy.sh
```

### Notes

- **Snippy version:** upgraded 4.0.2 → 4.6.0 — the old version misread
  samtools 1.20 as `<1.7` and refused to run. The `envs/apha_phylo.yml` carries
  4.6.0.
- Gubbins ran successfully (not the fallback) on all four serovars.
- The array is `-t 1-4` (one serovar per task), not one-per-isolate — the four
  serovars come from the typing table.

---

## Step-09c: SNP-distance matrices & outlier exclusion

Computes pairwise SNP distances per serovar and automatically flags/excludes
extreme outliers. A **login-node summary** script.

### Prerequisite

- Conda env (apha_phylo):
- Step-09b complete.

### What it does

- Runs `snp-dists` to produce raw and QC-filtered SNP-distance matrices.
- Auto-excludes any isolate with median SNP distance > 10,000 (logged).

### Outputs (`results/09_phylogeny/03_snp_distances/`)

- Raw and QC-filtered SNP-distance matrices per serovar.
- `03_excluded_isolates.tsv` — auto-excluded outliers, with their median SNP.

### Result on this dataset

| Serovar     | n    | SNP range | Note                                             |
| ----------- | ---- | --------- | ------------------------------------------------ |
| Typhimurium | 7    | 18–60     | tight cluster                                    |
| Enteritidis | 10   | 2–508     | ERR3843447 is an outlier                         |
| Infantis    | 9    | 3–305     | 746/758 a distant sublineage; pESI+ cluster 3–14 |
| Kentucky    | 5    | 11–80     | ERR9714962 auto-excluded (median 39,897)         |

The Infantis result is one line of the five-way pESI evidence: the two
pESI-negative isolates (SRR30149746, SRR30149758) form a distant sublineage,
while the pESI+ isolates cluster tightly at 3–14 SNP.

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
bash scripts/09c_phylo_summary.sh
```

### Notes

- Uses `set -uo pipefail` (not `-e`) — dropping `-e` avoids mid-loop aborts that
  the earlier version hit.
- Counts are forced to clean integers before arithmetic (`${var:-0}` +
  `tr -d '[:space:]'`) — a stray `"0\n0"` otherwise breaks `$(( ))`.
- Kentucky shows n=5 here (not 6) because ERR9714962 is auto-excluded at the
  SNP-distance stage; it may still appear as a flagged tip in the fig2 tree.

---



## Step-10: Long-read hybrid assembly

Hybrid-assembles the three Infantis isolates that have matched Nanopore reads, to
resolve the pESI megaplasmid that short reads alone leave fragmented (Step-06).
Runs as an **SGE array job** (`-t 1-3`, one long-read pair per task).

### Prerequisite

- Conda env (apha_longread):
- Step-00a ONT reads and Step-01a trimmed Illumina reads present.
- `samples/longread_pairs.tsv` (Step-00) with the three pairings.

### What it does

Per pair:

1. **awk read stats** — lightweight ONT read-length summary (replaces NanoPlot).
2. **Filtlong** — length/quality filter
   (`--min_length 1000 --keep_percent 90 --target_bases 500M`).
3. **Unicycler** — hybrid assembly (trimmed Illumina + filtered ONT).

### Outputs (`results/10_longread_hybrid/<short_acc>/`)

- `01_ont_readstats.txt` — ONT read metrics.
- `02_unicycler/assembly.fasta` — hybrid assembly.
- `<short_acc>.hybrid.fasta` — assembly copy for downstream typing.

### Run

```bash
# Load conda env
module load igmm/apps/anaconda/2023.03  		# or use a personal conda
conda env create -f envs/apha_longread.yml
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
qsub scripts/10_longread_hybrid.sh
```

### Result on this dataset

- **SRR30149748** — pESI reconstructed (MOB-suite: AC358, 314,988 bp, IncFIB,
  conjugative); assembly 5.04 Mb, 27 contigs. pESI is fragmented in the assembly
  but MOB-suite clusters the fragments into one replicon.
- **SRR30149746** — 0 plasmids; 4.70 Mb, 4 contigs.
- **SRR30149758** — 0 plasmids; 4.65 Mb, 17 contigs.

### Notes

- **VMEM fix:** NanoPlot was removed — its numpy/OpenBLAS backend reserved
  262–277 G of virtual memory on the large nodes and the job was killed at 137.
  Read metrics are done with awk instead, and the script exports
  `MALLOC_ARENA_MAX=2`, `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`,
  `MKL_NUM_THREADS=1` after activating the env.
- Resources: `h_vmem=32G × 2 slots`, `THREADS=2`, `h_rt=48:00:00`.
- ONT reads are single-file (`${RAW_DIR}/<long_acc>_1.fastq.gz`); Illumina are
  the trimmed pairs from Step-01a.

---

## Step-10a: Long-read assembly & plasmid summary

Summarises the three hybrid assemblies and their plasmid content, using
MOB-suite calls rather than contig-size heuristics. A **login-node summary**
script.

### Prerequisite

- Conda env (apha_longread):
- Step-10 complete.

### What it does

- Builds an assembly overview (contigs, size, plasmid count, largest replicon).
- Builds a per-plasmid detail table (cluster, size, replicon type, mobility).

### Outputs (`results/10_longread_hybrid/`)

- `01_assembly_overview.tsv` — isolate, long_acc, n_contigs, assembly_size_bp,
  n_plasmids, total_plasmid_bp, largest_replicon, remark (feeds fig4).
- `02_plasmids_detail.tsv` — isolate, plasmid_cluster, size_bp, rep_type,
  mobility.

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
bash scripts/10a_longread_summary.sh
```

### Notes

- Plasmid calls come from MOB-suite, not contig-size guesses — `total_plasmid_bp`
  is the sum across plasmids, while `size_bp` in the detail table is one plasmid.
- TSVs are pure numeric; the `remark` column is auto-generated with kb.

---

## Key finding: pESI carriage confirmed five independent ways

SRR30149748 carries the pESI megaplasmid; the other two long-read Infantis
isolates (SRR30149746, SRR30149758) do not. Five independent lines of evidence
agree:

1. **MOB-suite, short-read** (Step-06) — pESI cluster AC358 in 748, absent in
   746/758.
2. **Genome size** (Step-07) — 748 ≈ 5.0 Mb vs 746/758 ≈ 4.6–4.7 Mb.
3. **CDS count** (Step-07) — 748 ≈ 4662 vs 746/758 ≈ 4265–4342.
4. **Phylogeny** (Step-09c) — 746/758 form a distant Infantis sublineage.
5. **Hybrid assembly + MOB-suite** (Step-10) — AC358, ~315 kb, reconstructed in
   748, none in the others.

This convergence is the strongest single result in the project: a mobile,
ESBL-carrying megaplasmid identified and cross-validated by independent methods —
exactly the kind of call that matters for surveillance.

---



## Step-11: Data visualisation (R)

Renders the presentation figures from the Steps 04–10 outputs. A wrapper script
(`11_figures.sh`) runs five R scripts in turn on the **login node** (R is
single-threaded here). Each R script reads its inputs through the config
variables exported by `config.sh`.

### Prerequisite

- Conda env (apha_figures):
- Steps 04a, 05a, 06a, 07a, 08a, 09, 09b and 10a complete (the figures read
  their summary tables and treefiles).
- The five R scripts (`11a`–`11e`.R) present in `scripts/`.

### What it does

| Script                  | Figure | Shows                                              |
| ----------------------- | ------ | -------------------------------------------------- |
| `11a_tree_coregene.R`   | fig1   | whole-set core-gene tree, tips coloured by serovar |
| `11b_tree_perserovar.R` | fig2   | per-serovar SNP trees, outliers flagged            |
| `11c_amr_matrix.R`      | fig3   | AMR gene matrix, chromosome vs plasmid location    |
| `11d_pesi_evidence.R`   | fig4   | pESI evidence panel (size, CDS, plasmid content)   |
| `11e_pangenome_bar.R`   | fig5   | pan-genome structure + serovar-specific genes      |

The wrapper runs each script with `Rscript`, logging a status line per figure. If
one script fails it prints a `WARN` and continues, so a single broken input does
not lose the rest.

### Outputs (`results/11_figures/`)

- `fig1_coregene_tree.png` — core-gene phylogeny (32 isolates).
- `fig2_perserovar_trees.png` — four per-serovar SNP trees.
- `fig3_amr_matrix.png` — AMR carriage, chromosomal vs plasmid-borne.
- `fig4_pesi_evidence.png` — pESI carriage evidence (3 long-read Infantis).
- `fig5_pangenome.png` — pan-genome structure and serovar-specific genes.

### Inputs read (by figure)

- **fig1** — `09_phylogeny/01_core_gene_tree/core_tree.treefile`,
  `04_serotype_mlst/typing_summary.tsv` (serovar in column 5).
- **fig2** — `09_phylogeny/02_snippy/<Serovar>/iqtree/<Serovar>_snp.treefile`.
- **fig3** — `05_amr_virulence/03_amr_matrix_annotated.tsv`,
  `06_plasmid_mge/03_amr_gene_location.tsv` (location strip).
- **fig4** — `07_annotation/01_annotation_summary.tsv`,
  `10_longread_hybrid/01_assembly_overview.tsv`.
- **fig5** — `08_pangenome/gene_presence_absence.Rtab`,
  `08_pangenome/02_serovar_specific_genes.tsv`.

### Run

```bash
# Load conda env
module load igmm/apps/anaconda/2023.03
conda env create -f envs/apha_figures.yml
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script (runs all five R scripts in turn)
bash scripts/11_figures.sh
```

### Notes

- The R scripts read paths via `Sys.getenv()`, which only works because the
  wrapper sources `config.sh` first — do **not** run the R scripts standalone.
- **fig3 shows 19 isolates, not 32** — only isolates with ≥1 acquired AMR gene or
  resistance mutation appear; the 10 Enteritidis are clean and drop out. State
  this on the slide so the omission is not mistaken for missing data.
- **fig2 per-panel scales differ** — each serovar tree has its own x-axis
  (Typhimurium ~0.008, Infantis ~0.003, Kentucky ~0.03), so outlier branch
  lengths are not directly comparable across panels.
- **fig5 Infantis-71** — the serovar-specific count is deflated by the pESI split
  (Step-08a); flag it as a threshold artefact when the figure is shown.
- A Fontconfig warning on the headless node is harmless — the PNGs write fine.

---



## Step-12: Final aggregated QC (MultiQC)

Pulls every QC output across the run into one MultiQC report — FastQC (raw and
trimmed), fastp, QUAST and BUSCO — with the CheckM2 table copied in alongside.
A **single job** on the login node; run last, after everything else is done.

### Prerequisite

- Conda env (apha_wgs):
- Steps 01, 01a, 01b and 03a complete (the QC outputs to aggregate).

### What it does

- Scans the FastQC (raw + trimmed), fastp, QUAST and BUSCO output directories
  (skipping any that are absent, so a partial run does not abort).
- Runs MultiQC over them into one report.
- Copies the CheckM2 `quality_report.tsv` into the report directory (CheckM2 is
  not a MultiQC module, so it is attached rather than parsed).

### Outputs (`results/12_multiqc_final/`)

- `multiqc_final_report.html` — the aggregated QC report.
- `checkm2_quality_report.tsv` — CheckM2 completeness/contamination, alongside.

### Run

```bash
# navigate to the project directory
cd "${PROJECT_DIR}"
# run script
bash scripts/12_multiqc_final.sh
```

### Notes

- Expect FastQC to appear twice (raw and trimmed, split by filename), fastp once,
  the single combined QUAST report covering all 32 assemblies, and BUSCO across
  all 32.
- The scan-only-if-present loop means the report still builds if one QC stage was
  skipped — missing directories are logged as `WARN`, not fatal.

---



## Concluding remarks

This pipeline takes public *Salmonella* reads through to an AMR-characterised,
serovar-placed set of isolates — from raw FASTQ to the resistance profile of each
genome, where that resistance sits, and how the isolates relate within and across
serovars. The design follows a surveillance question rather than a generic
assembly workflow: identify the isolate, find its determinants, decide whether
they are chromosomal (stably inherited) or plasmid-borne (horizontally mobile),
and check each call against an independent method. The pESI result — a mobile,
ESBL-carrying megaplasmid identified and cross-validated five independent ways —
is the kind of finding that changes how a case is interpreted.

## Design principles

- **Single source of truth** — all paths, environment names and parameters live
  in `config/config.sh`; edit once, and every step follows.
- **Reproducibility** — each stage ships as a pinned conda environment
  (`envs/*.yml`), so the exact toolset rebuilds from the repo. The known-breaking
  dependencies are pinned deliberately (e.g. `pyhmmer=0.10.15` for Bakta).
- **Independent cross-checks** — serovar is called by both SISTR and SeqSero2;
  plasmids by both MOB-suite and PlasmidFinder; AMR genes by AMRFinderPlus and
  ABRicate. The pESI carriage call rests on five orthogonal signals. A single
  tool is never the sole basis for a conclusion.
- **What vs where** — Step-05 records *what* resistance is present; Step-06
  records *where* it sits. Keeping these separate is what lets the AMR figure
  distinguish chromosomal SGI1 from plasmid-borne pESI at a glance.
- **Honest reporting** — the coverage-failed isolate is documented, not silently
  dropped; the deflated Infantis serovar-specific count is flagged as a threshold
  artefact wherever it appears.

## Reproducing the full run

1. Edit `config/config.sh` (paths, conda env names).
2. Create the conda environments from `envs/*.yml`.
3. Run the setup steps once: **Step-00** (sample lists), **Step-00a** (raw
   reads), **Step-00b** (references), **Step-00c** (databases).
4. Submit Steps **01 → 12** in order, matching each array job's `-t 1-N` range to
   the isolate count (`wc -l < samples/short_read_ids.txt`). Test any array with
   `-t 1-2` first.
5. Run each step's summary helper (`a`/`b`/`c`) after its main step; run the
   figures (Step-11) and final MultiQC (Step-12) last.

## Adapting to your data

- **Different isolates** — replace `samples/accessions.txt` and re-run Step-00;
  nothing else needs changing. Update every `-t 1-N` range to the new count.
- **Different serovars** — swap the four `REF_<SEROVAR>` references in
  `config.sh` and re-run the phylogeny branch (Steps 09a–09c).
- **No long reads** — skip Step-10; the short-read pESI evidence (Steps 06, 07,
  09) still stands on its own.
- **Thresholds** — the coverage gate (30×), the CDS flag (4200) and the SNP
  outlier cutoff (10,000) are set at the top of their respective scripts.

## Citation & contact

If you use or adapt this pipeline, please cite the underlying tools:
FastQC, fastp, MultiQC, Shovill/SPAdes, QUAST, BUSCO, CheckM2, mlst, SISTR,
SeqSero2, AMRFinderPlus, ABRicate, PlasmidFinder, MOB-suite, Bakta, Panaroo,
Snippy, Gubbins, IQ-TREE, snp-dists, Filtlong and Unicycler — together with the
reference and database resources (ENA/NCBI, CARD, ResFinder, VFDB, PlasmidFinder,
enterobacterales_odb10).

**Author:** Md Ataul Goni Rabbani (<m.a.g.rabbani@roslin.ed.ac.uk>),
The Roslin Institute, University of Edinburgh, UK.