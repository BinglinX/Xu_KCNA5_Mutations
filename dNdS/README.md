# dN/dS Analysis of Potassium Channel Genes in Colorectal Adenocarcinoma

## Overview

This directory contains the pipeline for calculating dN/dS ratios for potassium
channel genes using COSMIC mutation data, to identify genes under positive selection
in colorectal adenocarcinoma (COAD). A permutation/binomial test assesses whether
KCNA5 mutations cluster in functionally important regions beyond random expectation.

## Directory structure

```
1_dNdS/
├── data/
│   ├── cosmic_LI_AC_filtered.tsv       # All columns, filtered to K+ channel genes only
│   ├── cosmic_LI_AC_for_dnds.tsv       # Whole genome, 5 columns for dndscv input
│   └── K_channels_LI.txt              # Potassium channel genes expressed in large intestine
│
├── filter.sh                           # Extracts 5 columns from COSMIC TSV for dndscv input
├── filter_gene.sh                      # Filters COSMIC TSV to genes in K_channels_LI.txt
│
├── cosmic_dnds_calculation.R           # Runs dndscv on COSMIC COAD data
├── cosmic_coad_extraction_dnds.ipynb   # Extracts and visualises K+ channel dN/dS results
│                                       # clustering in T1 domain
└── results/
    └── dnds_coad_cosmic.csv         # dndscv output: per-gene dN/dS estimates
```

## Data

Raw COSMIC data is not included due to licensing restrictions.

Download from: https://cancer.sanger.ac.uk/cosmic/download  
Select: Genome Screens Mutant → filter tissue with the filter: "large intestine", "include all", "carcinoma", "adenocarcinoma". Save as `data/cosmic_LI_AC.tsv`.

The three processed files are derived from this as described under Reproduction below.

## Reproduction

### Step 1 — Prepare input files

Extract the 5 columns required by dndscv (whole genome):
```bash
./filter.sh data/cosmic_LI_AC.tsv data/cosmic_LI_AC_for_dnds.tsv
```

Filter to potassium channel genes expressed in large intestine (for inspection):
```bash
./filter_gene.sh data/cosmic_LI_AC.tsv data/K_channels_LI.txt data/cosmic_LI_AC_filtered.tsv
```

### Step 2 — Run dndscv (R)

```r
source("cosmic_dnds_calculation.R")
```

This deduplicates mutations per sample, formats the input for dndscv, runs the
model against hg38, and writes per-gene dN/dS estimates to
`results/dnds_coad_cosmic.csv`.

Key parameters:
- Reference genome: hg38
- `max_muts_per_gene_per_sample = Inf` (no cap per gene)
- `max_coding_muts_per_sample = Inf` (no cap per sample)

### Step 3 — Extract and visualise results (Python)

Open `cosmic_coad_extraction_dnds.ipynb`.

This notebook:
- Loads `results/dnds_coad_cosmic.csv`
- Filters results to a panel of ~60 potassium channel genes
- Plots log2(dN/dS) per gene as a bar chart
- Highlights genes with significant missense selection (pmis_cv < 0.05)
- KCNA5: p = 1e-4; KCNK12: p = 5e-4


## Dependencies

### R
- dndscv
- dplyr

### Python
- pandas, numpy, matplotlib, seaborn, scipy

## Key output

`results/dnds_coad_cosmic.csv` — key columns:
- `gene_name`: Hugo symbol
- `n_syn` / `n_mis`: observed synonymous / missense mutation counts
- `wmis_cv`: missense dN/dS estimate
- `pmis_cv`: p-value for missense selection
