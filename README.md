# Master's Thesis - DRIPc-seq analysis of R-loops
![R](https://img.shields.io/badge/R-4.6.0-276DC?logo=r&logoColor=white)
![Bioconductor](https://img.shields.io/badge/Bioconductor-3.23-1f65b7)
![renv](https://img.shields.io/badge/reproducible-renv-2C8EBB)
![conda](https://img.shields.io/badge/conda-environment-44A833?logo=anaconda&logoColor=white)
![Genome](https://img.shields.io/badge/genome-hg38-orange)
![License](https://img.shields.io/badge/License-GPLv3-blue)

Differential analysis and genomic characterisation of R-loops (DRIPc-seq) in K562 cells under a control and two knockdown conditions (`siUAP56`, `siBRG1`). This repository holds the Shell and R pipeline used for the preprocessing of data, differential, annotation, features, chromosomal-distribution and metagene analyses, plus the figure generation, for my Master's thesis.

> The `data/` folder is **not** included in the repository because of its size. See [Data](#data) for the expected layout.

## Overview

- **`bash/`** holds the Shell scripts used in preprocessing the data and **saves its results to `results/`**.
- **`src/`** holds all the reusable functions. Nothing is computed there; it only defines functions that the scripts `source()`.
- **`scripts/`** holds all scripts used in R and **saves its results to `results/`**.

## Repository structure
 
```
├── bash/                  # Preprocessing (shell): QC, alignment,
│                          #   peak calling and generation of the stranded
│                          #   BAM/BED files consumed by the R pipeline
├── data/                  # NOT in the repo (too large) — see "Data"
├── results/               # Generated outputs (.rds, .csv) and figures (images/)
├── scripts/               # Analysis pipeline (run in order, see "Running")
├── src/                   # Reusable functions sourced by the scripts
├── .Rprofile              # Activates the renv environment on project open
├── .gitignore
├── LICENSE
├── README.md
├── TFM_DRIPc-seq.Rproj    # RStudio project (sets the working directory)
├── conda_TFM.txt          # Conda environment for the command-line tools (bash/)
└── renv.lock              # Exact R package versions (renv)
```

### `scripts/`
 
| Script | What it does |
|---|---|
| `01_differential_analysis.R` | Differential R-loop regions per condition and strand with DiffBind (DESeq2) vs the `siC` control; ChIPseeker annotation to genes; sense / antisense assignment. Saves the gene lists and annotation objects. |
| `02_features_analysis.R` | Gene features (GC content, gene length, exon number, biotype) via biomaRt; normality (Shapiro–Wilk), Wilcoxon tests and Cliff's delta effect sizes. |
| `03_physical_mapping.R` | Fraction of affected genes per chromosome, correlation with chromosomal GC content (Spearman), long-gene control and the intronic profile. |
| `04_metaplot.R` | Builds the metagene table: relative position of each R-loop across its host gene (TSS → TTS) with flanking regions. |
| `05_figures.R` | The only plotting script. Reads the saved results and renders every figure/panel: Venn, compositional, localisation/metagene, chromosomal, and long-gene specificity. |

### `src/`
 
| File | Contents |
|---|---|
| `config.R` | Paths, genome annotation (hg38 `TxDb` / `org.Hs.eg.db`) and the shared palette, group order and `theme_tfm()` used by every figure. Sourced first by every script. |
| `differential_annotation_functions.R` | DiffBind sample sheets and run, BED / DiffBind-result loading, ChIPseeker annotation, sense/antisense assignment and gene-list helpers. |
| `features_functions.R` | Statistical functions: Shapiro–Wilk, Wilcoxon-Mann-Whitney, and Cliff's delta. |
| `mapping_functions.R` | SYMBOL ↔ ENTREZ conversion, per-chromosome gene fractions, chromosomal GC content and long-gene annotation. |
| `metagene_functions.R` | Relative-position computation and assembly of the metagene table. |
| `figures_functions.R` | All plotting functions (feature panels, metagenes, chromosome plots and annotation bars). |

## Data
 
`data/` is not tracked (too large for GitHub). To reproduce the analysis,
regenerate it from the raw sequencing data with the `bash/` preprocessing. The R
pipeline expects, under `data/`:
 
```
data/
├── bam_stranded/   # stranded BAMs (fwd/rev), including the RNase H control (RNH)
├── peaks/          # broadPeak peak calls per sample and strand
└── bed/            # BED files of regions per condition (e.g. siC_fwd.bed)
```
 
The exact sub-paths are defined at the top of each script (`BAMDIR`, `PEAKDIR`,
`BEDDIR`, …) and ultimately under the `root` set in `src/config.R`.
 
> **Important:** `src/config.R` defines an absolute `root` path that points to my machine. Also, Shell scripts define the same absolute path. Edit them to your own location before running anything.

## Running

With the GEO data downloaded, run Shell and R scripts in order to reproduce the experiment after using the Conda environment list and renv.lock files to reproduce my enviroments.

## Author & license
 
José Livan Vargas Castro — Master's Thesis.
Released under the terms in [`LICENSE`](LICENSE).
