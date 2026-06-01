# Coral Reef Microbial Ecology Analysis Scripts

This repository contains analysis scripts used for coral reef microbial community analyses, including alpha diversity, ordination, co-occurrence networks, community assembly, environmental association, differential abundance, feature-ASV validation, Random Forest screening, and MAG/KO functional visualization.

The scripts are provided as reproducible analysis templates. Raw data and large generated results are intentionally excluded from the repository.

## Repository structure

```text
.
├── README.md
├── .gitignore
├── requirements_R.txt
├── requirements_python.txt
├── data/
│   └── README.md
├── results/
│   └── README.md
└── scripts/
    ├── R/
    └── python/
```

## Script index

| Script | Purpose |
|---|---|
| `scripts/R/01_alpha_diversity.R` | Alpha diversity analysis |
| `scripts/R/02_cca_envfit.R` | CCA and envfit analysis |
| `scripts/R/03_network_analysis.R` | Co-occurrence network analysis |
| `scripts/R/04_network_environment_correlation.R` | Environmental effects on network nodes and edges |
| `scripts/R/05_neutral_community_model.R` | Neutral community model |
| `scripts/R/06_bnti_qpen_assembly.R` | βNTI/RCbray community assembly |
| `scripts/R/07_rcbray_environment_effects.R` | βNTI/RCbray and environmental effects |
| `scripts/R/08_differential_abundance_volcano.R` | Differential abundance volcano plots |
| `scripts/R/09_ternary_plot.R` | Ternary plot and composition pies |
| `scripts/R/10_feature_asv_filtering.R` | Feature ASV filtering |
| `scripts/R/11_random_forest_asv_screening.R` | Random forest ASV screening and AUC validation |
| `scripts/R/12_feature_asv_health_prediction_validation.R` | Feature ASV health-state prediction validation |
| `scripts/R/13_environment_difference_near_off.R` | Near-reef vs off-reef environmental difference |
| `scripts/R/14_bray_curtis_piecewise_regression.R` | Bray-Curtis distance segmented regression |
| `scripts/R/15_ko_pathway_bubble_plot.R` | KO pathway bubble plot |
| `scripts/python/16_eigen_microstate_analysis.py` | Eigen microstate analysis |
| `scripts/python/17_ar1_entropy_analysis.py` | AR(1) entropy analysis |

## Data input format

Most scripts assume the following common microbiome data format:

- ASV/OTU abundance table: rows are ASVs/OTUs and columns are samples.
- Sample metadata table: rows are samples, with a `Group` column for health state or spatial group.
- Taxonomy table: rows are ASVs/OTUs and columns are taxonomic ranks.
- Phylogenetic tree: required for βNTI/iCAMP analyses.

Common file names used in the scripts include:

```text
asv_sampling.csv
sample_sampling.csv
tax_table.csv
filtered_tree.tree
otutable.csv
sample.csv
sample_Distance.csv
```

Before running a script, place the corresponding input files in your working directory, or edit the file paths inside the script.

## Installation

### R packages

Install the required R packages listed in `requirements_R.txt`. Many scripts automatically check and install missing CRAN packages, but specialized packages such as `DESeq2`, `phyloseq`, `ggtree`, `iCAMP`, and `ggClusterNet` may require Bioconductor or GitHub installation.

Example:

```r
install.packages(c("tidyverse", "vegan", "ggplot2", "ggpubr", "ranger", "caret", "pROC"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "phyloseq", "ggtree"))
```

### Python packages

Install Python dependencies with:

```bash
pip install -r requirements_python.txt
```

## Running the scripts

Each script is independent. In RStudio or an R terminal, set the working directory to the folder containing the required input data, then run the selected script.

Example:

```r
source("scripts/R/13_environment_difference_near_off.R")
```

For Python scripts:

```bash
python scripts/python/16_eigen_microstate_analysis.py
python scripts/python/17_ar1_entropy_analysis.py
```

## Main analysis modules

### Diversity and ordination

- `01_alpha_diversity.R`: alpha diversity analysis.
- `02_cca_envfit.R`: CCA and envfit analysis for environmental associations.

### Network analysis

- `03_network_analysis.R`: ggClusterNet-based co-occurrence network analysis, robustness/vulnerability, and Zi-Pi plots.
- `04_network_environment_correlation.R`: correlation heatmap between environmental factors and sample-level network nodes/edges.

### Community assembly

- `05_neutral_community_model.R`: neutral community model analysis.
- `06_bnti_qpen_assembly.R`: βNTI/RCbray community assembly analysis using iCAMP.
- `07_rcbray_environment_effects.R`: links assembly processes with environmental factors.

### Differential abundance and characteristic ASVs

- `08_differential_abundance_volcano.R`: DESeq2-based differential abundance and volcano plots.
- `09_ternary_plot.R`: ternary plots and regional composition pie charts.
- `10_feature_asv_filtering.R`: feature ASV filtering and enrichment analysis.
- `12_feature_asv_health_prediction_validation.R`: ASV-score-based health-state prediction validation.

### Machine learning

- `11_random_forest_asv_screening.R`: pairwise Random Forest screening, Top ASV ranking, health-state bias assessment, and AUC validation.

### Spatial and functional analyses

- `14_bray_curtis_piecewise_regression.R`: Bray-Curtis pairwise distance and near-distance segmented regression.
- `15_ko_pathway_bubble_plot.R`: KO pathway class level 2 bubble plot.

### Python analyses

- `16_eigen_microstate_analysis.py`: eigen microstate analysis based on transformed ASV/OTU abundance data.
- `17_ar1_entropy_analysis.py`: AR(1) analysis for mode contribution entropy.

## Notes for GitHub release

- Do not upload raw sequencing data, large intermediate matrices, or private local paths.
- Keep only small example input files if necessary.
- If a script contains a commented `setwd("path/to/your/data_directory")`, replace it locally before running, but avoid committing personal absolute paths.
- Output files such as `.csv`, `.pdf`, and `.png` are ignored by `.gitignore` by default.

## Citation

If these scripts are associated with a manuscript, cite the manuscript and relevant software packages used in each analysis.

## Contact

For questions, please contact the repository maintainer.
