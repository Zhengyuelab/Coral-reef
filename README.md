# Coral Reef Microbial Ecology Analysis Scripts

This repository contains analysis scripts used for coral reef microbial community analyses, including alpha diversity, ordination, co-occurrence networks, community assembly, environmental association, differential abundance, feature-ASV validation, Random Forest screening, segmented regression, eigen microstate analysis, and AR(1) analysis.

The scripts are provided as reproducible analysis templates. Raw data and large generated results are intentionally excluded from the repository.

## Repository structure

```text
.
├── README.md
├── Code/
│   ├── Alpha diversity_R.txt
│   ├── CCA_R.txt
│   ├── Network_Env_R.txt
│   ├── NCM_R.txt
│   ├── βNTI_R.txt
│   ├── RCbray_Env_R.txt
│   ├── Volcano_R.txt
│   ├── Ternary Plot_R.txt
│   ├── Feature ASV Filtering_R.txt
│   ├── Random forest_R.txt
│   ├── Feature ASV_R.txt
│   ├── Env_difference_Off_Near.txt
│   ├── Piecewise regression_R.txt
│   ├── EM_python_analysis.txt
│   └── AR1_python_analysis.txt
```

## Script index

| Script | Purpose |
|---|---|
| `Alpha diversity_R.txt` | Alpha diversity analysis |
| `CCA_R.txt` | CCA and envfit analysis |
| `Network_Env_R.txt` | Environmental effects on network nodes and edges |
| `NCM_R.txt` | Neutral community model analysis |
| `βNTI_R.txt` | βNTI-based community assembly analysis |
| `RCbray_Env_R.txt` | RCbray and environmental effects |
| `Volcano_R.txt` | Differential abundance analysis and volcano plots |
| `Ternary Plot_R.txt` | Ternary plot and composition visualization |
| `Feature ASV Filtering_R.txt` | Feature ASV filtering |
| `Random forest_R.txt` | Random Forest-based ASV screening and validation |
| `Feature ASV_R.txt` | Feature ASV-based health-state prediction validation |
| `Env_difference_Off_Near.txt` | Environmental difference analysis between near-reef and off-reef samples |
| `Piecewise regression_R.txt` | Bray-Curtis distance and segmented regression analysis |
| `EM_python_analysis.txt` | Eigen microstate analysis |
| `AR1_python_analysis.txt` | AR(1) analysis for entropy or mode-contribution trends |

## Data input format

Most scripts assume the following common microbiome data format:

- ASV/OTU abundance table: rows are ASVs/OTUs and columns are samples.
- Sample metadata table: rows are samples, with a `Group` column for health state or spatial group.
- Taxonomy table: rows are ASVs/OTUs and columns are taxonomic ranks.
- Phylogenetic tree: required for βNTI/iCAMP analyses.

Common file names used in the scripts may include:

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

Install the required R packages before running the R scripts. Commonly used packages include:

```r
install.packages(c(
  "tidyverse",
  "vegan",
  "ggplot2",
  "ggpubr",
  "ranger",
  "caret",
  "pROC",
  "segmented"
))
```

Some analyses may require additional packages such as `DESeq2`, `phyloseq`, `ggtree`, `iCAMP`, or `ggClusterNet`.

Example Bioconductor installation:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("DESeq2", "phyloseq", "ggtree"))
```

### Python packages

The Python scripts require common scientific-computing packages, including:

```text
numpy
pandas
scipy
matplotlib
```

Install with:

```bash
pip install numpy pandas scipy matplotlib
```

## Running the scripts

Each script is independent. In RStudio or an R terminal, set the working directory to the folder containing the required input data, then run the selected script.

Example for R:

```r
source("Code/Env_difference_Off_Near.txt")
```

Example for Python:

```bash
python Code/EM_python_analysis.txt
python Code/AR1_python_analysis.txt
```

If preferred, `.txt` files containing executable R or Python code can be renamed locally to `.R` or `.py` before running.

## Main analysis modules

### Diversity and ordination

- `Alpha diversity_R.txt`: alpha diversity analysis.
- `CCA_R.txt`: CCA and environmental fitting analysis.

### Network and environmental association

- `Network_Env_R.txt`: environmental association analysis for network nodes and edges.

### Community assembly

- `NCM_R.txt`: neutral community model analysis.
- `βNTI_R.txt`: βNTI-based community assembly analysis.
- `RCbray_Env_R.txt`: association between assembly-related metrics and environmental variables.

### Differential abundance and characteristic ASVs

- `Volcano_R.txt`: differential abundance analysis and volcano plots.
- `Ternary Plot_R.txt`: ternary plots and composition visualization.
- `Feature ASV Filtering_R.txt`: feature ASV filtering.
- `Feature ASV_R.txt`: feature ASV-based health-state prediction validation.

### Machine learning

- `Random forest_R.txt`: Random Forest-based ASV screening and validation.

### Spatial and time-series-related analyses

- `Env_difference_Off_Near.txt`: near-reef versus off-reef environmental factor analysis.
- `Piecewise regression_R.txt`: Bray-Curtis distance and segmented regression.
- `EM_python_analysis.txt`: eigen microstate analysis.
- `AR1_python_analysis.txt`: AR(1) analysis.

## Notes for GitHub release

- Do not upload raw sequencing data, large intermediate matrices, or private local paths.
- Keep only code files and small example input templates if necessary.
- If a script contains a commented `setwd("path/to/your/data_directory")`, replace it locally before running, but avoid committing personal absolute paths.
- Output files such as `.csv`, `.pdf`, and `.png` should generally be excluded from version control.

## Citation

If these scripts are associated with a manuscript, please cite the manuscript and the relevant software packages used in each analysis.

## Contact

For questions, please contact the repository maintainer.
