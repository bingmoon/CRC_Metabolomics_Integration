# Method-Optimized Integration of Fecal Metabolomics and Transcriptomics in CRC

This repository contains the data and analytical scripts required to reproduce the metabolomic findings presented in the manuscript:
**"Method-Optimized Reanalysis of Fecal Metabolomics and TCGA-GTEx Transcriptomics Identifies Parallel Dysregulation of Three Metabolic Axes in Colorectal Malignant Transformation."**

## 📁 Repository Structure

* `data/`: Contains the processed LC-MS feature matrices for both negative and positive ion modes (derived from Metabolomics Workbench Study ID: ST002787).
* `scripts/`: Contains the R pipelines used for statistical modeling, biomarker discovery, and ontology-guided stratification.
* `results/`: Contains the computationally extracted list of host-endogenous metabolites (`.txt`) used for downstream KEGG pathway enrichment in MetaboAnalyst 5.0.

## ⚙️ Prerequisites and Installation

To execute the pipelines, R (version 4.0.0 or higher) is required. Please ensure the following R packages are installed:
```R
install.packages(c("jsonlite", "dplyr", "stringr", "tibble", "ggplot2", "scales", "pROC", "ggVennDiagram", "ggrepel", "pheatmap", "readr"))
BiocManager::install(c("impute", "ropls"))
