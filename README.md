# Integrative Reanalysis of Cross-Cohort Multi-Omics in Colorectal Cancer

[![R version](https://img.shields.io/badge/R-v4.3.0+-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Data: Metabolomics](https://img.shields.io/badge/Data-ST002787-orange.svg)](https://www.metabolomicsworkbench.org/)
[![Data: TCGA-GTEx](https://img.shields.io/badge/Data-TCGA%20%7C%20GTEx-purple.svg)](https://portal.gdc.cancer.gov/)

Welcome to the official computational repository for our manuscript: **"Integrative Reanalysis of Cross-Cohort Multi-Omics Reveals Synchronized Dysregulation of Three Metabolic Axes and Yields a Robust Diagnostic Panel for Colorectal Cancer"** (Under Review at *Journal of Translational Medicine*).

This repository serves as the transparent computational backbone of our study. While the final high-resolution figures are provided as manuscript attachments, this repository contains the raw curated datasets and the exact R scripts required to reproduce our statistical analyses, machine learning models, and cross-omics mappings from scratch.

---

## 📂 Repository Structure

The repository is streamlined into three main directories: `data/`, `scripts/`, and `results/`.

```text
CRC_MultiOmics_Integration/
├── data/                                         # Curated input datasets
│   ├── st002787_negative_clean.txt               # Negative ion mode LC-MS data
│   ├── st002787_positive_clean.txt               # Positive ion mode LC-MS data
│   ├── InputData_03_TCGA_CRC_ASAH1_TPM_Stratified.csv # Transcriptomic expression matrix
│   └── InputData_04_TCGA_CRC_Clinical_Survival.csv    # Clinical metadata for Cox regression
├── scripts/                                      # Sequential R execution scripts
│   ├── 01_Metabolomics_Pipeline.R                # Data prep, PLS-DA, and differential analysis
│   ├── 02_Ontology_Pathway_Analysis.R            # HMDB stratification to prevent ecological fallacy
│   ├── 03_Clinical_Cox_Regression.R              # Multivariate survival & transcriptomic mapping
│   └── 04_ML_Diagnostic_Panel.R                  # Random Forest & ROC panel construction
├── results/                                      # Computational outputs
│   └── Step2_Host_Metabolites_for_MetaboAnalyst.txt # Filtered host-endogenous pool list
└── README.md
