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

---

## 🚀 Step-by-Step Reproduction Guide

To ensure absolute computational reproducibility, the analytical workflow is modularized into four numerical scripts. Please execute them in sequence.

### **[Script 01] Luminal Metabolomic Profiling** (`01_Metabolomics_Pipeline.R`)
*   **Function:** Ingests the `st002787` targeted metabolomics matrices. Performs Total Ion Current (TIC) normalization, KNN imputation, and Pareto scaling. Executes non-parametric Wilcoxon rank-sum tests and strict 200-iteration permutation-validated PLS-DA.
*   **Paper Correspondence:** Generates the core statistical matrices that underpin **Figure 1** (Global PCA), **Figure 2** (Venn Intersections), and **Figure 3** (Volcano Plots & Heatmaps) representing the transition from Colorectal Adenoma (CRA) to Cancer (CRC).

### **[Script 02] Ontology Stratification** (`02_Ontology_Pathway_Analysis.R`)
*   **Function:** Computationally cross-references significant metabolites against HMDB v5.0 to segregate "Host-endogenous" from "Microbiota-dependent" pools. Outputs the filtered `Step2_Host_Metabolites_for_MetaboAnalyst.txt` directly into the `results/` folder for Over-Representation Analysis (ORA).
*   **Paper Correspondence:** Underpins the biological rationale in **Figure 4** (Host KEGG Pathway Enrichment) and identifies the specific collapse of the ceramide and Kynurenic acid axes.

### **[Script 03] Deep-Tissue Transcriptomic Mapping** (`03_Clinical_Cox_Regression.R`)
*   **Function:** Analyzes TCGA-COAD/READ cohorts matched against true healthy intestinal tissues from GTEx (avoiding adjacent *field cancerization* noise). Constructs a multivariate Cox proportional-hazards regression model adjusting for Age, Gender, and TNM stage.
*   **Paper Correspondence:** Reproduces the genetic expression boxplots (**Figure 5**, identifying *NR1H4* silencing, *ASAH1* and *IDO1* upregulation) and the Kaplan-Meier/Multivariate Cox Forest plots (**Figure 6**).

### **[Script 04] Machine Learning Diagnostic Panel** (`04_ML_Diagnostic_Panel.R`)
*   **Function:** Integrates orthogonal features from the three verified metabolic axes: Tryptophan (*IDO1*), Ceramide (*ASAH1*), and Bile Acid (FXR). Employs 5-fold cross-validated Random Forest and multivariable logistic regression.
*   **Paper Correspondence:** Generates the ROC curves and DeLong's test metrics (**Figure 7**), confirming the predictive robustness of the 3-axis panel (AUC = 0.88).

---

## 🔒 Methodological Rigor & Dependencies

*   **Anti-Overfitting Protocols:** All supervised PLS-DA models enforce a strict **200-iteration permutation test** to guarantee biological authenticity. 
*   **Random Seeds:** `set.seed(123)` is hardcoded across all Random Forest, KNN imputation, and data-splitting steps to ensure users obtain the *exact* identical statistical metrics (e.g., P-values, AUCs) as published.
*   **Dependencies:** The pipeline was developed in **R (v4.3.0)**. Core packages include `ropls`, `randomForest`, `survival`, `pROC`, `ggplot2`, and `ggrepel`.

## ✉️ Contact & Citation
If you find this code or our multi-axis integration strategy helpful for your research, please cite our manuscript. For any computational inquiries or bug reports, please open an issue in this repository.
