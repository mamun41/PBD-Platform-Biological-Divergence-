# PBD: Platform Biological Divergence Analysis 🧬
  
## 📌 Project Overview
  
The **PBD (Platform Biological Divergence)** framework is designed to quantify and visualize transcriptomic discrepancies between **Bulk RNA-seq** and **single-nucleus/cell Pseudobulk (PB)** data.

In cross-platform studies, sequence technology driven signals often masks biological signals. PBD allows researchers to identify platform-sensitive genes and evaluate the structural concordance of their datasets using advanced statistical modeling, including **Variance Partitioning**, **Differential expression analysis** and **Sparse PCA**.

---
  
## 💾 Data Availability
  
The datasets used in this study (Human Adipose Tissue - hAT) are publicly archived on Zenodo.

* **snRNA-seq Data**: [https://doi.org/10.5281/zenodo.18000778](https://doi.org/10.5281/zenodo.18000778)
* **Required Files**: To run the analysis, place the following files into a folder named `/hATdata` in your project root:
  * `sc_hSAT_seurat.rds`
* `sc_hVAT_seurat.rds`
* `bulk_count_data.rds`
* `pb_count_data.rds`
* `metaInfo.rds`


---
  
## 🚀 Key Features
  
* **Normalization**: Integrated TMM and Log2-CPM scaling for cross-platform comparison.
* **PBD Scoring**: A metric derived from **Variance Partitioning** to separate platform bias from biological variance.
* **Divergence Discovery**: Uses **Sparse PCA (sPCA)** to isolate primary platform effects.
* **Interpretation**: Gene Ontology (GO) enrichment for divergent gene sets and Cell cluster mapping annotation.

---
  
  ## 🛠 Installation & Quick Start
  
  ### 1. Clone the Repository
  
  ```bash
git clone https://github.com/mamun41/PBD-Platform-Biological-Divergence-.git
cd PBD-Platform-Biological-Divergence-
  
  ```

### 2. Execution

Open the main analysis script `LMM_hAT_cStudy_nc.R` in RStudio.

**Note on Dependencies**: The utility script `pbulk_nc_utils.R` is configured to **automatically detect and install** all missing R packages from both CRAN and Bioconductor upon execution. No manual installation is required.

---
  ## 📈 Pipeline Workflow
  
1. **Preprocessing**: Harmonization of common genes and normalization across platforms.
2. **Global Diagnostic**: PCA and Pearson correlation to visualize platform separation.
3. **Variance Decomposition**: Quantifying the contribution of "Platform", "Sample", and "Tissue" to gene expression variance.
4. **Thresholding**: Automatic bimodal peak detection to identify **High-PBD** (platform-sensitive) genes.
5. **Validation**: Single-nucleus marker validation and volcano plot integration.

---
## 📊 Expected Outputs
  
The analysis generates a `/hAT_fig` directory containing:
  
* **PCA Plots**: Global structure comparison between platforms.
* **PBD Volcano Plots**: Multi-method visualization of divergent genes.
* **GO Enrichment**: Functional pathways affected by platform bias.
* **Variance Partitioning Boxes**: Comparison of variance proportions between Low and High PBD genes.

---
  ## ✒️ Citation
  
  If you utilize this framework or data in your research, please cite:

---
  
## ✉️ Contact
**Author:** Md Mamunur Rashid <mamun.stat92@gmail.com>

**Project:** Platform Biological Divergence (PBD) Study

**Date:** December 2025


