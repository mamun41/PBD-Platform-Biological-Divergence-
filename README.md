# PBD-Platform-Biological-Divergence-
# Project Name

# Bulk RNA-seq versus single-cell RNA-seq cros-platform divergence and similarity analysis
Author: Md Mamunur Rashid

This repository contains a complete, manuscript-ready R pipeline for integrating bulk RNA-seq and pseudo-bulk (single-nucleus RNA-seq–derived) expression data.
The workflow systematically quantifies platform-driven technical bias, disentangles it from biological variation, and validates findings at the single-cell level.
The pipeline was developed and applied to human adipose tissue (hSAT and hVAT) data and is fully generalizable to other matched bulk–single-cell datasets.

---

## 🚀 Overview

**Project Name** is a [brief explanation: e.g., *machine learning framework*, *bioinformatics pipeline*, *statistical method*] designed to:

- ✔ Key Objectives
- ✔ Target use case (e.g., bulk vs single-cell)
- ✔ Integrate matched bulk RNA-seq and pseudo-bulk RNA-seq data
- ✔ Quantify platform bias vs biological variance using variance partitioning
- ✔ Define a Platform Bias / Biological Divergence (PBD) score
- ✔ Identify platform-sensitive genes using data-driven peak detection
- ✔ Validate bias signatures assesment using: Correlation structure, PCA / sPCA, Differential expression (limma), Gene set enrichment (GSVA), mapping platform-biased genes back to single-cell resolution

---

## ✨ Key Features

- 🔬 Feature 1 (e.g., multi-omics integration)
- ⚡ Feature 2 (e.g., scalable to large datasets)
- 📊 Feature 3 (e.g., interpretable outputs)
- 🧩 Feature 4 (e.g., modular design)

---
snRNA-seq data can be obtained from: https://doi.org/10.5281/zenodo.18000778



# PBD: Platform Biological Divergence Analysis

This repository provides the code and framework for identifying **Platform Biological Divergence (PBD)** in adipose tissue transcriptomics, specifically comparing Bulk RNA-seq and single-nucleus Pseudobulk.

## 📊 Data Availability
The snRNA-seq Seurat objects and count matrices used in this analysis are available at:
**Zenodo**: [https://doi.org/10.5281/zenodo.18000778](https://doi.org/10.5281/zenodo.18000778)

## 🚀 Execution
1. Clone this repo.
2. Download the data from Zenodo into a folder named `hATdata/`.
3. Run `LMM_hAT_cStudy_nc.R`. The script will automatically install any missing R packages.

## ✒️ Author
**Md Mamunur Rashid** (2025)
** Email <mamun.stat92@gmail.com> **


### Install from GitHub

```bash
git clone https://github.com/USERNAME/REPO.git
cd REPO
