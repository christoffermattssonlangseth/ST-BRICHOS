# 🧠 BRICHOS Spatial Transcriptomics
<p align="center">
  <img src="assets/logo.png" alt="ST BRICHOS Logo" width="400"/>
</p>

## Project Overview

This repository contains analysis code and documentation for a spatial transcriptomics study investigating the effects of **BRICHOS treatment** in the **APP<sub>NLGF</sub> transgenic mouse model** of Alzheimer’s disease (AD). The project utilizes **SPATIAL TRANSCRIPTOMICS** technology to explore how BRICHOS modulates gene expression at a regional and cellular level across the mouse brain.

---

## 🧪 Experimental Design

| Group              | Treatment             | n |
| ------------------ | --------------------- | - |
| WT                 | Wild Type (untreated) | 1 |
| APP<sub>NLGF</sub> | PBS (non-treated)     | 3 |
| APP<sub>NLGF</sub> | BRICHOS (treated)     | 3 |

* **Treatment schedule**: 2×/week for 3 months
* **Platform**: Spatial Transcriptomics 10x Genomics
* **Tissue**: Coronal mouse brain sections, aligned to the Allen Mouse Brain Atlas

---

## 🎯 Objectives

1. **Baseline Gene Expression**
   Compare **WT vs. APP<sub>NLGF</sub> non-treated** to identify baseline AD-related transcriptional changes.

2. **Validation of AD Signatures**
   Cross-reference differentially expressed genes with published datasets to identify known **Alzheimer’s-associated genes**.

3. **Effect of BRICHOS Treatment**
   Compare **APP<sub>NLGF</sub> PBS vs. BRICHOS-treated** groups to identify BRICHOS-modulated genes globally and regionally.

4. **Intersection with AD Pathology**
   Determine which BRICHOS-affected genes overlap with validated AD-associated markers (e.g. **plaque-induced genes**, microglial activation, etc.).

---

## 📁 Repository Structure

```
├── data/                    # Raw and processed data (e.g. .h5ad files)
├── notebooks/               # Jupyter Notebooks for analysis
│   ├── preprocessing.ipynb  # Sample merging, QC, spatial setup
│   ├── clustering.ipynb     # Dimensionality reduction and Leiden clustering
│   ├── differential_expression.ipynb
│   └── visualization.ipynb
├── utils/                   # Helper Python functions (e.g., plotting, scoring)
└── README.md
```

---

