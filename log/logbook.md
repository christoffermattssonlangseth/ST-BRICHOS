# ST BRICHOS
## 2025-06-17
Today we focused on processing and analyzing spatial transcriptomics data from multiple mouse brain sections to investigate the effects of BRICHOS treatment in an Alzheimer’s disease model. We began by importing SpaceRanger outputs, handling variability in tissue position metadata, and converting filtered feature-barcode matrices into individual .h5ad files using Scanpy. Each sample was then tagged with its corresponding sample ID and spatial metadata.

Once preprocessed, we merged all individual .h5ad files into a single combined AnnData object. We ensured proper harmonization of the ```.uns["spatial"]``` key and corrected for axis flipping by modifying the spatial coordinate columns. Sample-specific spatial plots were generated for quality control, and Leiden clustering was performed after standard normalization, PCA, and neighborhood graph construction. 

## 2025-06-18
We then ran differential expression analyses to compare wild-type mice versus PBS-treated APPNLGF mice (baseline pathology) and BRICHOS-treated versus PBS-treated mice (treatment effect). Volcano plots were generated, and the most significant genes were annotated. Differentially expressed genes were also cross-referenced with known plaque-induced genes from literature to interpret treatment relevance.

We developed multiple visualizations, including spatial gene expression maps, treatment-specific cluster distributions, dot plots, and volcano plots with highlighted AD-related genes. To streamline the workflow, we modularized reusable functions into a Python script (analysis_utils.py) and configured autoreload in Jupyter notebooks. The working directory was organized with separate folders for notebooks, utility functions, and results.

We also generated project documentation in the form of a README.md that includes experimental background, design, and objectives. 

Work concluded with preparation for regional DE comparisons, module scoring for known gene sets, and further modularization of the analysis codebase.

I added the regional annotations that were made by Cloé. 