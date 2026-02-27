# Single-Cell RNA-Seq Analysis: AD00201

This repository provides a pipeline to extract and visualize DNA Damage Response (DDR) signaling from Alzheimer’s Disease (AD) brain samples (e.g., AD00201) hosted on the **ssREAD** database. 

## Overview
We map cell-type-specific genomic stress across major neural populations including Excitatory/Inhibitory Neurons, Glia, and Vasculature. The analysis specifically targets the expression of DDR machinery:
* **Key Markers:** ATM, TP53, CHEK2, BRCA1, PARP1, CDKN1A, CDKN2A, RELA, ABL1.

## Visualizations
| Cell Type Clusters | DDR Gene Expression |
|---|---|
|![Cell Types](./plots/AD00201_human_celltypes.png)|
|![DDR Expression](./plots/AD00201_human_genes.png)|

## Methods & Pipeline
* **Hybrid Bridge:** Uses `rpy2` to bridge Seurat (R) and Pandas (Python) for seamless data transition.
* **Extraction:** Loads `.qsave` files, performing normalization, PCA, and UMAP (top 20 PCs) in R via the `qs` package.
* **Species Logic:** Filters Human vs. Mouse cells by comparing expression sums of species-specific gene nomenclature (e.g., all-caps `ATM` vs. title-case `Atm`).
* **Visualization:** Exports processed CSVs and UMAP coordinates for plotting with `matplotlib` and `seaborn`.

## Quick Start
1. **Environment:** Ensure R (with `Seurat`, `qs`) and Python are installed.
2. **Install Dependencies:** ```bash
   pip install rpy2 pandas matplotlib seaborn
3. Run Extraction: ```bash
   python extract_ssREAD.py
