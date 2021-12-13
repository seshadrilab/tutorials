# Flow_Analysis_Tutorial

A tutorial for analyzing flow cytometry data in R using COMPASS (Combinatorial Polyfunctionality Analysis of Single Cells). 

# Introduction

This repository contains a basic workflow for analyzing flow cytometry data using COMPASS. COMPASS is a statistical framework that allows for unbiased analysis of antigen-specific T-cell subsets. COMPASS uses a Bayesian hierarchial framework to model all obeserved cell-subsets and select the most likely to be antigen-specific while regularizing the small cell counts that often arise in multi-parameter space. The model gives a posterior probability of specificity for each cell-subset and each sample, which can be used to profile a subject's immune response to external stimuli (e.g., infection or vaccination). 

For this tutorial, we will be analyzing a batch of data from the HAARVI cohort. For more information about COMPASS, please refer to the original manuscript and documentation:
* https://www.nature.com/articles/nbt.3187
* https://bioconductor.org/packages/COMPASS/

# Installation

The following R packages are required for this tutorial:
* [flowCore](	https://bioconductor.org/packages/flowCore/)
* [flowWorkspace](https://bioconductor.org/packages/flowWorkspace/)
* [COMPASS](https://bioconductor.org/packages/COMPASS/)

To install the R packages, open an R session and enter the following command lines:
```R
install.packages("BiocManager")
BiocManager::install("flowCore")
BiocManager::install("flowWorkspace")
BiocManager::install("COMPASS")
```

# Workflow Overview
1. Load data
2. Create a GatingSet
3. Create a COMPASSContainer
4. Run COMPASS
5. Visualize
6. Perform UMAP
7. Troubleshooting
