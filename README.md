# Flow_Analysis_Tutorial

A tutorial for analyzing flow cytometry data in R using COMPASS (Combinatorial Polyfunctionality Analysis of Single Cells). 

# Introduction

This repository contains a basic workflow for analyzing flow cytometry data using COMPASS. COMPASS is a statistical framework that allows for unbiased analysis of antigen-specific T-cell subsets. COMPASS uses a Bayesian hierarchial framework to model all obeserved cell-subsets and select the most likely to be antigen-specific while regularizing the small cell counts that often arise in multi-parameter space. The model gives a posterior probability of specificity for each cell-subset and each sample, which can be used to profile a subject's immune response to external stimuli (e.g., infection or vaccination). 

For this tutorial, we will be analyzing a batch of data from the HAARVI cohort. For more information about COMPASS, please refer to the original manuscript and documentation:
* https://www.nature.com/articles/nbt.3187
* https://bioconductor.org/packages/COMPASS/

# Installation

The following R packages are required for this tutorial:
* [here](https://cran.r-project.org/package=here) (not necessary but makes file referencing easier)
* [CytoML](https://bioconductor.org/packages/CytoML/)
* [flowCore](https://bioconductor.org/packages/flowCore/)
* [flowWorkspace](https://bioconductor.org/packages/flowWorkspace/)
* [COMPASS](https://bioconductor.org/packages/COMPASS/)

To install the R packages, open an R session and enter the following command lines:
```R
install.packages("here")
install.packages("BiocManager")
BiocManager::install("CytoML")
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

## Load data
```R
library(here)
library(CytoML)
library(flowCore)
library(flowWorkspace)
library(COMPASS)
```
The .xml and .fcs files for this dataset are stored in the Seshadri Lab shared drive (LSR Fortessa/2021 Summer HAARVIVAC/20210902 HAARVIVAC Batch 4/).
Download the .xml file and folder containing the .fcs files and drag into the "data" folder of this project directory.
```R
# Location of XML file
xml_path <- here::here("data/20211014_HAARVIVAC_B4V3_JP.xml")

# Location of .fcs files
fcs_subfolder <- here::here("data/20210902_HAARVIVAC_FCS_B4/")
```

Create a flowjo_workspace object with the function open_flowjo_xml().
```{r}
ws <- open_flowjo_xml(xml_path)
```

## Create a GatingSet
### Set-up
A GatingSet holds a set of GatingHierarchy objects, representing a set of samples and the gating scheme associated with each.
Look at the workspace metadata to choose which keywords to extract into the GatingSet. The flowjo_to_gatingset() function parses a flowJo Workspace to generate a GatingSet object.
```R
names(fj_ws_get_keywords(ws, 117)) 
keywords2import <- c("EXPERIMENT NAME",
                       "$DATE",
                       "SAMPLE ID",
                       "PATIENT ID",
                       "Stim",
                       "WELL ID",
                       "PLATE NAME") 
sampleGroup <- "Samples"

gs <- flowjo_to_gatingset(ws,                                    
                          name = sampleGroup, 
                          keywords = keywords2import,
                          path = fcs_subfolder, 
                          extend_val = -10000)
```
