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

### QC
Make sure that the gating trees are consistent for all samples.
```R
pop_lists <- lapply(gs, gh_get_pop_paths)
unique(pop_lists)
```
Remove channels from flow data that are not used by gates.
```R
gs <- gs_remove_redundant_channels(gs) # drop SSC-H, V655-A, V570-A
```
Add names to all channels or fix their names.
```R
dput(unname(pData(parameters(gh_pop_get_data(gs[[1]])))[,2]))
markernames <- c("Time", "FSC-A", "FSC-H", "SSC-A", "CD8b", "TNFa", "CD107a",
                 "CD154", "CD3", "IL2", "CD4", "IL17a", "IL4_5_13", "CD14_19",
                 "CCR7", "CD38", "LD", "IFNg", "CD45RA", "HLADR")
names(markernames) <- pData(parameters(gh_pop_get_data(gs[[1]])))[,1]
markernames(gs) <- markernames
pData(parameters(gh_pop_get_data(gs[[1]])))[,c(1,2)]
```
Plot gating tree.
```R
plot(gs, fontsize=15, bool=T)
```

### Save GatingSet
This can take a while.
```R
if(!dir.exists(here::here("out/GatingSet"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/GatingSet")))
  dir.create(here::here("out/GatingSet"), recursive = T)
}
save_gs(gs, here::here("out/GatingSet"))
```

## Create a COMPASSContainer
A COMPASSContainer is the data structure used to hold data from an ICS experiment. The input for this code is a GatingSet or a GatingSetList. Counts, metadata, and 
single cell data are extracted and fed into the COMPASSContainer constructor.
```R
# A regular expression to match a single node in the gating tree
parent_node <- "CD4+"

# A character identifying the subject id column in the GatingSet metadata
id <- "SAMPLE ID"

# markermap contains the output of markernames(gs)
markernames(gs) 
markermap <- list("IL2", "IL4_5_13", "IFNg", "TNFa", "IL17a", "CD154", "CD107a")
names(markermap) <- paste0(parent_node, "/", c("IL2+", "IL4_5_13+", "IFNg+",
                                               "TNFa+", "IL17a+", "CD154+", 
                                               "CD107a+"))
```
Construct the COMPASSContainer. If the number of parent cells is less than countFilterThreshold, we drop that file (default is 5000 cells).
```R
CC <- COMPASSContainerFromGatingSet(gs,
                                    node = parent_node,
                                    individual_id = id,
                                    mp = markermap,
                                    countFilterThreshold = 5000)
```
Look at some basic info about our COMPASSContainer.
```R
CC
```

## Run COMPASS
Fit the COMPASS model using the COMPASSContainer. To fit the COMPASS model, we need to specify how to identify the samples that are our treatment condition and our control condition. Here, we will run COMPASS on the samples stimmed by spike 1 with DMSO as our negative control. For now, let's just do 100 iterations for speed.
```R
fit <- COMPASS(CC,
               treatment = Stim == "S1",
               control = Stim == "DMSO",
               iterations = 100)
```
Save the COMPASS run output.
```R
if(!dir.exists(here::here("out/COMPASSResult"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/COMPASSResult")))
  dir.create(here::here("out/COMPASSResult"), recursive = T)
}
saveRDS(fit, file.path("out/COMPASSResult", "COMPASSResult.rds"))
```
Save the Functionality and Polyfunctionality Scores.
```R
FS <- FunctionalityScore(fit)
FS_df <- data.frame(tmp = names(FS), FS = FS)
colnames(FS_df) <- c("SAMPLE ID", "FS")

PFS <- PolyfunctionalityScore(fit)
PFS_df <- data.frame(tmp = names(PFS), PFS = PFS)
colnames(PFS_df) <- c("SAMPLE ID", "PFS")

FS_PFS_df <- merge(FS_df, PFS_df, by = "SAMPLE ID")
write.table(FS_PFS_df,
            file = file.path("out/COMPASSResult", "FS_PFS.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
FS_PFS_df
```

## Visualize
Plot a heatmap of the mean probability of response.
```R
plot(fit, show_rownames = TRUE)
```
