# Flow_Analysis_COMPASS

A tutorial for analyzing flow cytometry data in R using COMPASS (Combinatorial Polyfunctionality Analysis of Single Cells). 

# Introduction

This tutorial contains a basic workflow for analyzing flow cytometry data using COMPASS. COMPASS is a statistical framework that allows for unbiased analysis of antigen-specific T-cell subsets. COMPASS uses a Bayesian hierarchial framework to model all obeserved cell-subsets and select the most likely to be antigen-specific while regularizing the small cell counts that often arise in multi-parameter space. The model gives a posterior probability of specificity for each cell-subset and each sample, which can be used to profile a subject's immune response to external stimuli (e.g., infection or vaccination). 

For this tutorial, we will be analyzing a batch of data from the HAARVI cohort. For more information about COMPASS, please refer to the original manuscript and documentation:
* https://www.nature.com/articles/nbt.3187
* https://bioconductor.org/packages/COMPASS/

# Installation

This pipeline should be completed in R and RStudio. The following R packages are required for this tutorial:
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
# Directory Structure
To directly use the code in this tutorial, you should create an RStudio project and set up your project directory as follows. The `data/` folder should contain all `.xml` and `.fcs` files and the `out/` folder will contain the GatingSet all COMPASS outputs.
![image](https://user-images.githubusercontent.com/89667908/147301852-f5c1d505-cb04-4841-bdbe-981b0d4bc6f9.png)

You can achieve this directory structure by running the following command lines:
```R
if(!dir.exists(here::here("data"))) {
  cat(sprintf("Creating folder %s\n", here::here("data")))
  dir.create(here::here("data"), recursive = T)
}
if(!dir.exists(here::here("out"))) {
  cat(sprintf("Creating folder %s\n", here::here("out")))
  dir.create(here::here("out"), recursive = T)
}
if(!dir.exists(here::here("out/GatingSet"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/GatingSet")))
  dir.create(here::here("out/GatingSet"), recursive = T)
}
if(!dir.exists(here::here("out/COMPASSResult"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/COMPASSResult")))
  dir.create(here::here("out/COMPASSResult"), recursive = T)
}
```

# Workflow Overview
1. Load data
2. Create a GatingSet
3. Create a COMPASSContainer
4. Run COMPASS
5. Visualize

## Load data
```R
# Load libraries into your current R session
library(here)
library(CytoML)
library(flowCore)
library(flowWorkspace)
library(COMPASS)
```
The .xml and .fcs files for this dataset are stored in the Seshadri Lab shared drive (LSR Fortessa/2021 Summer HAARVIVAC/20210902 HAARVIVAC Batch 4/).
Download the .xml file and folder containing the associated .fcs files and drag them into the "data" folder of the project directory.

This data has been gated in FlowJo v9 and will be parsed using flowWorkspace.
***FYI:*** *When using FlowJo v9, the FlowJo workspace must be exported as an .xml file to create a flowjo_workspace object with the function open_flowjo_xml(). However, when using FlowJo v10, the FlowJo workspace can be loaded directly as a .wsp file using the same function open_flow_xml().*
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

### Quality Control (QC)
Make sure that the gating trees are consistent for all samples.
```R
pop_lists <- lapply(gs, gh_get_pop_paths)
unique(pop_lists)
```
```R
## [[1]]
##  [1] "root"                                                         
##  [2] "/Time"                                                        
##  [3] "/Time/CD3+"                                                   
##  [4] "/Time/CD3+/CD14-CD19-"                                        
##  [5] "/Time/CD3+/CD14-CD19-/Lymphocytes"                            
##  [6] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet"                    
##  [7] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live"               
##  [8] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+"          
##  [9] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CCR7+"    
## [10] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CD45RA+"  
## [11] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CD107a+"  
## [12] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/CD154+"   
## [13] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/IFNg+"    
## [14] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/IL2+"     
## [15] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/IL4_5_13+"
## [16] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/IL17a+"   
## [17] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD4+/TNFa+"    
## [18] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+"          
## [19] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CCR7+"    
## [20] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CD45RA+"  
## [21] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CD107a+"  
## [22] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/CD154+"   
## [23] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/IFNg+"    
## [24] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/IL2+"     
## [25] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/IL4_5_13+"
## [26] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/IL17a+"   
## [27] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD8+/TNFa+"    
## [28] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/CD38+"         
## [29] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/DN"            
## [30] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/DP"            
## [31] "/Time/CD3+/CD14-CD19-/Lymphocytes/Singlet/Live/HLADR+"
```
Remove channels from flow data that are not used by gates.
```R
gs <- gs_remove_redundant_channels(gs) 
```
```R
## drop SSC-H, V655-A, V570-A
```
Add names to all channels or fix their names.
```R
dput(unname(pData(parameters(gh_pop_get_data(gs[[1]])))[,2]))
```
```R
## structure(c(NA, NA, NA, NA, "CD8b", "TNFa", "CD107a", "CD154", 
## "CD3", "IL2", "CD4", "IL17a", "IL4_5_13", "CD14_19", "CCR7", 
## "CD38", "LD", "INFg", "CD45RA", "HLADR"), class = "AsIs")
```
```R
markernames <- c("Time", "FSC-A", "FSC-H", "SSC-A", "CD8b", "TNFa", "CD107a",
                 "CD154", "CD3", "IL2", "CD4", "IL17a", "IL4_5_13", "CD14_19",
                 "CCR7", "CD38", "LD", "IFNg", "CD45RA", "HLADR")
names(markernames) <- pData(parameters(gh_pop_get_data(gs[[1]])))[,1]
markernames(gs) <- markernames
pData(parameters(gh_pop_get_data(gs[[1]])))[,c(1,2)]
```
```R
##          name     desc
## $P1      Time     Time
## $P2     FSC-A    FSC-A
## $P3     FSC-H    FSC-H
## $P4     SSC-A    SSC-A
## $P6  <B710-A>     CD8b
## $P7  <B515-A>     TNFa
## $P8  <G780-A>   CD107a
## $P9  <G660-A>    CD154
## $P10 <G610-A>      CD3
## $P11 <G575-A>      IL2
## $P12 <R780-A>      CD4
## $P13 <R710-A>    IL17a
## $P14 <R660-A> IL4_5_13
## $P15 <V780-A>  CD14_19
## $P16 <V710-A>     CCR7
## $P18 <V610-A>     CD38
## $P20 <V510-A>       LD
## $P21 <V450-A>     IFNg
## $P22 <U730-A>   CD45RA
## $P23 <U395-A>    HLADR
```
Plot gating tree.
```R
plot(gs, fontsize=15, bool=T)
```
![image](https://user-images.githubusercontent.com/89667908/145899923-a9269520-6f99-4e70-a170-6a2094fd7554.png)


### Save GatingSet
**Note:** this can take a while.
```R
save_gs(gs, here::here("out/GatingSet"), overwrite = TRUE)
```

## Create a COMPASSContainer
A COMPASSContainer is the data structure used to hold data from an ICS experiment. The input for this code is a GatingSet or a GatingSetList. Counts, metadata, and 
single cell data are extracted and fed into the COMPASSContainer constructor.
```R
# Set the seed
set.seed(123)

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
```R
## A COMPASSContainer with 80 samples from 16 individuals, containing data across 7 markers.
```

## Run COMPASS
Fit the COMPASS model using the COMPASSContainer. To fit the COMPASS model, we need to specify how to identify the samples that are our treatment condition and our control condition based on the metadata. Here, we will run COMPASS on the samples stimmed by spike 1 with DMSO as our negative control. For now, let's just do 100 iterations for speed.
```R
fit <- COMPASS(CC,
               treatment = Stim == "S1",
               control = Stim == "DMSO",
               iterations = 100)
```
Save the COMPASS run output.
```R
saveRDS(fit, file.path(here::here("out/COMPASSResult"), "COMPASSResult.rds"))
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
            file = file.path(here::here("out/COMPASSResult"), "FS_PFS.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
FS_PFS_df
```
```R
##    SAMPLE ID          FS          PFS
## 1     188C-1 0.010866142 0.0083877551
## 2     188C-2 0.057244094 0.0307482993
## 3     191C-1 0.025039370 0.0169455782
## 4     191C-2 0.033385827 0.0195102041
## 5     210C-1 0.005511811 0.0031462585
## 6     210C-2 0.047559055 0.0253061224
## 7     235C-1 0.001181102 0.0006122449
## 8     235C-2 0.020236220 0.0089965986
## 9      27H-1 0.004960630 0.0025204082
## 10     27H-2 0.049921260 0.0269795918
## 11     35H-1 0.013937008 0.0108095238
## 12     35H-2 0.019763780 0.0104387755
## 13     38H-1 0.008740157 0.0075034014
## 14     38H-2 0.046299213 0.0247346939
## 15     44H-1 0.040944882 0.0224319728
## 16     44H-2 0.037874016 0.0225510204
```

## Visualize
Plot a heatmap of the mean probability of response.
```R
plot(fit, show_rownames = TRUE)
```
```R
## The 'threshold' filter has removed 8 categories:
## IL2&!IL4_5_13&!IFNg&!TNFa&!IL17a&!CD154&!CD107a, !IL2&IL4_5_13&!IFNg&!TNFa&!IL17a&!CD154&!CD107a, !IL2&!IL4_5_13&IFNg&!TNFa&!IL17a&!CD154&!CD107a, !IL2&!IL4_5_13&!IFNg&TNFa&!IL17a&!CD154&!CD107a, !IL2&!IL4_5_13&!IFNg&!TNFa&IL17a&!CD154&!CD107a, !IL2&!IL4_5_13&!IFNg&!TNFa&!IL17a&!CD154&CD107a, IL2&!IL4_5_13&!IFNg&TNFa&!IL17a&!CD154&!CD107a, IL2&!IL4_5_13&!IFNg&!TNFa&!IL17a&CD154&!CD107a
```
![image](https://user-images.githubusercontent.com/89667908/145901258-89ba874c-77db-4c3e-b46c-565a9ac0be53.png)

