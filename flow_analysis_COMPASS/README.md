Flow Analysis COMPASS Tutorial
================
Jolie Phan, Emma Bishop
version December 04, 2024

# Directory structure

To directly use the code in this tutorial, you should create an RStudio
project and set up your project directory as follows.

The `data/` folder should contain all `.xml` and `.fcs` files and the
`out/` folder will contain the GatingSet and all COMPASS outputs.
![image](https://private-user-images.githubusercontent.com/46635347/392623484-ddd02c3a-6b4b-41a5-8fa1-ad0d6801efbe.png)

You can achieve this directory structure by running the following
commands:

``` r
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
```

    ## Creating folder /home/emmabishop/workspace/tutorials/out/GatingSet

``` r
if(!dir.exists(here::here("out/COMPASSResult"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/COMPASSResult")))
  dir.create(here::here("out/COMPASSResult"), recursive = T)
}
```

    ## Creating folder /home/emmabishop/workspace/tutorials/out/COMPASSResult

# Load libraries

``` r
library(here)
library(CytoML)
library(flowCore)
library(flowWorkspace)
library(COMPASS)
library(tidyverse)
```

# Load data

Download the .xml file and folder containing the associated .fcs files
(see the README for updated directions). Then, drag them into the “data”
folder of the project directory.

This data has been gated in FlowJo v9 and will be parsed using
flowWorkspace.

**FYI:** When using FlowJo v9, the FlowJo workspace must be exported as
an .xml file to create a flowjo_workspace object with the function
open_flowjo_xml(). However, when using FlowJo v10, the FlowJo workspace
can be loaded directly as a .wsp file using the same function
open_flow_xml().

``` r
# Location of XML file
xml_path <- here::here("data/20200605_COVID_ICS-B3-trunc.xml")

# Location of .fcs files
fcs_subfolder <- here::here("data/20200605_COVID_ICS-B3-FCS-trunc/")
```

Create a flowjo_workspace object with the function open_flowjo_xml().

``` r
ws <- open_flowjo_xml(xml_path)
```

# Create a GatingSet

## Set-up

A GatingSet holds a set of GatingHierarchy objects, representing a set
of samples and the gating scheme associated with each.

Look at the workspace metadata to choose which keywords to extract into
the GatingSet. The flowjo_to_gatingset() function parses a flowJo
Workspace to generate a GatingSet object.

``` r
# Look at all of the keywords
names(fj_ws_get_keywords(ws, 117)) 
```

    ##   [1] "$BEGINANALYSIS"               "$BEGINDATA"                   "$BEGINSTEXT"                  "$BTIM"                        "$BYTEORD"                     "$DATATYPE"                   
    ##   [7] "$DATE"                        "$ENDANALYSIS"                 "$ENDDATA"                     "$ENDSTEXT"                    "$ETIM"                        "$FIL"                        
    ##  [13] "$INST"                        "$MODE"                        "$NEXTDATA"                    "$PAR"                         "$P1B"                         "$P1E"                        
    ##  [19] "$P1G"                         "$P1N"                         "$P1R"                         "$P2B"                         "$P2E"                         "$P2G"                        
    ##  [25] "$P2N"                         "$P2R"                         "$P2V"                         "$P3B"                         "$P3E"                         "$P3G"                        
    ##  [31] "$P3N"                         "$P3R"                         "$P3V"                         "$P4B"                         "$P4E"                         "$P4G"                        
    ##  [37] "$P4N"                         "$P4R"                         "$P4V"                         "$P5B"                         "$P5E"                         "$P5G"                        
    ##  [43] "$P5N"                         "$P5R"                         "$P5V"                         "$P6B"                         "$P6E"                         "$P6G"                        
    ##  [49] "$P6N"                         "$P6R"                         "$P6S"                         "$P6V"                         "$P7B"                         "$P7E"                        
    ##  [55] "$P7G"                         "$P7N"                         "$P7R"                         "$P7S"                         "$P7V"                         "$P8B"                        
    ##  [61] "$P8E"                         "$P8G"                         "$P8N"                         "$P8R"                         "$P8S"                         "$P8V"                        
    ##  [67] "$P9B"                         "$P9E"                         "$P9G"                         "$P9N"                         "$P9R"                         "$P9S"                        
    ##  [73] "$P9V"                         "$P10B"                        "$P10E"                        "$P10G"                        "$P10N"                        "$P10R"                       
    ##  [79] "$P10S"                        "$P10V"                        "$P11B"                        "$P11E"                        "$P11G"                        "$P11N"                       
    ##  [85] "$P11R"                        "$P11S"                        "$P11V"                        "$P12B"                        "$P12E"                        "$P12G"                       
    ##  [91] "$P12N"                        "$P12R"                        "$P12S"                        "$P12V"                        "$P13B"                        "$P13E"                       
    ##  [97] "$P13G"                        "$P13N"                        "$P13R"                        "$P13S"                        "$P13V"                        "$P14B"                       
    ## [103] "$P14E"                        "$P14G"                        "$P14N"                        "$P14R"                        "$P14S"                        "$P14V"                       
    ## [109] "$P15B"                        "$P15E"                        "$P15G"                        "$P15N"                        "$P15R"                        "$P15S"                       
    ## [115] "$P15V"                        "$P16B"                        "$P16E"                        "$P16G"                        "$P16N"                        "$P16R"                       
    ## [121] "$P16S"                        "$P16V"                        "$P17B"                        "$P17E"                        "$P17G"                        "$P17N"                       
    ## [127] "$P17R"                        "$P17S"                        "$P17V"                        "$P18B"                        "$P18E"                        "$P18G"                       
    ## [133] "$P18N"                        "$P18R"                        "$P18S"                        "$P18V"                        "$P19B"                        "$P19E"                       
    ## [139] "$P19G"                        "$P19N"                        "$P19R"                        "$P19S"                        "$P19V"                        "$P20B"                       
    ## [145] "$P20E"                        "$P20G"                        "$P20N"                        "$P20R"                        "$P20S"                        "$P20V"                       
    ## [151] "$P21B"                        "$P21E"                        "$P21G"                        "$P21N"                        "$P21R"                        "$P21S"                       
    ## [157] "$P21V"                        "$SRC"                         "$SYS"                         "$TIMESTEP"                    "$TOT"                         "AUTOBS"                      
    ## [163] "CREATOR"                      "CST BASELINE DATE"            "CST BEADS EXPIRED"            "CST BEADS LOT ID"             "CST PERFORMANCE EXPIRED"      "CST REGULATORY STATUS"       
    ## [169] "CST SETUP DATE"               "CST SETUP STATUS"             "CYTOMETER CONFIG CREATE DATE" "CYTOMETER CONFIG NAME"        "EXPERIMENT NAME"              "EXPORT TIME"                 
    ## [175] "EXPORT USER NAME"             "FJ_$TIMESTEP"                 "FSC ASF"                      "GUID"                         "PLATE ID"                     "PLATE NAME"                  
    ## [181] "P1BS"                         "P1MS"                         "P2BS"                         "P2DISPLAY"                    "P2MS"                         "P3BS"                        
    ## [187] "P3DISPLAY"                    "P3MS"                         "P4BS"                         "P4DISPLAY"                    "P4MS"                         "P5BS"                        
    ## [193] "P5DISPLAY"                    "P5MS"                         "P6BS"                         "P6DISPLAY"                    "P6MS"                         "P7BS"                        
    ## [199] "P7DISPLAY"                    "P7MS"                         "P8BS"                         "P8DISPLAY"                    "P8MS"                         "P9BS"                        
    ## [205] "P9DISPLAY"                    "P9MS"                         "P10BS"                        "P10DISPLAY"                   "P10MS"                        "P11BS"                       
    ## [211] "P11DISPLAY"                   "P11MS"                        "P12BS"                        "P12DISPLAY"                   "P12MS"                        "P13BS"                       
    ## [217] "P13DISPLAY"                   "P13MS"                        "P14BS"                        "P14DISPLAY"                   "P14MS"                        "P15BS"                       
    ## [223] "P15DISPLAY"                   "P15MS"                        "P16BS"                        "P16DISPLAY"                   "P16MS"                        "P17BS"                       
    ## [229] "P17DISPLAY"                   "P17MS"                        "P18BS"                        "P18DISPLAY"                   "P18MS"                        "P19BS"                       
    ## [235] "P19DISPLAY"                   "P19MS"                        "P20BS"                        "P20DISPLAY"                   "P20MS"                        "P21BS"                       
    ## [241] "P21DISPLAY"                   "P21MS"                        "SETTINGS"                     "SPILL"                        "STIM"                         "THRESHOLD"                   
    ## [247] "TUBE NAME"                    "WELL ID"                      "WINDOW EXTENSION"

``` r
# Choose which keywords to keep
keywords2import <- c("EXPERIMENT NAME",
                       "$DATE",
                       "SAMPLE ID",
                       "PATIENT ID",
                       "STIM",
                       "WELL ID",
                       "PLATE NAME") 
sampleGroup <- "Samples"

gs <- flowjo_to_gatingset(ws,                                    
                          name = sampleGroup, 
                          keywords = keywords2import,
                          path = fcs_subfolder, 
                          extend_val = -10000)
```

## QC

Make sure that the gating trees are consistent for all samples.

``` r
pop_lists <- lapply(gs, gh_get_pop_paths)
unique(pop_lists)
```

    ## [[1]]
    ##  [1] "root"                                      "/Time"                                     "/Time/LD-3+"                               "/Time/LD-3+/1419-3+"                      
    ##  [5] "/Time/LD-3+/1419-3+/S"                     "/Time/LD-3+/1419-3+/S/Lymph"               "/Time/LD-3+/1419-3+/S/Lymph/4+"            "/Time/LD-3+/1419-3+/S/Lymph/4+/107a"      
    ##  [9] "/Time/LD-3+/1419-3+/S/Lymph/4+/154"        "/Time/LD-3+/1419-3+/S/Lymph/4+/CCR7+"      "/Time/LD-3+/1419-3+/S/Lymph/4+/CD45RA+"    "/Time/LD-3+/1419-3+/S/Lymph/4+/IFNG"      
    ## [13] "/Time/LD-3+/1419-3+/S/Lymph/4+/IL2"        "/Time/LD-3+/1419-3+/S/Lymph/4+/IL17"       "/Time/LD-3+/1419-3+/S/Lymph/4+/IL4513"     "/Time/LD-3+/1419-3+/S/Lymph/4+/TNF"       
    ## [17] "/Time/LD-3+/1419-3+/S/Lymph/8+"            "/Time/LD-3+/1419-3+/S/Lymph/8+/IFNG"       "/Time/LD-3+/1419-3+/S/Lymph/CD38+"         "/Time/LD-3+/1419-3+/S/Lymph/HLADR+"       
    ## [21] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+"         "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/107a"    "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/154"     "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/CCR7+"  
    ## [25] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/CD45RA+" "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IFNG"    "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL2"     "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL17"   
    ## [29] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL4513"  "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/TNF"

Remove channels from flow data that are not used by gates.

``` r
gs <- gs_remove_redundant_channels(gs) # drop SSC-H
```

    ## drop SSC-H

Add names to all channels or change their names.

``` r
dput(unname(pData(parameters(gh_pop_get_data(gs[[1]])))[,2]))
```

    ## structure(c(NA, NA, NA, NA, "CD8b BB700", "TNFa FITC", "CD107a PE-Cy7", 
    ## "CD154 PE-Cy5", "CD3 ECD", "IL2 PE", "CD4 APC-H7", "IL17a Ax700", 
    ## "IL4/5/13 APC", "CD14/CD19 BV785", "CCR7 BV711", "CD38 BV605", 
    ## "L/D", "IFNg V450", "CD45RA BUV737", "HLADR BUV395"), class = "AsIs")

``` r
markernames <- c("Time", "FSC-A", "FSC-H", "SSC-A", "CD8b", "TNFa", "CD107a",
                 "CD154", "CD3", "IL2", "CD4", "IL17a", "IL4_5_13", "CD14_19",
                 "CCR7", "CD38", "LD", "IFNg", "CD45RA", "HLADR")
names(markernames) <- pData(parameters(gh_pop_get_data(gs[[1]])))[,1]
markernames(gs) <- markernames
pData(parameters(gh_pop_get_data(gs[[1]])))[,c(1,2)]
```

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
    ## $P17 <V610-A>     CD38
    ## $P18 <V510-A>       LD
    ## $P19 <V450-A>     IFNg
    ## $P20 <U730-A>   CD45RA
    ## $P21 <U395-A>    HLADR

In this tutorial, we will only run COMPASS on CD4+ T cells, but what if
we wanted to also run COMPASS on CD8+ T cells? Currently, most of the
gates are missing from the “8+” node. Let’s grab the missing gates from
the “NOT4+” node and add them under the “8+” node.

``` r
gates_to_copy <- c("/Time/LD-3+/1419-3+/S/Lymph/NOT4+/107a",
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/154", 
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL2",
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL17",
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL4513",
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/TNF",
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/CCR7+",
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/CD45RA+")

for(path in gates_to_copy) {
  gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y=path),
             parent = "/Time/LD-3+/1419-3+/S/Lymph/8+")
}
recompute(gs, "/Time/LD-3+/1419-3+/S/Lymph/8+")
```

    ## done!

Next, let’s get rid of the “NOT4+” node and all of its descendants.

``` r
gs_pop_remove(gs, "/Time/LD-3+/1419-3+/S/Lymph/NOT4+")
```

Plot the gating tree.

``` r
plot(gs, fontsize=15, bool=T)
```

![](/home/emmabishop/workspace/tutorials/flow_analysis_COMPASS/README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

## Save GatingSet

**Note:** this can take a while.

``` r
save_gs(gs, here::here("out/GatingSet"))
```

    ## Done
    ## To reload it, use 'load_gs' function

# Get tabular data from the GatingSet

The following steps were inspired from advice in the COMPASS [issue
\#81](https://github.com/RGLab/COMPASS/issues/81#issuecomment-2484533481)
thread.

``` r
# Load GatingSet if needed
gs <- load_gs(here::here("out/GatingSet"))
```

## Get matrices of expression data

Extract metadata and take a look

``` r
metadata <- pData(gs) %>%
  mutate(sampleName = rownames(.))

head(metadata)
```

    ##                         name       EXPERIMENT NAME       $DATE SAMPLE ID PATIENT ID    STIM WELL ID PLATE NAME        sampleName
    ## 114716.fcs_445737 114716.fcs 20200605_COVID_ICS-B3 05-JUN-2020        23         23    DMSO     G02         P3 114716.fcs_445737
    ## 114718.fcs_426995 114718.fcs 20200605_COVID_ICS-B3 05-JUN-2020        23         23     SEB     G04         P3 114718.fcs_426995
    ## 114720.fcs_339914 114720.fcs 20200605_COVID_ICS-B3 05-JUN-2020        23         23    VEMP     G06         P3 114720.fcs_339914
    ## 114722.fcs_383957 114722.fcs 20200605_COVID_ICS-B3 05-JUN-2020        23         23 Spike 1     G08         P3 114722.fcs_383957
    ## 114724.fcs_273422 114724.fcs 20200605_COVID_ICS-B3 05-JUN-2020        23         23 Spike 2     G10         P3 114724.fcs_273422
    ## 114726.fcs_373199 114726.fcs 20200605_COVID_ICS-B3 05-JUN-2020        23         23    NCAP     G12         P3 114726.fcs_373199

Separate the Spike 1 stimmed and DMSO samples

``` r
sp1_samples <- metadata %>%
  filter(STIM == "Spike 1") %>%
  pull(sampleName)
spike1 <- gs[sp1_samples]

dmso_samples <- metadata %>%
  filter(STIM == "DMSO") %>%
  pull(sampleName)
dmso <- gs[dmso_samples]
```

Extract data for the “CD4+” population only. These are ‘cytoset’
objects, which are a collection of ’cytoframe’s.

``` r
spike1_cd4 <- gs_pop_get_data(spike1, "4+")  # Using "+4" node name
dmso_cd4 <- gs_pop_get_data(dmso, "4+")
```

We can see the sample names in the object like so:

``` r
sampleNames(dmso_cd4)
```

    ## [1] "114716.fcs_445737" "114776.fcs_846625" "114752.fcs_489936" "114788.fcs_371622" "114764.fcs_435350"

Let’s take a look at one of these ‘cytoframes’, which are similar to a
matrix of expression values.

``` r
samplename <- "114716.fcs_445737"
cytoframe1 <- dmso_cd4[[samplename]]
head(cytoframe1)
```

    ##       Time     FSC-A  FSC-H    SSC-A <B710-A> <B515-A>  <G780-A> <G660-A> <G610-A>  <G575-A> <R780-A>  <R710-A> <R660-A>  <V780-A> <V710-A> <V610-A> <V510-A> <V450-A>  <U730-A>  <U395-A>
    ## [1,] 2.441 122933.04 111005 22292.88 912.8039 605.4547 220.27991 560.6566 3060.245  395.8604 1715.682  577.2411 731.5930  654.3718 2341.164 1162.096 1066.730 757.4949 1392.2936  691.1917
    ## [2,] 2.441 126327.96 110188 34016.13 597.5212 732.2466 190.20494 873.4103 2982.977  410.8742 1868.221 1142.5469 436.4252 1132.8197 1025.377 1015.737 1221.067 837.0206  779.2865  578.2829
    ## [3,] 2.442  93704.20  80479 30572.67 897.9616 617.5473  78.86456 440.2558 3033.829  995.1781 2014.507  783.3704 659.7129 1403.7660 2366.746 1613.703 1262.292 405.5316 2563.0171  568.5731
    ## [4,] 2.447 109862.56  95799 24374.79 316.1867 428.8226 893.86353 997.1701 3053.399  976.4500 1855.618  381.2457 654.4659  889.6058 2180.418 2253.460 1215.396 645.0524 2731.7053 1013.0549
    ## [5,] 2.449  82438.72  70628 40428.90 766.8693 484.2879 911.77679 915.0909 3070.636  630.3830 2046.922 1543.0740 401.9236  890.8912 2353.598 1283.655 1216.928 919.3452 2743.6641  768.1965
    ## [6,] 2.453  80187.60  68863 36685.29 469.5395 537.9316 433.60690 789.3820 3119.851 1252.5963 2124.699 1185.0879 831.7593  580.8347 1753.941 1515.809 1387.849 771.2418  884.6489 1112.8729

Note that we have gates we don’t need for our boolean combinations, like
time and forward/side scatter (Time, FSC-A, FSC-H, SSC-A). We will get
rid of these below.

Get a list of expression matrices for stim and unstim from the
cytoset/cytoframe objects.

``` r
markernames <- c("Time", "FSC-A", "FSC-H", "SSC-A", "CD8b", "TNFa", "CD107a",
                 "CD154", "CD3", "IL2", "CD4", "IL17a", "IL4_5_13", "CD14_19",
                 "CCR7", "CD38", "LD", "IFNg", "CD45RA", "HLADR")

markermap <- list("IL2", "IL4_5_13", "IFNg", "TNFa", "IL17a", "CD154", "CD107a")

get_mtx_list <- function(cytoset_obj) {
  outlist = list()
  for (i in sampleNames(cytoset_obj)) {
    sub_cf <- get_cytoframe_from_cs(cytoset_obj, i)  # Extract cytoframe from cytoset
    out <- exprs(sub_cf)  # Get expression values as a matrix
    colnames(out) <- markernames  # Rename columns (e.g. <V450-A> becomes IFNg)
    out[out < 0] <- 0  # Replace negative numbers with 0 for CellCounts() later
    out <- out[, colnames(out) %in% markermap]  # Only keep markers we want boolean combos of
    outlist[[i]] <- out
  }
  return(outlist)
}

spikelist <- get_mtx_list(spike1_cd4)
dmsolist <- get_mtx_list(dmso_cd4)
```

Now we have named lists of matrices of expression data, using sample
names as indices

``` r
names(dmsolist)
```

    ## [1] "114716.fcs_445737" "114776.fcs_846625" "114752.fcs_489936" "114788.fcs_371622" "114764.fcs_435350"

``` r
head(dmsolist[[1]])
```

    ##          TNFa    CD107a    CD154       IL2     IL17a IL4_5_13     IFNg
    ## [1,] 605.4547 220.27991 560.6566  395.8604  577.2411 731.5930 757.4949
    ## [2,] 732.2466 190.20494 873.4103  410.8742 1142.5469 436.4252 837.0206
    ## [3,] 617.5473  78.86456 440.2558  995.1781  783.3704 659.7129 405.5316
    ## [4,] 428.8226 893.86353 997.1701  976.4500  381.2457 654.4659 645.0524
    ## [5,] 484.2879 911.77679 915.0909  630.3830 1543.0740 401.9236 919.3452
    ## [6,] 537.9316 433.60690 789.3820 1252.5963 1185.0879 831.7593 771.2418

Save

``` r
saveRDS(spikelist, here::here("out/spikelist.rds"))
saveRDS(dmsolist, here::here("out/dmsolist.rds"))
```

## Create boolean combinations

``` r
spikelist <- readRDS(here::here("out/spikelist.rds"))
dmsolist <- readRDS(here::here("out/dmsolist.rds"))
```

Can take a little while

``` r
# Seven markers
n_s <- CellCounts(spikelist, Combinations(7)) # Stimmed data
n_u <- CellCounts(dmsolist, Combinations(7)) # Unstimmed data
```

Save count matrices and metadata

``` r
saveRDS(n_s, here::here("out/n_s.rds"))
saveRDS(n_u, here::here("out/n_u.rds"))
saveRDS(metadata, here::here("out/metadata.rds"))
```

# Run SimpleCOMPASS

``` r
n_s <- readRDS(here::here("out/n_s.rds"))
n_u <- readRDS(here::here("out/n_u.rds"))
metadata <- readRDS(here::here("out/metadata.rds"))
```

``` r
# Use subject ID as rownames
rownames(n_s) <- c("23", "25", "133C", "142C", "150C")
rownames(n_u) <- c("23", "25", "133C", "142C", "150C")
```

Fits a polyfunctionality response model to our data.

``` r
fit = COMPASS::SimpleCOMPASS(n_s = n_s,
                             n_u = n_u,
                             meta = metadata,
                             individual_id = "SAMPLE ID")
```

    ## Setting the seed to  100

Save the COMPASS run output.

``` r
saveRDS(fit, file.path(here::here("out/COMPASSResult"), "COMPASSResult.rds"))
# fit <- readRDS(file.path(here::here("out/COMPASSResult"), "COMPASSResult.rds"))
```

Save the Functionality and Polyfunctionality Scores.

``` r
FS <- FunctionalityScore(fit)
FS_df <- data.frame(tmp = names(FS), FS = FS)
colnames(FS_df) <- c("sampleName", "FS")

PFS <- PolyfunctionalityScore(fit)
PFS_df <- data.frame(tmp = names(PFS), PFS = PFS)
colnames(PFS_df) <- c("sampleName", "PFS")

FS_PFS_df <- merge(FS_df, PFS_df, by = "sampleName")
write.table(FS_PFS_df,
            file = file.path(here::here("out/COMPASSResult"), "FS_PFS.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
FS_PFS_df
```

    ##   sampleName          FS        PFS
    ## 1       133C 0.010739370 0.04174320
    ## 2       142C 0.019429134 0.07537653
    ## 3       150C 0.014395276 0.05564915
    ## 4         23 0.007925984 0.03066837
    ## 5         25 0.031719685 0.12269048

# Visualize

Plot a heatmap of the mean probability of response.

``` r
plot(fit, show_rownames = TRUE)
```

![](/home/emmabishop/workspace/tutorials/flow_analysis_COMPASS/README_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->
