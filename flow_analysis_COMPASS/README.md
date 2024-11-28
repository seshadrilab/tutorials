# Directory structure

To directly use the code in this tutorial, you should create an RStudio
project and set up your project directory as follows.

The `data/` folder should contain all `.xml` and `.fcs` files and the
`out/` folder will contain the GatingSet and all COMPASS outputs.

<div class="grViz html-widget html-fill-item" id="htmlwidget-6b9f3c6e769063ea8822" style="width:672px;height:480px;"></div>
<script type="application/json" data-for="htmlwidget-6b9f3c6e769063ea8822">{"x":{"diagram":"digraph flowchart {\n  node [fontname = arial, shape = retangle, fontsize=8, height=0.1]\n  dir1 [label = \"Rproject/\"]\n  dir2 [label = \"data/\"]\n  dir3 [label = \"out/\"]\n  dir4 [label = \".xml\"]\n  dir5 [label = \".fcs files\"]\n  dir6 [label = \"GatingSet/\"]\n  dir7 [label = \"COMPASSResult/\"]\n  dir1 -> dir2,dir3;\n  dir2 -> dir4,dir5;\n  dir3 -> dir6,dir7;\n  }","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}</script>

You can achieve this directory structure by running the following
command lines:

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

# Load libraries

    library(here)

    ## here() starts at /home/emmabishop/workspace/tutorials

    library(CytoML)
    library(flowCore)
    library(flowWorkspace)

    ## As part of improvements to flowWorkspace, some behavior of
    ## GatingSet objects has changed. For details, please read the section
    ## titled "The cytoframe and cytoset classes" in the package vignette:
    ## 
    ##   vignette("flowWorkspace-Introduction", "flowWorkspace")

    library(COMPASS)
    library(tidyverse)

    ## ── Attaching core tidyverse packages ────────────────────────────────────────────────────────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2

    ## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks flowCore::filter(), stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

# Load data

Download the .xml file and folder containing the associated .fcs files
(see the README for updated directions). Then, drag them into the “data”
folder of the project directory.

This data has been gated in FlowJo v9 and will be parsed using
flowWorkspace.

**FYI:** *When using FlowJo v9, the FlowJo workspace must be exported as
an .xml file to create a flowjo\_workspace object with the function
open\_flowjo\_xml(). However, when using FlowJo v10, the FlowJo
workspace can be loaded directly as a .wsp file using the same function
open\_flow\_xml().*

    # Location of XML file
    xml_path <- here::here("data/20200605_COVID_ICS-B3-trunc.xml")

    # Location of .fcs files
    fcs_subfolder <- here::here("data/20200605_COVID_ICS-B3-FCS-trunc/")

Create a flowjo\_workspace object with the function open\_flowjo\_xml().

    ws <- open_flowjo_xml(xml_path)

# Create a GatingSet

## Set-up

A GatingSet holds a set of GatingHierarchy objects, representing a set
of samples and the gating scheme associated with each.

Look at the workspace metadata to choose which keywords to extract into
the GatingSet. The flowjo\_to\_gatingset() function parses a flowJo
Workspace to generate a GatingSet object.

    # Look at all of the keywords
    names(fj_ws_get_keywords(ws, 117)) 

    ##   [1] "$BEGINANALYSIS"               "$BEGINDATA"                   "$BEGINSTEXT"                  "$BTIM"                       
    ##   [5] "$BYTEORD"                     "$DATATYPE"                    "$DATE"                        "$ENDANALYSIS"                
    ##   [9] "$ENDDATA"                     "$ENDSTEXT"                    "$ETIM"                        "$FIL"                        
    ##  [13] "$INST"                        "$MODE"                        "$NEXTDATA"                    "$PAR"                        
    ##  [17] "$P1B"                         "$P1E"                         "$P1G"                         "$P1N"                        
    ##  [21] "$P1R"                         "$P2B"                         "$P2E"                         "$P2G"                        
    ##  [25] "$P2N"                         "$P2R"                         "$P2V"                         "$P3B"                        
    ##  [29] "$P3E"                         "$P3G"                         "$P3N"                         "$P3R"                        
    ##  [33] "$P3V"                         "$P4B"                         "$P4E"                         "$P4G"                        
    ##  [37] "$P4N"                         "$P4R"                         "$P4V"                         "$P5B"                        
    ##  [41] "$P5E"                         "$P5G"                         "$P5N"                         "$P5R"                        
    ##  [45] "$P5V"                         "$P6B"                         "$P6E"                         "$P6G"                        
    ##  [49] "$P6N"                         "$P6R"                         "$P6S"                         "$P6V"                        
    ##  [53] "$P7B"                         "$P7E"                         "$P7G"                         "$P7N"                        
    ##  [57] "$P7R"                         "$P7S"                         "$P7V"                         "$P8B"                        
    ##  [61] "$P8E"                         "$P8G"                         "$P8N"                         "$P8R"                        
    ##  [65] "$P8S"                         "$P8V"                         "$P9B"                         "$P9E"                        
    ##  [69] "$P9G"                         "$P9N"                         "$P9R"                         "$P9S"                        
    ##  [73] "$P9V"                         "$P10B"                        "$P10E"                        "$P10G"                       
    ##  [77] "$P10N"                        "$P10R"                        "$P10S"                        "$P10V"                       
    ##  [81] "$P11B"                        "$P11E"                        "$P11G"                        "$P11N"                       
    ##  [85] "$P11R"                        "$P11S"                        "$P11V"                        "$P12B"                       
    ##  [89] "$P12E"                        "$P12G"                        "$P12N"                        "$P12R"                       
    ##  [93] "$P12S"                        "$P12V"                        "$P13B"                        "$P13E"                       
    ##  [97] "$P13G"                        "$P13N"                        "$P13R"                        "$P13S"                       
    ## [101] "$P13V"                        "$P14B"                        "$P14E"                        "$P14G"                       
    ## [105] "$P14N"                        "$P14R"                        "$P14S"                        "$P14V"                       
    ## [109] "$P15B"                        "$P15E"                        "$P15G"                        "$P15N"                       
    ## [113] "$P15R"                        "$P15S"                        "$P15V"                        "$P16B"                       
    ## [117] "$P16E"                        "$P16G"                        "$P16N"                        "$P16R"                       
    ## [121] "$P16S"                        "$P16V"                        "$P17B"                        "$P17E"                       
    ## [125] "$P17G"                        "$P17N"                        "$P17R"                        "$P17S"                       
    ## [129] "$P17V"                        "$P18B"                        "$P18E"                        "$P18G"                       
    ## [133] "$P18N"                        "$P18R"                        "$P18S"                        "$P18V"                       
    ## [137] "$P19B"                        "$P19E"                        "$P19G"                        "$P19N"                       
    ## [141] "$P19R"                        "$P19S"                        "$P19V"                        "$P20B"                       
    ## [145] "$P20E"                        "$P20G"                        "$P20N"                        "$P20R"                       
    ## [149] "$P20S"                        "$P20V"                        "$P21B"                        "$P21E"                       
    ## [153] "$P21G"                        "$P21N"                        "$P21R"                        "$P21S"                       
    ## [157] "$P21V"                        "$SRC"                         "$SYS"                         "$TIMESTEP"                   
    ## [161] "$TOT"                         "AUTOBS"                       "CREATOR"                      "CST BASELINE DATE"           
    ## [165] "CST BEADS EXPIRED"            "CST BEADS LOT ID"             "CST PERFORMANCE EXPIRED"      "CST REGULATORY STATUS"       
    ## [169] "CST SETUP DATE"               "CST SETUP STATUS"             "CYTOMETER CONFIG CREATE DATE" "CYTOMETER CONFIG NAME"       
    ## [173] "EXPERIMENT NAME"              "EXPORT TIME"                  "EXPORT USER NAME"             "FJ_$TIMESTEP"                
    ## [177] "FSC ASF"                      "GUID"                         "PLATE ID"                     "PLATE NAME"                  
    ## [181] "P1BS"                         "P1MS"                         "P2BS"                         "P2DISPLAY"                   
    ## [185] "P2MS"                         "P3BS"                         "P3DISPLAY"                    "P3MS"                        
    ## [189] "P4BS"                         "P4DISPLAY"                    "P4MS"                         "P5BS"                        
    ## [193] "P5DISPLAY"                    "P5MS"                         "P6BS"                         "P6DISPLAY"                   
    ## [197] "P6MS"                         "P7BS"                         "P7DISPLAY"                    "P7MS"                        
    ## [201] "P8BS"                         "P8DISPLAY"                    "P8MS"                         "P9BS"                        
    ## [205] "P9DISPLAY"                    "P9MS"                         "P10BS"                        "P10DISPLAY"                  
    ## [209] "P10MS"                        "P11BS"                        "P11DISPLAY"                   "P11MS"                       
    ## [213] "P12BS"                        "P12DISPLAY"                   "P12MS"                        "P13BS"                       
    ## [217] "P13DISPLAY"                   "P13MS"                        "P14BS"                        "P14DISPLAY"                  
    ## [221] "P14MS"                        "P15BS"                        "P15DISPLAY"                   "P15MS"                       
    ## [225] "P16BS"                        "P16DISPLAY"                   "P16MS"                        "P17BS"                       
    ## [229] "P17DISPLAY"                   "P17MS"                        "P18BS"                        "P18DISPLAY"                  
    ## [233] "P18MS"                        "P19BS"                        "P19DISPLAY"                   "P19MS"                       
    ## [237] "P20BS"                        "P20DISPLAY"                   "P20MS"                        "P21BS"                       
    ## [241] "P21DISPLAY"                   "P21MS"                        "SETTINGS"                     "SPILL"                       
    ## [245] "STIM"                         "THRESHOLD"                    "TUBE NAME"                    "WELL ID"                     
    ## [249] "WINDOW EXTENSION"

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

## QC

Make sure that the gating trees are consistent for all samples.

    pop_lists <- lapply(gs, gh_get_pop_paths)
    unique(pop_lists)

    ## [[1]]
    ##  [1] "root"                                      "/Time"                                    
    ##  [3] "/Time/LD-3+"                               "/Time/LD-3+/1419-3+"                      
    ##  [5] "/Time/LD-3+/1419-3+/S"                     "/Time/LD-3+/1419-3+/S/Lymph"              
    ##  [7] "/Time/LD-3+/1419-3+/S/Lymph/4+"            "/Time/LD-3+/1419-3+/S/Lymph/4+/107a"      
    ##  [9] "/Time/LD-3+/1419-3+/S/Lymph/4+/154"        "/Time/LD-3+/1419-3+/S/Lymph/4+/CCR7+"     
    ## [11] "/Time/LD-3+/1419-3+/S/Lymph/4+/CD45RA+"    "/Time/LD-3+/1419-3+/S/Lymph/4+/IFNG"      
    ## [13] "/Time/LD-3+/1419-3+/S/Lymph/4+/IL2"        "/Time/LD-3+/1419-3+/S/Lymph/4+/IL17"      
    ## [15] "/Time/LD-3+/1419-3+/S/Lymph/4+/IL4513"     "/Time/LD-3+/1419-3+/S/Lymph/4+/TNF"       
    ## [17] "/Time/LD-3+/1419-3+/S/Lymph/8+"            "/Time/LD-3+/1419-3+/S/Lymph/8+/IFNG"      
    ## [19] "/Time/LD-3+/1419-3+/S/Lymph/CD38+"         "/Time/LD-3+/1419-3+/S/Lymph/HLADR+"       
    ## [21] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+"         "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/107a"   
    ## [23] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/154"     "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/CCR7+"  
    ## [25] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/CD45RA+" "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IFNG"   
    ## [27] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL2"     "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL17"   
    ## [29] "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL4513"  "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/TNF"

Remove channels from flow data that are not used by gates.

    gs <- gs_remove_redundant_channels(gs) # drop SSC-H

    ## drop SSC-H

Add names to all channels or change their names.

    dput(unname(pData(parameters(gh_pop_get_data(gs[[1]])))[,2]))

    ## structure(c(NA, NA, NA, NA, "CD8b BB700", "TNFa FITC", "CD107a PE-Cy7", 
    ## "CD154 PE-Cy5", "CD3 ECD", "IL2 PE", "CD4 APC-H7", "IL17a Ax700", 
    ## "IL4/5/13 APC", "CD14/CD19 BV785", "CCR7 BV711", "CD38 BV605", 
    ## "L/D", "IFNg V450", "CD45RA BUV737", "HLADR BUV395"), class = "AsIs")

    markernames <- c("Time", "FSC-A", "FSC-H", "SSC-A", "CD8b", "TNFa", "CD107a",
                     "CD154", "CD3", "IL2", "CD4", "IL17a", "IL4_5_13", "CD14_19",
                     "CCR7", "CD38", "LD", "IFNg", "CD45RA", "HLADR")
    names(markernames) <- pData(parameters(gh_pop_get_data(gs[[1]])))[,1]
    markernames(gs) <- markernames
    pData(parameters(gh_pop_get_data(gs[[1]])))[,c(1,2)]

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

    ## done!

Next, let’s get rid of the “NOT4+” node and all of its descendants.

    gs_pop_remove(gs, "/Time/LD-3+/1419-3+/S/Lymph/NOT4+")

Plot the gating tree.

    plot(gs, fontsize=15, bool=T)

![](/home/emmabishop/workspace/tutorials/flow_analysis_COMPASS/README_files/figure-markdown_strict/unnamed-chunk-11-1.png)

## Save GatingSet

**Note:** this can take a while.

    save_gs(gs, here::here("out/GatingSet"))

    ## Done
    ## To reload it, use 'load_gs' function

# Get tabular data from the GatingSet

The following steps were inspired from advice in the COMPASS [issue
\#81](https://github.com/RGLab/COMPASS/issues/81#issuecomment-2484533481)
thread.

    # Load GatingSet if needed
    gs <- load_gs(here::here("out/GatingSet"))

## Get matrices of expression data

Extract metadata and take a look

    metadata <- pData(gs) %>%
      mutate(sampleName = rownames(.))

    head(metadata)

    ##                         name       EXPERIMENT NAME       $DATE SAMPLE ID PATIENT ID    STIM WELL ID PLATE NAME        sampleName
    ## 114716.fcs_445737 114716.fcs 20200605_COVID_ICS-B3 05-JUN-2020        23         23    DMSO     G02         P3 114716.fcs_445737
    ## 114718.fcs_426995 114718.fcs 20200605_COVID_ICS-B3 05-JUN-2020        23         23     SEB     G04         P3 114718.fcs_426995
    ## 114720.fcs_339914 114720.fcs 20200605_COVID_ICS-B3 05-JUN-2020        23         23    VEMP     G06         P3 114720.fcs_339914
    ## 114722.fcs_383957 114722.fcs 20200605_COVID_ICS-B3 05-JUN-2020        23         23 Spike 1     G08         P3 114722.fcs_383957
    ## 114724.fcs_273422 114724.fcs 20200605_COVID_ICS-B3 05-JUN-2020        23         23 Spike 2     G10         P3 114724.fcs_273422
    ## 114726.fcs_373199 114726.fcs 20200605_COVID_ICS-B3 05-JUN-2020        23         23    NCAP     G12         P3 114726.fcs_373199

Separate the Spike 1 stimmed and DMSO samples

    sp1_samples <- metadata %>%
      filter(STIM == "Spike 1") %>%
      pull(sampleName)
    spike1 <- gs[sp1_samples]

    dmso_samples <- metadata %>%
      filter(STIM == "DMSO") %>%
      pull(sampleName)
    dmso <- gs[dmso_samples]

Extract data for the “CD4+” population only. These are ‘cytoset’
objects, which are a collection of ’cytoframe’s.

    spike1_cd4 <- gs_pop_get_data(spike1, "4+")  # Using "+4" node name
    dmso_cd4 <- gs_pop_get_data(dmso, "4+")

We can see the sample names in the object like so:

    sampleNames(dmso_cd4)

    ## [1] "114716.fcs_445737" "114776.fcs_846625" "114752.fcs_489936" "114788.fcs_371622" "114764.fcs_435350"

Let’s take a look at one of these ‘cytoframes’, which are similar to a
matrix of expression values.

    samplename <- "114716.fcs_445737"
    cytoframe1 <- dmso_cd4[[samplename]]
    head(cytoframe1)

    ##       Time     FSC-A  FSC-H    SSC-A <B710-A> <B515-A>  <G780-A> <G660-A> <G610-A>  <G575-A> <R780-A>  <R710-A> <R660-A>
    ## [1,] 2.441 122933.04 111005 22292.88 912.8039 605.4547 220.27991 560.6566 3060.245  395.8604 1715.682  577.2411 731.5930
    ## [2,] 2.441 126327.96 110188 34016.13 597.5212 732.2466 190.20494 873.4103 2982.977  410.8742 1868.221 1142.5469 436.4252
    ## [3,] 2.442  93704.20  80479 30572.67 897.9616 617.5473  78.86456 440.2558 3033.829  995.1781 2014.507  783.3704 659.7129
    ## [4,] 2.447 109862.56  95799 24374.79 316.1867 428.8226 893.86353 997.1701 3053.399  976.4500 1855.618  381.2457 654.4659
    ## [5,] 2.449  82438.72  70628 40428.90 766.8693 484.2879 911.77679 915.0909 3070.636  630.3830 2046.922 1543.0740 401.9236
    ## [6,] 2.453  80187.60  68863 36685.29 469.5395 537.9316 433.60690 789.3820 3119.851 1252.5963 2124.699 1185.0879 831.7593
    ##       <V780-A> <V710-A> <V610-A> <V510-A> <V450-A>  <U730-A>  <U395-A>
    ## [1,]  654.3718 2341.164 1162.096 1066.730 757.4949 1392.2936  691.1917
    ## [2,] 1132.8197 1025.377 1015.737 1221.067 837.0206  779.2865  578.2829
    ## [3,] 1403.7660 2366.746 1613.703 1262.292 405.5316 2563.0171  568.5731
    ## [4,]  889.6058 2180.418 2253.460 1215.396 645.0524 2731.7053 1013.0549
    ## [5,]  890.8912 2353.598 1283.655 1216.928 919.3452 2743.6641  768.1965
    ## [6,]  580.8347 1753.941 1515.809 1387.849 771.2418  884.6489 1112.8729

Note that we have gates we don’t need for our boolean combinations, like
time and forward/side scatter (Time, FSC-A, FSC-H, SSC-A). We will get
rid of these below.

Get a list of expression matrices for stim and unstim from the
cytoset/cytoframe objects.

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

Now we have named lists of matrices of expression data, using sample
names as indices

    names(dmsolist)

    ## [1] "114716.fcs_445737" "114776.fcs_846625" "114752.fcs_489936" "114788.fcs_371622" "114764.fcs_435350"

    head(dmsolist[[1]])

    ##          TNFa    CD107a    CD154       IL2     IL17a IL4_5_13     IFNg
    ## [1,] 605.4547 220.27991 560.6566  395.8604  577.2411 731.5930 757.4949
    ## [2,] 732.2466 190.20494 873.4103  410.8742 1142.5469 436.4252 837.0206
    ## [3,] 617.5473  78.86456 440.2558  995.1781  783.3704 659.7129 405.5316
    ## [4,] 428.8226 893.86353 997.1701  976.4500  381.2457 654.4659 645.0524
    ## [5,] 484.2879 911.77679 915.0909  630.3830 1543.0740 401.9236 919.3452
    ## [6,] 537.9316 433.60690 789.3820 1252.5963 1185.0879 831.7593 771.2418

Save

    saveRDS(spikelist, here::here("out/spikelist.rds"))
    saveRDS(dmsolist, here::here("out/dmsolist.rds"))

## Create boolean combinations

    spikelist <- readRDS(here::here("out/spikelist.rds"))
    dmsolist <- readRDS(here::here("out/dmsolist.rds"))

Can take a little while

    k = 7  # Seven markers

    n_s <- CellCounts(spikelist, Combinations(k)) # Stimmed data

    n_u <- CellCounts(dmsolist, Combinations(k)) # Unstimmed data

Save count matrices and metadata

    saveRDS(n_s, here::here("out/n_s.rds"))
    saveRDS(n_u, here::here("out/n_u.rds"))
    saveRDS(metadata, here::here("out/metadata.rds"))

# Run SimpleCOMPASS

Fits a polyfunctionality response model to our data.

    fit = COMPASS::SimpleCOMPASS(n_s = n_s,
                                 n_u = n_u,
                                 meta = metadata,
                                 individual_id = "sampleName")

    ## Ordering meta, n_s and n_u by individual_id since this wasn't done.
    ## If you think this is an error, check your data and rerun the code.

    ## Setting the seed to  100

    ## Initializing parameters...

    ## Computing initial parameter estimates...

    ## Iteration 1000 of 10000.

    ## Iteration 2000 of 10000.

    ## Iteration 3000 of 10000.

    ## Iteration 4000 of 10000.

    ## Iteration 5000 of 10000.

    ## Iteration 6000 of 10000.

    ## Iteration 7000 of 10000.

    ## Iteration 8000 of 10000.

    ## Iteration 9000 of 10000.

    ## Iteration 10000 of 10000.

    ## Keeping 10000 iterations. We'll thin every 8 iterations.

    ## Burnin for 10000 iterations...

    ## Iteration 1000 of 90000.

    ## Iteration 2000 of 90000.

    ## Iteration 3000 of 90000.

    ## Iteration 4000 of 90000.

    ## Iteration 5000 of 90000.

    ## Iteration 6000 of 90000.

    ## Iteration 7000 of 90000.

    ## Iteration 8000 of 90000.

    ## Iteration 9000 of 90000.

    ## Iteration 10000 of 90000.

    ## Sampling 80000 iterations...

    ## Iteration 11000 of 90000.

    ## Iteration 12000 of 90000.

    ## Iteration 13000 of 90000.

    ## Iteration 14000 of 90000.

    ## Iteration 15000 of 90000.

    ## Iteration 16000 of 90000.

    ## Iteration 17000 of 90000.

    ## Iteration 18000 of 90000.

    ## Iteration 19000 of 90000.

    ## Iteration 20000 of 90000.

    ## Iteration 21000 of 90000.

    ## Iteration 22000 of 90000.

    ## Iteration 23000 of 90000.

    ## Iteration 24000 of 90000.

    ## Iteration 25000 of 90000.

    ## Iteration 26000 of 90000.

    ## Iteration 27000 of 90000.

    ## Iteration 28000 of 90000.

    ## Iteration 29000 of 90000.

    ## Iteration 30000 of 90000.

    ## Iteration 31000 of 90000.

    ## Iteration 32000 of 90000.

    ## Iteration 33000 of 90000.

    ## Iteration 34000 of 90000.

    ## Iteration 35000 of 90000.

    ## Iteration 36000 of 90000.

    ## Iteration 37000 of 90000.

    ## Iteration 38000 of 90000.

    ## Iteration 39000 of 90000.

    ## Iteration 40000 of 90000.

    ## Iteration 41000 of 90000.

    ## Iteration 42000 of 90000.

    ## Iteration 43000 of 90000.

    ## Iteration 44000 of 90000.

    ## Iteration 45000 of 90000.

    ## Iteration 46000 of 90000.

    ## Iteration 47000 of 90000.

    ## Iteration 48000 of 90000.

    ## Iteration 49000 of 90000.

    ## Iteration 50000 of 90000.

    ## Iteration 51000 of 90000.

    ## Iteration 52000 of 90000.

    ## Iteration 53000 of 90000.

    ## Iteration 54000 of 90000.

    ## Iteration 55000 of 90000.

    ## Iteration 56000 of 90000.

    ## Iteration 57000 of 90000.

    ## Iteration 58000 of 90000.

    ## Iteration 59000 of 90000.

    ## Iteration 60000 of 90000.

    ## Iteration 61000 of 90000.

    ## Iteration 62000 of 90000.

    ## Iteration 63000 of 90000.

    ## Iteration 64000 of 90000.

    ## Iteration 65000 of 90000.

    ## Iteration 66000 of 90000.

    ## Iteration 67000 of 90000.

    ## Iteration 68000 of 90000.

    ## Iteration 69000 of 90000.

    ## Iteration 70000 of 90000.

    ## Iteration 71000 of 90000.

    ## Iteration 72000 of 90000.

    ## Iteration 73000 of 90000.

    ## Iteration 74000 of 90000.

    ## Iteration 75000 of 90000.

    ## Iteration 76000 of 90000.

    ## Iteration 77000 of 90000.

    ## Iteration 78000 of 90000.

    ## Iteration 79000 of 90000.

    ## Iteration 80000 of 90000.

    ## Iteration 81000 of 90000.

    ## Iteration 82000 of 90000.

    ## Iteration 83000 of 90000.

    ## Iteration 84000 of 90000.

    ## Iteration 85000 of 90000.

    ## Iteration 86000 of 90000.

    ## Iteration 87000 of 90000.

    ## Iteration 88000 of 90000.

    ## Iteration 89000 of 90000.

    ## Iteration 90000 of 90000.

    ## Done!

Save the COMPASS run output.

    saveRDS(fit, file.path(here::here("out/COMPASSResult"), "COMPASSResult.rds"))

Save the Functionality and Polyfunctionality Scores.

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

    ##          sampleName          FS        PFS
    ## 1 114722.fcs_383957 0.007885827 0.03062500
    ## 2 114758.fcs_429530 0.011088189 0.04310816
    ## 3 114770.fcs_392388 0.013247244 0.05133401
    ## 4 114782.fcs_941925 0.031700000 0.12266922
    ## 5 114794.fcs_318768 0.020082677 0.07757687

# Visualize

Plot a heatmap of the mean probability of response.

    plot(fit, show_rownames = TRUE)

    ## The 'threshold' filter has removed 122 categories:
    ## TNFa&CD107a&CD154&IL2&IL17a&IL4_5_13&IFNg, !TNFa&!CD107a&CD154&IL2&IL17a&IL4_5_13&IFNg, !TNFa&CD107a&!CD154&IL2&IL17a&IL4_5_13&IFNg, TNFa&!CD107a&!CD154&IL2&IL17a&IL4_5_13&IFNg, !TNFa&!CD107a&!CD154&IL2&IL17a&IL4_5_13&IFNg, !TNFa&CD107a&CD154&!IL2&IL17a&IL4_5_13&IFNg, TNFa&!CD107a&CD154&!IL2&IL17a&IL4_5_13&IFNg, !TNFa&!CD107a&CD154&!IL2&IL17a&IL4_5_13&IFNg, TNFa&CD107a&!CD154&!IL2&IL17a&IL4_5_13&IFNg, !TNFa&CD107a&!CD154&!IL2&IL17a&IL4_5_13&IFNg, TNFa&!CD107a&!CD154&!IL2&IL17a&IL4_5_13&IFNg, !TNFa&!CD107a&!CD154&!IL2&IL17a&IL4_5_13&IFNg, TNFa&CD107a&CD154&IL2&!IL17a&IL4_5_13&IFNg, !TNFa&CD107a&CD154&IL2&!IL17a&IL4_5_13&IFNg, TNFa&!CD107a&CD154&IL2&!IL17a&IL4_5_13&IFNg, !TNFa&!CD107a&CD154&IL2&!IL17a&IL4_5_13&IFNg, TNFa&CD107a&!CD154&IL2&!IL17a&IL4_5_13&IFNg, !TNFa&CD107a&!CD154&IL2&!IL17a&IL4_5_13&IFNg, TNFa&!CD107a&!CD154&IL2&!IL17a&IL4_5_13&IFNg, !TNFa&!CD107a&!CD154&IL2&!IL17a&IL4_5_13&IFNg, TNFa&CD107a&CD154&!IL2&!IL17a&IL4_5_13&IFNg, !TNFa&CD107a&CD154&!IL2&!IL17a&IL4_5_13&IFNg, TNFa&!CD107a&CD154&!IL2&!IL17a&IL4_5_13&IFNg, !TNFa&!CD107a&CD154&!IL2&!IL17a&IL4_5_13&IFNg, TNFa&CD107a&!CD154&!IL2&!IL17a&IL4_5_13&IFNg, !TNFa&CD107a&!CD154&!IL2&!IL17a&IL4_5_13&IFNg, TNFa&!CD107a&!CD154&!IL2&!IL17a&IL4_5_13&IFNg, !TNFa&!CD107a&!CD154&!IL2&!IL17a&IL4_5_13&IFNg, !TNFa&CD107a&CD154&IL2&IL17a&!IL4_5_13&IFNg, TNFa&!CD107a&CD154&IL2&IL17a&!IL4_5_13&IFNg, !TNFa&!CD107a&CD154&IL2&IL17a&!IL4_5_13&IFNg, TNFa&CD107a&!CD154&IL2&IL17a&!IL4_5_13&IFNg, !TNFa&CD107a&!CD154&IL2&IL17a&!IL4_5_13&IFNg, TNFa&!CD107a&!CD154&IL2&IL17a&!IL4_5_13&IFNg, !TNFa&!CD107a&!CD154&IL2&IL17a&!IL4_5_13&IFNg, TNFa&CD107a&CD154&!IL2&IL17a&!IL4_5_13&IFNg, !TNFa&CD107a&CD154&!IL2&IL17a&!IL4_5_13&IFNg, TNFa&!CD107a&CD154&!IL2&IL17a&!IL4_5_13&IFNg, !TNFa&!CD107a&CD154&!IL2&IL17a&!IL4_5_13&IFNg, TNFa&CD107a&!CD154&!IL2&IL17a&!IL4_5_13&IFNg, !TNFa&CD107a&!CD154&!IL2&IL17a&!IL4_5_13&IFNg, TNFa&!CD107a&!CD154&!IL2&IL17a&!IL4_5_13&IFNg, !TNFa&!CD107a&!CD154&!IL2&IL17a&!IL4_5_13&IFNg, TNFa&CD107a&CD154&IL2&!IL17a&!IL4_5_13&IFNg, !TNFa&CD107a&CD154&IL2&!IL17a&!IL4_5_13&IFNg, TNFa&!CD107a&CD154&IL2&!IL17a&!IL4_5_13&IFNg, !TNFa&!CD107a&CD154&IL2&!IL17a&!IL4_5_13&IFNg, TNFa&CD107a&!CD154&IL2&!IL17a&!IL4_5_13&IFNg, !TNFa&CD107a&!CD154&IL2&!IL17a&!IL4_5_13&IFNg, TNFa&!CD107a&!CD154&IL2&!IL17a&!IL4_5_13&IFNg, !TNFa&!CD107a&!CD154&IL2&!IL17a&!IL4_5_13&IFNg, TNFa&CD107a&CD154&!IL2&!IL17a&!IL4_5_13&IFNg, !TNFa&CD107a&CD154&!IL2&!IL17a&!IL4_5_13&IFNg, TNFa&!CD107a&CD154&!IL2&!IL17a&!IL4_5_13&IFNg, !TNFa&!CD107a&CD154&!IL2&!IL17a&!IL4_5_13&IFNg, TNFa&CD107a&!CD154&!IL2&!IL17a&!IL4_5_13&IFNg, !TNFa&CD107a&!CD154&!IL2&!IL17a&!IL4_5_13&IFNg, TNFa&!CD107a&!CD154&!IL2&!IL17a&!IL4_5_13&IFNg, !TNFa&!CD107a&!CD154&!IL2&!IL17a&!IL4_5_13&IFNg, TNFa&CD107a&CD154&IL2&IL17a&IL4_5_13&!IFNg, !TNFa&CD107a&CD154&IL2&IL17a&IL4_5_13&!IFNg, TNFa&!CD107a&CD154&IL2&IL17a&IL4_5_13&!IFNg, !TNFa&!CD107a&CD154&IL2&IL17a&IL4_5_13&!IFNg, TNFa&CD107a&!CD154&IL2&IL17a&IL4_5_13&!IFNg, !TNFa&CD107a&!CD154&IL2&IL17a&IL4_5_13&!IFNg, TNFa&!CD107a&!CD154&IL2&IL17a&IL4_5_13&!IFNg, !TNFa&!CD107a&!CD154&IL2&IL17a&IL4_5_13&!IFNg, TNFa&CD107a&CD154&!IL2&IL17a&IL4_5_13&!IFNg, !TNFa&CD107a&CD154&!IL2&IL17a&IL4_5_13&!IFNg, TNFa&!CD107a&CD154&!IL2&IL17a&IL4_5_13&!IFNg, !TNFa&!CD107a&CD154&!IL2&IL17a&IL4_5_13&!IFNg, TNFa&CD107a&!CD154&!IL2&IL17a&IL4_5_13&!IFNg, !TNFa&CD107a&!CD154&!IL2&IL17a&IL4_5_13&!IFNg, TNFa&!CD107a&!CD154&!IL2&IL17a&IL4_5_13&!IFNg, !TNFa&!CD107a&!CD154&!IL2&IL17a&IL4_5_13&!IFNg, TNFa&CD107a&CD154&IL2&!IL17a&IL4_5_13&!IFNg, !TNFa&CD107a&CD154&IL2&!IL17a&IL4_5_13&!IFNg, TNFa&!CD107a&CD154&IL2&!IL17a&IL4_5_13&!IFNg, !TNFa&!CD107a&CD154&IL2&!IL17a&IL4_5_13&!IFNg, TNFa&CD107a&!CD154&IL2&!IL17a&IL4_5_13&!IFNg, !TNFa&CD107a&!CD154&IL2&!IL17a&IL4_5_13&!IFNg, TNFa&!CD107a&!CD154&IL2&!IL17a&IL4_5_13&!IFNg, !TNFa&!CD107a&!CD154&IL2&!IL17a&IL4_5_13&!IFNg, TNFa&CD107a&CD154&!IL2&!IL17a&IL4_5_13&!IFNg, !TNFa&CD107a&CD154&!IL2&!IL17a&IL4_5_13&!IFNg, TNFa&!CD107a&CD154&!IL2&!IL17a&IL4_5_13&!IFNg, !TNFa&!CD107a&CD154&!IL2&!IL17a&IL4_5_13&!IFNg, TNFa&CD107a&!CD154&!IL2&!IL17a&IL4_5_13&!IFNg, !TNFa&CD107a&!CD154&!IL2&!IL17a&IL4_5_13&!IFNg, TNFa&!CD107a&!CD154&!IL2&!IL17a&IL4_5_13&!IFNg, !TNFa&!CD107a&!CD154&!IL2&!IL17a&IL4_5_13&!IFNg, TNFa&CD107a&CD154&IL2&IL17a&!IL4_5_13&!IFNg, !TNFa&CD107a&CD154&IL2&IL17a&!IL4_5_13&!IFNg, TNFa&!CD107a&CD154&IL2&IL17a&!IL4_5_13&!IFNg, !TNFa&!CD107a&CD154&IL2&IL17a&!IL4_5_13&!IFNg, TNFa&CD107a&!CD154&IL2&IL17a&!IL4_5_13&!IFNg, !TNFa&CD107a&!CD154&IL2&IL17a&!IL4_5_13&!IFNg, TNFa&!CD107a&!CD154&IL2&IL17a&!IL4_5_13&!IFNg, !TNFa&!CD107a&!CD154&IL2&IL17a&!IL4_5_13&!IFNg, TNFa&CD107a&CD154&!IL2&IL17a&!IL4_5_13&!IFNg, !TNFa&CD107a&CD154&!IL2&IL17a&!IL4_5_13&!IFNg, TNFa&!CD107a&CD154&!IL2&IL17a&!IL4_5_13&!IFNg, !TNFa&!CD107a&CD154&!IL2&IL17a&!IL4_5_13&!IFNg, TNFa&CD107a&!CD154&!IL2&IL17a&!IL4_5_13&!IFNg, !TNFa&CD107a&!CD154&!IL2&IL17a&!IL4_5_13&!IFNg, TNFa&!CD107a&!CD154&!IL2&IL17a&!IL4_5_13&!IFNg, !TNFa&!CD107a&!CD154&!IL2&IL17a&!IL4_5_13&!IFNg, TNFa&CD107a&CD154&IL2&!IL17a&!IL4_5_13&!IFNg, !TNFa&CD107a&CD154&IL2&!IL17a&!IL4_5_13&!IFNg, TNFa&!CD107a&CD154&IL2&!IL17a&!IL4_5_13&!IFNg, !TNFa&!CD107a&CD154&IL2&!IL17a&!IL4_5_13&!IFNg, TNFa&CD107a&!CD154&IL2&!IL17a&!IL4_5_13&!IFNg, !TNFa&CD107a&!CD154&IL2&!IL17a&!IL4_5_13&!IFNg, TNFa&!CD107a&!CD154&IL2&!IL17a&!IL4_5_13&!IFNg, !TNFa&!CD107a&!CD154&IL2&!IL17a&!IL4_5_13&!IFNg, TNFa&CD107a&CD154&!IL2&!IL17a&!IL4_5_13&!IFNg, !TNFa&CD107a&CD154&!IL2&!IL17a&!IL4_5_13&!IFNg, TNFa&!CD107a&CD154&!IL2&!IL17a&!IL4_5_13&!IFNg, !TNFa&!CD107a&CD154&!IL2&!IL17a&!IL4_5_13&!IFNg, TNFa&CD107a&!CD154&!IL2&!IL17a&!IL4_5_13&!IFNg, !TNFa&CD107a&!CD154&!IL2&!IL17a&!IL4_5_13&!IFNg, TNFa&!CD107a&!CD154&!IL2&!IL17a&!IL4_5_13&!IFNg

![](/home/emmabishop/workspace/tutorials/flow_analysis_COMPASS/README_files/figure-markdown_strict/unnamed-chunk-29-1.png)
