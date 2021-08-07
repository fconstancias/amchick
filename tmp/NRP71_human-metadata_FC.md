---
title: "NRP72 -metadata - humann - chiken "
author: "Florentin Constancias"
date: "August 04, 2021"
output: 
  html_document: 
    toc: yes
    keep_md: yes

---



#### Load required packages


```r
library(tidyverse)
library(phyloseq)
library(microbiome)
library(here)
options(getClass.msg=FALSE) # https://github.com/epurdom/clusterExperiment/issues/66
#this fixes an error message that pops up because the class 'Annotated' is defined in two different packages
```


# Import phyloseq object


```r
ps = "data/raw/metabarcoding/ps_silva_dada2_human-chicken.RDS"

ps %>% 
  here::here() %>%
  readRDS()  -> physeq

physeq$physeq -> physeq

physeq
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1155 taxa and 600 samples ]
## tax_table()   Taxonomy Table:    [ 1155 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1155 tips and 1154 internal nodes ]
## refseq()      DNAStringSet:      [ 1155 reference sequences ]
```


```r
physeq %>% 
  sample_names() %>% 
  sort() %>% 
  head()
```

```
## [1] "AB24-1-S55"  "AB24-1E-S38" "AB24-2-S48"  "AB24-2E-S23" "AB24-3-S71" 
## [6] "AB24-3E-S25"
```



```r
here::here("data/raw/metabarcoding/merged_chicken_human_04.08.2021.tsv") %>%
                        readr::read_tsv() %>% 
  pull("sample") %>% 
  sort() %>% 
  head()
```

```
## 
## ── Column specification ────────────────────────────────────────────────────────
## cols(
##   .default = col_double(),
##   sample = col_character(),
##   Sample_description = col_character(),
##   I7_Index_ID = col_character(),
##   index = col_character(),
##   I5_Index_ID = col_character(),
##   index2 = col_character(),
##   Description2 = col_character(),
##   Experiment = col_character(),
##   Reactor = col_character(),
##   Treatment = col_character(),
##   Enrichment = col_character(),
##   Phase = col_character(),
##   Treatment2 = col_character(),
##   Date = col_character(),
##   Paul = col_character(),
##   Reactor_Treatment = col_character(),
##   Model = col_character(),
##   Antibiotic = col_character()
## )
## ℹ Use `spec()` for the full column specifications.
```

```
## [1] "AB24-1-S55"  "AB24-1E-S38" "AB24-2-S48"  "AB24-2E-S23" "AB24-3-S71" 
## [6] "AB24-3E-S25"
```


```r
source("https://raw.githubusercontent.com/fconstancias/metabarcodingRpipeline/dev/scripts/functions_export_simplified.R")

# ps@sam_data = NULL

physeq %>%
  physeq_add_metadata(physeq = .,
                      metadata = here::here("data/raw/metabarcoding/merged_chicken_human_04.08.2021.tsv") %>%
                        readr::read_tsv(), sample_column = "sample") -> physeq_meta
```

```
## 
## ── Column specification ────────────────────────────────────────────────────────
## cols(
##   .default = col_double(),
##   sample = col_character(),
##   Sample_description = col_character(),
##   I7_Index_ID = col_character(),
##   index = col_character(),
##   I5_Index_ID = col_character(),
##   index2 = col_character(),
##   Description2 = col_character(),
##   Experiment = col_character(),
##   Reactor = col_character(),
##   Treatment = col_character(),
##   Enrichment = col_character(),
##   Phase = col_character(),
##   Treatment2 = col_character(),
##   Date = col_character(),
##   Paul = col_character(),
##   Reactor_Treatment = col_character(),
##   Model = col_character(),
##   Antibiotic = col_character()
## )
## ℹ Use `spec()` for the full column specifications.
```

```r
physeq;physeq_meta
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1155 taxa and 600 samples ]
## tax_table()   Taxonomy Table:    [ 1155 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1155 tips and 1154 internal nodes ]
## refseq()      DNAStringSet:      [ 1155 reference sequences ]
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1155 taxa and 600 samples ]
## sample_data() Sample Data:       [ 600 samples by 59 sample variables ]
## tax_table()   Taxonomy Table:    [ 1155 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1155 tips and 1154 internal nodes ]
## refseq()      DNAStringSet:      [ 1155 reference sequences ]
```



```r
source("https://raw.githubusercontent.com/fconstancias/metabarcodingRpipeline/dev/scripts/functions_export_simplified.R")

# ps@sam_data = NULL

physeq %>%
  physeq_add_metadata(physeq = .,
                      metadata = here::here("data/raw/metabarcoding/merged_chicken_human_04.08.2021.tsv") %>%
                        readr::read_tsv() %>% 
    mutate(Reactor_Treatment_Dose = ifelse(!is.na(`Antibiotic_mg/mL`),
       paste0(Reactor_Treatment, `Antibiotic_mg/mL`), Reactor_Treatment),
    Treatment_Dose = ifelse(!is.na(`Antibiotic_mg/mL`),
       paste0(Treatment, `Antibiotic_mg/mL`), Reactor_Treatment))
  , sample_column = "sample") -> physeq_meta
```

```
## 
## ── Column specification ────────────────────────────────────────────────────────
## cols(
##   .default = col_double(),
##   sample = col_character(),
##   Sample_description = col_character(),
##   I7_Index_ID = col_character(),
##   index = col_character(),
##   I5_Index_ID = col_character(),
##   index2 = col_character(),
##   Description2 = col_character(),
##   Experiment = col_character(),
##   Reactor = col_character(),
##   Treatment = col_character(),
##   Enrichment = col_character(),
##   Phase = col_character(),
##   Treatment2 = col_character(),
##   Date = col_character(),
##   Paul = col_character(),
##   Reactor_Treatment = col_character(),
##   Model = col_character(),
##   Antibiotic = col_character()
## )
## ℹ Use `spec()` for the full column specifications.
```

```r
# 
# meta %>% 
#     mutate(Reactor_Treatment_Dose = ifelse(!is.na(`Antibiotic_mg/mL`), 
#        paste0(Reactor_Treatment, `Antibiotic_mg/mL`), Reactor_Treatment)
#     Treatment_Dose = ifelse(!is.na(`Antibiotic_mg/mL`), 
#        paste0(Treatment, `Antibiotic_mg/mL`), Reactor_Treatment))

# ifelse(!is.na(sample_data(physeq_meta)$Antibiotic_mg.mL), 
#        sample_data(physeq_meta)$Antibiotic_mg.mL, "") %>% 
# mutate(gradebook, Pass.Fail = ifelse(grade > 60, "Pass", "Fail"))



physeq;physeq_meta
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1155 taxa and 600 samples ]
## tax_table()   Taxonomy Table:    [ 1155 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1155 tips and 1154 internal nodes ]
## refseq()      DNAStringSet:      [ 1155 reference sequences ]
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1155 taxa and 600 samples ]
## sample_data() Sample Data:       [ 600 samples by 61 sample variables ]
## tax_table()   Taxonomy Table:    [ 1155 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1155 tips and 1154 internal nodes ]
## refseq()      DNAStringSet:      [ 1155 reference sequences ]
```

```r
physeq_meta %>% 
  saveRDS(here::here("data/raw/metabarcoding/ps_silva_dada2_human_chicken_meta.RDS"))
```

600 samples initially and 471 with metadata ????


```r
intersect(
  sample_names(physeq_meta),
  sample_names(physeq)) %>% 
  length()
```

```
## [1] 600
```


```r
difference <- function(x, y) {
c(setdiff(x, y), setdiff(y, x))
}

difference(
  sample_names(physeq),
  sample_names(physeq_meta))
```

```
## character(0)
```

Above missing from the metadata but present in the initiall phyloseq object.




```r
# physeq_meta %>% 
#   add_phylogeny_to_phyloseq(export = FALSE) %>% 
#   saveRDS(here::here("data/raw/metabarcoding/ps_silva_dada2_humanonly_phylo_meta.RDS"))
```


```r
sessionInfo()
```

```
## R version 4.0.2 (2020-06-22)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS High Sierra 10.13.6
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] here_1.0.1           microbiome_1.10.0    phyloseq_1.34.0     
##  [4] forcats_0.5.0        stringr_1.4.0        dplyr_1.0.4         
##  [7] purrr_0.3.4          readr_1.4.0          tidyr_1.1.2         
## [10] tibble_3.0.6         ggplot2_3.3.3        tidyverse_1.3.0.9000
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-149        fs_1.5.0            lubridate_1.7.9    
##  [4] progress_1.2.2      httr_1.4.2          rprojroot_2.0.2    
##  [7] tools_4.0.2         backports_1.2.1     R6_2.5.0           
## [10] vegan_2.5-7         DBI_1.1.1           BiocGenerics_0.34.0
## [13] mgcv_1.8-32         colorspace_2.0-0    permute_0.9-5      
## [16] ade4_1.7-16         withr_2.4.1         tidyselect_1.1.0   
## [19] prettyunits_1.1.1   compiler_4.0.2      cli_2.3.0          
## [22] rvest_0.3.6         Biobase_2.50.0      xml2_1.3.2         
## [25] scales_1.1.1        digest_0.6.27       rmarkdown_2.4      
## [28] XVector_0.28.0      pkgconfig_2.0.3     htmltools_0.5.1.1  
## [31] dbplyr_1.4.4        rlang_0.4.10        readxl_1.3.1       
## [34] rstudioapi_0.13     generics_0.1.0      jsonlite_1.7.2     
## [37] magrittr_2.0.1      biomformat_1.7.0    Matrix_1.2-18      
## [40] Rcpp_1.0.6          munsell_0.5.0       S4Vectors_0.26.1   
## [43] Rhdf5lib_1.10.1     ape_5.4-1           lifecycle_1.0.0    
## [46] stringi_1.5.3       yaml_2.2.1          MASS_7.3-52        
## [49] zlibbioc_1.34.0     Rtsne_0.15          rhdf5_2.32.4       
## [52] plyr_1.8.6          grid_4.0.2          blob_1.2.1         
## [55] parallel_4.0.2      crayon_1.4.1        lattice_0.20-41    
## [58] Biostrings_2.56.0   haven_2.3.1         splines_4.0.2      
## [61] multtest_2.44.0     hms_1.0.0           knitr_1.31         
## [64] pillar_1.4.7        igraph_1.2.6        reshape2_1.4.4     
## [67] codetools_0.2-16    stats4_4.0.2        reprex_0.3.0       
## [70] glue_1.4.2          evaluate_0.14       data.table_1.13.6  
## [73] modelr_0.1.8        vctrs_0.3.6         foreach_1.5.1      
## [76] cellranger_1.1.0    gtable_0.3.0        assertthat_0.2.1   
## [79] xfun_0.21           broom_0.7.2         survival_3.2-3     
## [82] iterators_1.0.13    IRanges_2.22.2      cluster_2.1.0      
## [85] ellipsis_0.3.1
```

