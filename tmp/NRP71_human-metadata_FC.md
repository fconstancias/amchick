---
title: "NRP72 -metadata - humann - chiken "
author: "Florentin Constancias"
date: "December 01, 2021"
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
# here::here("data/raw/metabarcoding/merged_chicken_human_04.08.2021.tsv") %>%
here::here("data/raw/metabarcoding/25.11.2021_metadata_updated.tsv") %>%
                        readr::read_tsv() %>% 
  pull("sample") %>% 
  sort() %>% 
  head()
```

```
## Rows: 604 Columns: 61
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (19): sample, Sample_description, I7_Index_ID, index, I5_Index_ID, index...
## dbl (42): input, filtered, denoisedF, denoisedR, merged, tabled, filtered_pc...
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
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
                      metadata = here::here("data/raw/metabarcoding/25.11.2021_metadata_updated.tsv") %>%
                        readr::read_tsv(), sample_column = "sample") -> physeq_meta
```

```
## Rows: 604 Columns: 61
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (19): sample, Sample_description, I7_Index_ID, index, I5_Index_ID, index...
## dbl (42): input, filtered, denoisedF, denoisedR, merged, tabled, filtered_pc...
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
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
## sample_data() Sample Data:       [ 600 samples by 61 sample variables ]
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
## Rows: 600 Columns: 61
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (19): sample, Sample_description, I7_Index_ID, index, I5_Index_ID, index...
## dbl (42): input, filtered, denoisedF, denoisedR, merged, tabled, filtered_pc...
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
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
## sample_data() Sample Data:       [ 600 samples by 63 sample variables ]
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
## R version 4.1.2 (2021-11-01)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Mojave 10.14.6
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] here_1.0.1           microbiome_1.14.0    phyloseq_1.36.0     
##  [4] forcats_0.5.1        stringr_1.4.0        dplyr_1.0.7         
##  [7] purrr_0.3.4          readr_2.1.0          tidyr_1.1.4         
## [10] tibble_3.1.6         ggplot2_3.3.5        tidyverse_1.3.1.9000
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-153           bitops_1.0-7           fs_1.5.0              
##  [4] bit64_4.0.5            lubridate_1.8.0        httr_1.4.2            
##  [7] rprojroot_2.0.2        GenomeInfoDb_1.28.4    tools_4.1.2           
## [10] backports_1.4.0        bslib_0.3.1            vegan_2.5-7           
## [13] utf8_1.2.2             R6_2.5.1               mgcv_1.8-38           
## [16] DBI_1.1.1              BiocGenerics_0.38.0    colorspace_2.0-2      
## [19] permute_0.9-5          rhdf5filters_1.4.0     ade4_1.7-18           
## [22] withr_2.4.2            tidyselect_1.1.1       bit_4.0.4             
## [25] compiler_4.1.2         cli_3.1.0              rvest_1.0.2           
## [28] Biobase_2.52.0         xml2_1.3.2             sass_0.4.0            
## [31] scales_1.1.1           digest_0.6.28          rmarkdown_2.11        
## [34] XVector_0.32.0         pkgconfig_2.0.3        htmltools_0.5.2       
## [37] dbplyr_2.1.1           fastmap_1.1.0          rlang_0.4.12          
## [40] readxl_1.3.1           rstudioapi_0.13        jquerylib_0.1.4       
## [43] generics_0.1.1         jsonlite_1.7.2         vroom_1.5.6           
## [46] RCurl_1.98-1.5         magrittr_2.0.1         GenomeInfoDbData_1.2.6
## [49] biomformat_1.20.0      Matrix_1.3-4           Rcpp_1.0.7            
## [52] munsell_0.5.0          S4Vectors_0.30.2       Rhdf5lib_1.14.2       
## [55] fansi_0.5.0            ape_5.5                lifecycle_1.0.1       
## [58] stringi_1.7.5          yaml_2.2.1             MASS_7.3-54           
## [61] zlibbioc_1.38.0        Rtsne_0.15             rhdf5_2.36.0          
## [64] plyr_1.8.6             grid_4.1.2             parallel_4.1.2        
## [67] crayon_1.4.2           lattice_0.20-45        splines_4.1.2         
## [70] Biostrings_2.60.2      haven_2.4.3            multtest_2.48.0       
## [73] hms_1.1.1              knitr_1.36             pillar_1.6.4          
## [76] igraph_1.2.6           reshape2_1.4.4         codetools_0.2-18      
## [79] stats4_4.1.2           reprex_2.0.1           glue_1.5.0            
## [82] evaluate_0.14          data.table_1.14.2      modelr_0.1.8          
## [85] vctrs_0.3.8            tzdb_0.2.0             foreach_1.5.1         
## [88] cellranger_1.1.0       gtable_0.3.0           assertthat_0.2.1      
## [91] xfun_0.28              broom_0.7.10           survival_3.2-13       
## [94] iterators_1.0.13       IRanges_2.26.0         cluster_2.1.2         
## [97] ellipsis_0.3.2
```

