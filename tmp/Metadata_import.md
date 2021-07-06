---
title: "metadata_import"
author: "Sneha & FC"
date: "June 08, 2021"
output: 
  html_document: 
    toc: yes
    keep_md: yes
---





```r
library(tidyverse)
library(phyloseq)
library(readxl)
library(here)
```


First import the phyloseq object.

```r
#get the file path of 'phyloseq.RDS' 
ps <- "data/raw/metabarcoding/physeq.RDS"

#load the phyloseq object
ps %>% 
  here::here() %>% 
  readRDS() -> physeq
```


Import the excel file containing the metadata.

```r
#get file path of the metadata file (it is an excel file)
md <- "data/raw/metabarcoding/p654_Alessia_Pennacchia_sample_sheet_file_11.01.2021.xlsx"
#load the metadata file, set na = 'NA' so that it does not read the NAs in the excel sheet as characters
md %>% 
  here::here() %>% 
  read_excel(na='NA') -> metadata  
```


We need to add the new columns in `metadata` to the existing columns in `sample_data(physeq)`. 

First alphabetically arrange `metadata` by `Sample_name` to make sure we merge the dataframes correctly

```r
metadata %>% arrange(Sample_name) -> metadata

#metadata
```


Get the `sample_data(physeq)` as a dataframe for merging

```r
physeq %>%  
  sample_data() %>% 
  as.matrix() %>%  
  as.data.frame() %>% 
  rownames_to_column(var="Sample_name") %>% 
  arrange(Sample_name)-> sample.data 

#sample.data
```

Check if the sample names match before merging


```r
all.equal(metadata$Sample_name,sample.data$Sample_name)
```

```
## [1] TRUE
```



Merge the two dataframes (we add columns 2 to 25 of `sample.data` to `metadata`)

```r
metadata %>% 
  mutate(sample.data[,2:25],) -> metadata

#check dimensions
dim(metadata)
```

```
## [1] 384  51
```


`sample_data` had `Sample_name` as rownames . In order to match the format, we make sure that `metadata` also has rownames.

```r
metadata %>% 
          column_to_rownames(var="Sample_name") -> metadata
```

Overwriting  `sample_data(physeq)` of the phyloseq object with `metadata`

```r
sample_data(physeq) <- metadata
```


Add metabolite data to sample data of phyloseq object

```r
here("data/raw/hplc Fermentation (Salvato automaticamente).xlsx") %>%
      readxl::read_xlsx(sheet = "All total") -> metabolites


physeq@sam_data %>%
  data.frame() %>%
  rownames_to_column('id') %>%
  left_join(
    metabolites,
    by = c("Sample_description" = "Sample_Id")) %>%
  column_to_rownames('id') %>% 
  sample_data() -> physeq@sam_data
```


Creating factor variables out of all the metadata columns so that the levels are in the order we specify. 


```r
sample_data(physeq)$Experiment <- fct_relevel(sample_data(physeq)$Experiment, "NTC","Mock", "HV292.1", "CCUG59168", "Cecum","Continuous","Batch") 


sample_data(physeq)$Treatment <- fct_relevel(sample_data(physeq)$Treatment,"negative","positive","STRAIN", "DONOR", "UNTREATED", "CTX", "CTX+HV292.1", "HV292.1", "VAN", "VAN+CCUG59168", "CCUG59168")


sample_data(physeq)$Reactor <- fct_relevel(sample_data(physeq)$Reactor,"negative","positive","STRAIN", "DONOR", "IR1", "IR2", "CR", "TR1", "TR2", "TR3", "TR4", "TR5", "TR6","Batch3","Batch4")


sample_data(physeq)$Reactor_Treatment <- fct_relevel(sample_data(physeq)$`Reactor_Treatment`,"negative","positive","STRAIN","DONOR", "IR1_UNTREATED","IR2_UNTREATED","CR_UNTREATED","CR_CTX","CR_VAN","TR1_CTX+HV292.1", "TR2_CTX", "TR3_HV292.1", "TR4_VAN","TR5_VAN+CCUG59168","TR6_CCUG59168","Batch3_UNTREATED","Batch3_CTX+HV292.1","Batch3_CTX","Batch3_HV292.1","Batch4_UNTREATED","Batch4_VAN","Batch4_VAN+CCUG59168","Batch4_CCUG59168")

sample_data(physeq)$Enrichment <- fct_relevel(sample_data(physeq)$Enrichment, "NotEnriched","Enriched") 

#sample_data(physeq)$Phase  
sample_data(physeq)$Phase<-factor(x=sample_data(physeq)$Phase,levels=c( "STRAIN","DONOR","Stab","Treat"),exclude=NULL)


sample_data(physeq)$Treatment2<-factor(x=sample_data(physeq)$Treatment2,levels=c("negative","positive","STRAIN", "DONOR","UNTREATED","AB","AB+E. coli","E. coli","AB+E. faecium","E. faecium"),exclude=NULL)

sample_data(physeq)$Batch_hour<-factor(x=sample_data(physeq)$Batch_hour,levels=c( "STRAIN","DONOR","0","7","24"),exclude=NULL)


sample_data(physeq)$plate <- fct_relevel(sample_data(physeq)$plate, "P1","P2","P3","P4") 
```


Check if all is well 


```r
 as(sample_data(physeq),"matrix") %>%
         data.frame(check.names=FALSE) %>%
         rownames_to_column("Sample_ID") %>% dim()
```

```
## [1] 384  64
```

Saving the new updated physeq object

```r
#store the path where you want the new physeq object
out.path <- "data/processed"
physeq_update <- physeq

saveRDS(physeq_update,file = file.path(here::here(out.path),"physeq_update_16S_chicken_08-06-21.RDS"))
```



