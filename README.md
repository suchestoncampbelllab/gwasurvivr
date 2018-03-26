---
title: "gwasurvivr Vignette"
author:
- name: Abbas Rizvi
  affiliation: The Ohio State University, Columbus, OH
- name: Ezgi Karaesmen
  affiliation: The Ohio State University, Columbus, OH
- name: Martin Morgan
  affiliation: Roswell Park Comprehensive Cancer Center, Buffalo, NY    
- name: Lara Sucheston-Campbell  
  affiliation: The Ohio State University, Columbus, OH
package: gwasurvivr
date: 'March 26, 2018'
output:
        BiocStyle::html_document2
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEncoding{UTF-8}          
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---



# Introduction
`gwasurvivr` can perform survival analysis on full chromosomes on imputed genotypes from Sanger/Michigan imputation servers or chromosome chunks imputed genotypes from IMPUTE2. This vignette is a tutorial on how to perform survival analysis (Cox proportional hazard models) on imputed GWAS data from either using the imputation software IMPUTE2 or an imputation server (Sanger/Michigan imputation servers). This package can be conveniently batched on a high performing computing cluster. 

# Getting started
Install `gwasurvivr` from the [Sucheston-Campbell Lab Github repository](http://github.com/suchestoncampbelllab/gwasurvivr) using `devtools`.


```r
devtools::install_github("suchestoncampbelllab/gwasurvivr")
```

## Dependencies
**Note**:  This package depends on `GWASTools` which uses netcdf framework on linux. Therefore, for linux users, please install `libnetcdf-dev` and `netcdf-bin` before installing `gwasurvivr`. 

CRAN packages:  
1. `ncdf4`  
2. `matrixStats`  
3. `parallel`    
4. `survival` 


```r
install.packages(c("ncdf4", "matrixStats", "parallel", "survival"))
```

Bioconductor packages:  
1. `GWASTools`  
2. `VariantAnnotation` 


```r
source("https://bioconductor.org/biocLite.R")
biocLite("GWASTools")
biocLite("VariantAnnotation")
```

## Load package

```r
library(gwasurvivr)
```

# Sanger Imputation Server
[Sanger Imputation Server](https://imputation.sanger.ac.uk/) pre-phases using either SHAPEIT or EAGLE and then imputes genotypes using PBWT algorithm. The output from Sanger imputation are `.vcf.gz` files, which are needed as input for `gwasurvivr`. The output comes in full chromosomes, which can easily be run using `gwasurvivr`. `sangerCoxSurv` is the function needed to perform survival. `sangerCoxSurv` can filter info score (imputation quality metric) and on minor allele frequency. `sangerCoxSurv` grabs the reference panel allele frequency from the VCF file (`RefPanelAF`) and also computes a minor allele frequency (MAF) for the sample genotype. The `RefPanelAF` can be filtered as the input argument `maf.filter`.

Users can select which samples that they would like to include in the analysis by providing a vector of `sample.ids`. Since the output from Sanger imputation server returns the samples as `SAMP1, ..., SAMPN`, where `N` is the total number of samples, the user must create an alternate ID for the samples of interest that correspond to the order that the samples were in their input to the Sanger imputation server (this sample order can also be found in the `.fam` file from the genotyping chip). 

To run a prespecified number of cores to use in parallel while fitting the Cox model can be set:


```r
options("gwasurvivr.cores"=4)
```
* Note: setting number of cores is not required. 

## R Session Example

```r
vcf.file <- system.file(package="gwasurvivr",
                        "extdata", 
                        "sanger.pbwt_reference_impute.vcf.gz")
pheno.fl <- system.file(package="gwasurvivr",
                        "extdata", 
                        "simulated_pheno.txt")
pheno.file <- read.table(pheno.fl, sep=" ", header=TRUE, stringsAsFactors = FALSE)
head(pheno.file)
```

```
##   ID_1  ID_2 event  time   age bmiOVWT    sex        group
## 1    1 SAMP1     0 12.00 33.93       0   male      control
## 2    2 SAMP2     1  7.61 58.71       1   male experimental
## 3    3 SAMP3     0 12.00 39.38       0 female      control
## 4    4 SAMP4     0  4.30 38.85       0   male      control
## 5    5 SAMP5     0 12.00 43.58       0   male experimental
## 6    6 SAMP6     1  2.60 57.74       0   male      control
```
Categorical variables need to be converted into indicator (dummy) variables. In this example, we will select samples from the `experimental` group will test survival only on these genotypes. The first column in the `pheno.file` needs to be sample IDs that we want to match on. 

We leverage the `tidyverse` and `magrittr` packages in this example solely for data preparation/manipulation purposes.


```r
library(tidyverse)
library(magrittr)
```


```r
# recode sex column and remove first column 
pheno.file <- pheno.file %>%
        mutate(SexFemale=case_when(sex=="female"~1L,
                                   sex=="male"~0L)) %>%
        dplyr::select(-ID_1)

# select only experimental group sample.ids
sample.ids <- pheno.file %>%
        filter(group=="experimental") %$%
        ID_2 
head(sample.ids)
```

```
## [1] "SAMP2"  "SAMP5"  "SAMP7"  "SAMP9"  "SAMP11" "SAMP12"
```

Now run survival using `sangerCoxSurv`. 


```r
sangerCoxSurv(vcf.file=vcf.file,
              pheno.file=pheno.file,
              time.to.event="time",
              event="event",
              covariates=c("age", "SexFemale", "bmiOVWT"),
              sample.ids=sample.ids,
              output.name="sanger_example",
              chunk.size=10000,
              info.filter=0.7,
              maf.filter=0.01,
              verbose=TRUE)
```

```
## Analysis started on 2018-03-26 at 11:37:40
```

```
## Analysis running for 253 samples.
```

```
## Covariates included in the models are: age, bmiOVWT, SexFemale
```

```
## If your covariates of interest are not included in the model
## please stop the analysis and make sure user defined covariates
## match the column names in the pheno.file
```

```
## Analyzing chunk 0-10000
```

```
## Analysis completed on 2018-03-26 at 11:38:01
```

```
## 7254 SNPs were removed from the analysis for not meeting the threshold criteria.
```

```
## List of removed SNPs can be found in sanger_example.MAF_INFO_removed
```

## Peek at results
We can look at the results that just came out of that run by loading the data.

```r
sanger_res <- read_tsv("sanger_example.coxph")
```


```r
dim(sanger_res)
```

```
## [1] 736  18
```

```r
sanger_res %>% 
        select(RSID, TYPED, PVALUE, COEF, HR) %>%
        arrange(PVALUE) %>% 
        head() %>%
        kable()
```



|RSID       |TYPED |    PVALUE|     COEF|       HR|
|:----------|:-----|---------:|--------:|--------:|
|rs56178399 |FALSE | 0.0001774| 1.646147| 5.186956|
|rs80175559 |FALSE | 0.0002293| 1.569031| 4.801994|
|rs74833732 |FALSE | 0.0008713| 1.368206| 3.928297|
|rs12586922 |TRUE  | 0.0009214| 1.098080| 2.998402|
|rs34610474 |FALSE | 0.0009248| 1.189763| 3.286301|
|rs4643237  |FALSE | 0.0010024| 1.130329| 3.096675|

# IMPUTE2 Imputation
## File structure
IMPUTE2 outputs 6 files, but only 2 of these files required for analysis using `gwasurvivr`:

  - Genotype file (`.impute`)  
  - Sample file (`.sample`)  
  
[More information can be read about these formats](http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html)

## R Session Example
First the user needs to find correct paths of the files and load up the covariate file. The covariate file needs to be a `data.frame`

```r
impute.file <- system.file(package="gwasurvivr",
                           "extdata",
                           "impute_example.impute2")
sample.file <- system.file(package="gwasurvivr",
                           "extdata", 
                           "impute_example.impute2_sample")
chr <- 14
covariate.file <- system.file(package="gwasurvivr", 
                       "extdata",
                       "simulated_pheno.txt")
covariate.file <- read_delim(covariate.file, delim=" ")
covariate.file <- covariate.file %>% 
        mutate(SexFemale=case_when(sex=="female"~1L,
                                   sex=="male"~0L)) %>%
        dplyr::select(-ID_1)
covariate.file %>% head
```

```
## # A tibble: 6 x 8
##   ID_2  event  time   age bmiOVWT sex    group        SexFemale
##   <chr> <int> <dbl> <dbl>   <int> <chr>  <chr>            <int>
## 1 SAMP1     0 12.0   33.9       0 male   control              0
## 2 SAMP2     1  7.61  58.7       1 male   experimental         0
## 3 SAMP3     0 12.0   39.4       0 female control              1
## 4 SAMP4     0  4.30  38.8       0 male   control              0
## 5 SAMP5     0 12.0   43.6       0 male   experimental         0
## 6 SAMP6     1  2.60  57.7       0 male   control              0
```

```r
sample.ids <- covariate.file %>%
        filter(group=="experimental") %$%
        ID_2 
```

Now we can run survival

```r
impute2CoxSurv(impute.file=impute.file,
               sample.file=sample.file,
               chr=14,
               covariate.file=covariate.file,
               sample.ids=sample.ids,
               time.to.event="time",
               event="event",
               covariates=c("age", "SexFemale", "bmiOVWT"),
               out.file="impute_example",
               maf.filter=0.01,
               info.filter=0.7,
               flip.dosage=TRUE,
               verbose=TRUE)
```

```
## Analysis started on 2018-03-26 at 11:38:02
```

```
## Determining number of SNPs and samples...
```

```
## Including all SNPs.
```

```
## scan.df not given. Assigning scanIDs automatically.
```

```
## Reading sample file...
```

```
## Reading genotype file...
```

```
## Block 1 of 2
```

```
## Block 2 of 2
```

```
## Writing annotation...
```

```
## Compressing...
```

```
## Clean up the fragments of GDS file:
##     open the file '/private/var/folders/1w/bb5rrzjn4v9bvq_hlptzq4bc0000gn/T/Rtmp3WN0Ay/994178d79fa0.gds' (13.9M)
##     # of fragments: 30
##     save to '/private/var/folders/1w/bb5rrzjn4v9bvq_hlptzq4bc0000gn/T/Rtmp3WN0Ay/994178d79fa0.gds.tmp'
##     rename '/private/var/folders/1w/bb5rrzjn4v9bvq_hlptzq4bc0000gn/T/Rtmp3WN0Ay/994178d79fa0.gds.tmp' (230.5K, reduced: 13.6M)
##     # of fragments: 14
```

```
## 5888 SNPs were removed from the analysis for not meeting the given MAF < 0.01
```

```
## List of removed SNPs are saved to impute_example.snps_removed
```

```
## Covariates included in the models are: age, SexFemale, bmiOVWT
```

```
## 253 samples are included in the analysis
```

```
## Analysis completed on 2018-03-26 at 11:38:22
```

## Peek at results

```r
impute2_res <- read_tsv("impute_example.coxph")
dim(impute2_res)
```

```
## [1] 782  17
```


```r
impute2_res %>%
        select(RSID, TYPED, PVALUE, COEF, HR) %>%
        arrange(PVALUE) %>%
        head() %>%
        kable() 
```



|RSID       |TYPED      |    PVALUE|     COEF|       HR|
|:----------|:----------|---------:|--------:|--------:|
|rs56178399 |---        | 0.0001774| 1.646147| 5.186955|
|rs80175559 |---        | 0.0002293| 1.569031| 4.801994|
|rs74833732 |---        | 0.0008713| 1.368206| 3.928297|
|rs12586922 |rs12586922 | 0.0009214| 1.098080| 2.998402|
|rs34610474 |---        | 0.0009248| 1.189763| 3.286301|
|rs4643237  |---        | 0.0010024| 1.130329| 3.096675|

The results look the same.

# Michigan Imputation Server
[Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html) outputs gunzipped variant call format files (`.vcf.gz`). These files are typically stored in full chromosome chunks. Imputation is done using the Minimac3 imputation engine. Pre-phasing is done using either HAPI-UR, SHAPEIT, or EAGLE (default is EAGLE2). 

## R Session Example
The michigan example will be added soon!

# Batch Examples

## Batch Example sangerCoxSurv
Batch jobs for multiple analyses and different subsets are easy to implement using `gwasurvivr`. This is facilitated by the package `batch`, which can internalize R variables from bash. First write an R script (e.g. `mysurvivalscript.R`) to pass in bash.


```r
## mysurvivalscript.R
library(gwasurvivr)
library(batch)
parseCommandArgs(evaluate=TRUE)

options("gwasurvivr.cores"=4)
# recode sex column and remove first column
pheno.file <- pheno.file %>%
        mutate(SexFemale=case_when(sex=="female"~1L,
                                   sex=="male"~0L)) %>%
        dplyr::select(-ID_1)

# select only experimental group
sample.ids <- pheno.file %>%
        filter(group=="experimental") %$%
        ID_2 

## -- unlist the covariates 
## (refer below to the shell script as to why we are doing this)
covariates <- covariates %>%
        str_split("_") %>%
        unlist()

sangerCoxSurv(vcf.file,
              pheno.file,,
              time.to.event,
              event,
              covariates,
              sample.ids,
              chunk.size,
              info.filter,
              maf.filter,
              output.name,
              verbose=TRUE)
```

Now we can run a shell script. This can be used well with manifest files to set up multiple runs with different outcomes and different subsets of samples. The covariates are separated by an underscore (`"_"`). This is so it can be passed properly, and also why we used `str_split` to split the covariates. 


```bash
#!/bin/bash
PATH=/path/to/dir/impute_chr

R --script ${PATH}/survival/code/mysurvivalscript.R -q --args \
        vcf.file ${PATH}/chr14.vcf.gz \
        pheno.file ${PATH}/phenotype_data/pheno.txt \
        time.to.event time \
        event event \
        covariates age_SexFemale_bmiOVWT \
        sample.ids sample.ids \
        chunk.size 10000 \
        info.filter 0.7 \
        maf.filter 0.01 \
        output.name ${PATH}/survival/results/sanger_example_output
```

The file paths above are completely arbitrary and were just used as an example of how file structure may be and where desirable output would be stored.

## Batch Example impute2CoxSurv
Exactly the same as for `sangerCoxSurv` but this time with the input arguments for `impute2CoxSurv`. 

## Batch Example michiganCoxSurv

The michigan example will be added soon!

# Session info
Here is the output of `sessionInfo()` on the system that this document was compiled:

```
## R version 3.4.3 (2017-11-30)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: OS X El Capitan 10.11.6
## 
## Matrix products: default
## BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] bindrcpp_0.2        magrittr_1.5        forcats_0.2.0      
##  [4] stringr_1.2.0       dplyr_0.7.4         purrr_0.2.4        
##  [7] readr_1.1.1         tidyr_0.8.0         tibble_1.4.2       
## [10] ggplot2_2.2.1       tidyverse_1.2.1     knitr_1.19         
## [13] gwasurvivr_0.99.0   GWASTools_1.22.0    Biobase_2.36.2     
## [16] BiocGenerics_0.22.0
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_1.3-2           rsconnect_0.8.5           
##   [3] rprojroot_1.2              DNAcopy_1.50.1            
##   [5] XVector_0.16.0             GenomicRanges_1.28.4      
##   [7] rstudioapi_0.7             mice_2.30                 
##   [9] roxygen2_6.0.1             MatrixModels_0.4-1        
##  [11] bit64_0.9-7                AnnotationDbi_1.38.2      
##  [13] lubridate_1.7.1            xml2_1.1.1                
##  [15] splines_3.4.3              ncdf4_1.16                
##  [17] mnormt_1.5-5               jsonlite_1.5              
##  [19] logistf_1.22               Rsamtools_1.28.0          
##  [21] broom_0.4.3                compiler_3.4.3            
##  [23] httr_1.3.1                 backports_1.1.0           
##  [25] assertthat_0.2.0           Matrix_1.2-12             
##  [27] lazyeval_0.2.0             cli_1.0.0                 
##  [29] htmltools_0.3.6            quantreg_5.33             
##  [31] tools_3.4.3                gtable_0.2.0              
##  [33] glue_1.2.0                 GenomeInfoDbData_0.99.0   
##  [35] reshape2_1.4.3             Rcpp_0.12.15              
##  [37] cellranger_1.1.0           Biostrings_2.44.2         
##  [39] nlme_3.1-131               rtracklayer_1.36.4        
##  [41] lmtest_0.9-35              psych_1.7.8               
##  [43] rvest_0.3.2                devtools_1.13.3           
##  [45] XML_3.98-1.9               zlibbioc_1.22.0           
##  [47] MASS_7.3-47                zoo_1.8-1                 
##  [49] scales_0.5.0               BiocStyle_2.4.1           
##  [51] BSgenome_1.44.0            VariantAnnotation_1.22.3  
##  [53] hms_0.3                    gdsfmt_1.12.0             
##  [55] SummarizedExperiment_1.6.3 sandwich_2.4-0            
##  [57] SparseM_1.77               yaml_2.1.16               
##  [59] curl_2.8.1                 memoise_1.1.0             
##  [61] biomaRt_2.32.1             rpart_4.1-11              
##  [63] stringi_1.1.6              RSQLite_2.0               
##  [65] highr_0.6                  S4Vectors_0.14.3          
##  [67] desc_1.1.1                 GenomicFeatures_1.28.4    
##  [69] BiocParallel_1.10.1        GenomeInfoDb_1.12.2       
##  [71] rlang_0.1.6                pkgconfig_2.0.1           
##  [73] commonmark_1.2             matrixStats_0.53.1        
##  [75] bitops_1.0-6               evaluate_0.10.1           
##  [77] lattice_0.20-35            bindr_0.1                 
##  [79] GenomicAlignments_1.12.1   bit_1.1-12                
##  [81] plyr_1.8.4                 R6_2.2.2                  
##  [83] IRanges_2.10.2             DelayedArray_0.2.7        
##  [85] DBI_0.7                    pillar_1.1.0              
##  [87] haven_1.1.0                foreign_0.8-69            
##  [89] withr_2.0.0                mgcv_1.8-22               
##  [91] GWASExactHW_1.01           survival_2.41-3           
##  [93] RCurl_1.95-4.8             nnet_7.3-12               
##  [95] modelr_0.1.1               crayon_1.3.4              
##  [97] utf8_1.1.3                 rmarkdown_1.8             
##  [99] grid_3.4.3                 readxl_1.0.0              
## [101] quantsmooth_1.42.0         blob_1.1.0                
## [103] git2r_0.19.0               digest_0.6.13             
## [105] stats4_3.4.3               munsell_0.4.3
```

# References
1. Terry M. Therneau and Patricia M. Grambsch (2000). Modeling Survival Data: Extending the Cox Model. Springer, New York. ISBN 0-387-98784-3.  

2. Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2017). SummarizedExperiment: SummarizedExperiment container. R package version 1.6.3.  

3. Gogarten SM, Bhangale T, Conomos MP, Laurie CA, McHugh CP, Painter I, Zheng X, Crosslin DR, Levine D, Lumley T, Nelson SC, Rice K, Shen J, Swarnkar R, Weir BS and Laurie CC (2012). “GWASTools: an R/Bioconductor package for quality control and analysis of genome-wide association studies.” Bioinformatics, 28(24), pp. 3329-3331. doi: 10.1093/bioinformatics/bts610.  

4. B. N. Howie, P. Donnelly and J. Marchini (2009) A flexible and accurate genotype imputation method for the next generation of genome-wide association studies. PLoS Genetics 5(6): e1000529  

5. Das S, Forer L, Schönherr S, Sidore C, Locke AE, Kwong A, Vrieze S, Chew EY, Levy S, McGue M, Schlessinger D, Stambolian D, Loh PR, Iacono WG, Swaroop A, Scott LJ, Cucca F, Kronenberg F, Boehnke M, Abecasis GR, Fuchsberger C. **Next-generation genotype imputation service and methods.** Nature Genetics 48, 1284–1287 (2016). [27571263](https://www.ncbi.nlm.nih.gov/pubmed/27571263)  

6. Efficient haplotype matching and storage using the Positional Burrows-Wheeler Transform (PBWT)", Richard Durbin Bioinformatics 30:1266-72 (2014).   


