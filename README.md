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
date: 'May 02, 2018'
output: 
  BiocStyle::html_document:
    toc_float: true
vignette: >
  %\VignetteIndexEntry{gwasurvivr Vignette}
  %\VignetteEngine{rmarkdown::render}
  %\VignetteEncoding{UTF-8}
---

# Introdcution
`gwasurvivr` can be used to perform survival analyses of imputed genotypes from Sanger and Michigan imputation servers and IMPUTE2 software. This vignette is a tutorial on how to perform these analyses. This package can be run locally on a Linux, Mac OS X, Windows or conveniently batched on a high performing computing cluster. `gwasurvivr` iteratively processes the data in chunks and therefore intense memory requirements are not necessary.      
`gwasurvivr` package comes with three main functions to perform survival analyses using Cox proportional hazard (Cox PH) models depending on the imputation method used to generate the genotype data:    

1. `michiganCoxSurv`: Performs survival analysis on imputed genetic data stored in compressed VCF files generated via Michigan imputation server.    
2. `sangerCoxSurv`:  Performs survival analysis on imputed genetic data stored in compressed VCF files generated via Sanger imputation server.    
3. `impute2CoxSurv`: Performs survival analysis on imputed genetic data from IMPUTE2 output.        
4. `gdsCoxSurv`:  For files that are already in GDS format (originally in IMPUTE2 format), users can provide a path to their GDS file and perform survival analysis and avoid having to recompress their files each run.
5. `plinkCoxSurv`: For directly typed data (or imputed data that is thresholded in plink) that are plink format (.bed, .bim, .fam files), users can can perform survival analysis.

All functions fit a Cox PH model to each SNP including other user defined covariates and will save the results as a text file directly to disk that contains survival analysis results. `gwasurvivr` functions can also test for interaction of SNPs with a given covariate. See examples for further details. 

<<<<<<< HEAD
# Installation
This package is currently available on [Bioconductor devel branch](https://bioconductor.org/packages/devel/bioc/html/gwasurvivr.html) or by using `devtools` library for R >= 3.4 and going to the Sucheston Campbell Lab GitHub repository (this page).  If using R 3.5, use `BiocManager` to install the package, if using R >= 3.4, `BiocInstaller` or `biocLite` can be used.
=======
## Main input arguments

All three functions require the following main arguments:

* `vcf.file`: A character string giving the path to genotype data file (`impute.file` for IMPUTE2)
* `covariate.file`: A data frame comprising sample IDs (that match to the genotype data), phenotype (time, event) and additional covariate data
* `id.column`: A character string providing exact match to sample ID column in covariate.file
* `time.to.event`: A character string that matches time column name in covariate.file
* `event`: Character string that matches event column name in covariate.file
* `covariates`: Character vector with matching column names in covariate.file of covariates of interest
* `out.file`: A character string giving output name

Further arguments can be passed depending on the user preference. For example, user can define minor allele frequency (MAF) or info score threshold to filter out SNPs that have low MAF or info score to avoid false-positive signals and to reduce computing time. User can also define a subset of samples to be analyzed by providing the set of sample IDs. Users can also control how chunk size -- the number of rows (SNPs) to include for each iteration. 

**IMPORTANT: In the `covariate.file`, categorical variables need to be converted to indicator (dummy) variables and be of class `numeric`. Ordinal variables represented as characters, ie "low", "medium" and "high" should be converted to the appropriate numeric values as well.**

## Main output format
The output for the 3 main functions in `gwasurvivr` are very similar but with subtle differences. In general the output includes the following main fields: RSID, TYPED, CHR, POS, REF, ALT, Allele Frequencies\*, INFO\*, PVALUES, HRs, HR confidence intervals, coefficient estimates, standard errors, Z-statistic, N, and NEVENT. Allele frequencies and INFO differ by the input imputation software. 

**Note: Invoking the `inter.term` argument for any of the functions will make the PVALUE and HRs and HR confidence intervals represent the INTERACTION term and not the SNP alone.**

The non-software specific fields are summarized below:


|Column        |Description                                           |
|:-------------|:-----------------------------------------------------|
|RSID          |SNP ID                                                |
|CHR           |Chromosome number                                     |
|POS           |Genomic Position (BP)                                 |
|REF           |Reference Allele                                      |
|ALT           |Alternate Allele                                      |
|SAMP_FREQ_ALT |Alternate Allele frequency in sample being tested     |
|SAMP_MAF      |Minor allele frequency in sample being tested         |
|PVALUE        |P-value of single SNP or interaction term             |
|HR            |Hazard Ratio (HR)                                     |
|HR_lowerCI    |Lower bound 95% CI of HR                              |
|HR_upperCI    |Upper bound 95% CI of HR                              |
|COEF          |Estimated coefficient of SNP                          |
|SE.COEF       |Standard error of coefficient estimate                |
|Z             |Z-statistic                                           |
|N             |Number of individuals in sample being tested          |
|NEVENT        |Number of events that occurred in sample being tested |

The software specific fields are summarized below:  
1. `michiganCoxSurv` unique output columns are AF, MAF, INFO, ER2. They are summarized below.  


|Column |Description                                                   |
|:------|:-------------------------------------------------------------|
|TYPED  |Imputation status: TRUE (SNP IS TYPED)/FALSE (SNP IS IMPUTED) |
|AF     |Minimac3 output Alternate Allele Frequency                    |
|MAF    |Minimac3 output of Minor Allele Frequency                     |
|INFO   |Imputation INFO score (minimac3 $R^2$)                        |
|ER2    |Minimac3 ouput empirical $R^2$                                |

Please see [Minimac3 Info File](https://genome.sph.umich.edu/wiki/Minimac3_Info_File) for details on output

2. `sangerCoxSurv`  


|Column     |Description                                                   |
|:----------|:-------------------------------------------------------------|
|TYPED      |Imputation status: TRUE (SNP IS TYPED)/FALSE (SNP IS IMPUTED) |
|RefPanelAF |HRC Reference Panel Allele Frequency                          |
|INFO       |Imputation INFO score from PBWT                               |


3. `impute2CoxSurv`  


|Column      |Description                                                                                      |
|:-----------|:------------------------------------------------------------------------------------------------|
|TYPED       |`---` is imputed, repeated RSID is typed                                                         |
|A0          |Allele coded 0 in IMPUTE2                                                                        |
|A1          |Allele coded 1 in IMPUTE2                                                                        |
|exp_freq_A1 |Expected allele frequency of alelle code A1                                                      |
|INFO        |Imputation INFO score computed based on ratio of empircal and expected variance in dosage format |


More statistics can be printed out by invoking the `print.covs` argument and setting it to `print.covs=all` (single SNP/SNP\*covariate interaction) or `print.covs=some` (SNP\*covariate ineraction). These options are available mainly for modeling purposes (covariate selection) and aren't suggested for very large analyses as it will greatly increase the number of columns in the output, depending on how many covariates users are adjusting for. 

# Getting started
Install `gwasurvivr` from the [Sucheston-Campbell Lab Github repository](http://github.com/suchestoncampbelllab/gwasurvivr) using `devtools`.


```r
devtools::install_github("suchestoncampbelllab/gwasurvivr")
```

## Dependencies
**Note**:  This package depends on `GWASTools` which uses netcdf framework on linux. Therefore, for linux users, please install `libnetcdf-dev` and `netcdf-bin` before installing `gwasurvivr`. These linux libraries may already installed on an academic computing cluster. 

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
3. `SummarizedExperiment`  


```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("GWASTools")
BiocManager::install("VariantAnnotation")
BiocManager::install("SummarizedExperiment")
```

Load `gwasurvivr`.


```r
library(gwasurvivr)
```

## User settings: parallelization setup

`gwasurvivr` uses `parallel` package for its internal parallelization to fit the Cox PH models. Users are not required to define a parallelization setup, by default `gwasurvivr` functions will detect the user's operating system and set the cluster object to `FORK` if the platform is Linux/OS X and to `SOCK` if Windows. However, parallelization settings can be modified by the user if needed. Users are given two ways to define their cluster settings:  

**1. Setting the number of cores to be used:**

Linux/OS X users can run analyses on a prespecified number of cores by setting the option in the R session as shown below. This option should be defined in the R session before running any of the `gwasurvivr` functions. This option is not available to Windows users.   


```r
options("gwasurvivr.cores"=4)
```

**2. Providing a user defined cluster object**

To modify more settings, users can also provide a "cluster object" to any of the `gwasurvivr` functions. The cluster object can be generated via `makeCluster`, `makePSOCKcluster`, `makeForkCluster` functions from `parallel` package or similar cluster object functions from `snow` or `snowfall` packages. This method can be applied on any operating system. User can create a cluster object before running any of the functions and pass the cluster object to the `clusterObj` argument as shown below. For further details see `??parallel::parallel`.


```r
library(parallel)
cl <- makeCluster(detectCores())

impute2CoxSurv(..., clusterObj=cl)
sangerCoxSurv(..., clusterObj=cl)
michiganCoxSurv(..., clusterObj=cl)
```

# R Session Examples
While we use `tidyverse` and `magrittr` packages in these example for data preparation/manipulation purposes, any of the data pre-processing steps can be done in base R.


```r
library(tidyverse)
library(magrittr)
```

##  Michigan Imputation Server
[Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html) pre-phases typed genotypes using HAPI-UR, SHAPEIT, or EAGLE (default is EAGLE2), imputes using Minimac3 imputation engine and outputs Blocked GNU Zip Format VCF files (`.vcf.gz`). Just as with Sanger these `.vcf.gz` files are used as input for `gwasurvivr`. The function, `michiganCoxSurv`  uses a modification of Cox proportional hazard regression from the R library `survival`.  Minimac uses slightly different metrics to assess imputation quality ($R^2$ versus info score) and complete details as to minimac output are available on the [Minimac3 Wikipage](https://genome.sph.umich.edu/wiki/Minimac3_Imputation_Cookbook).

The function, `michiganCoxSurv`  uses a modification of cox proportional hazard regression from the R library `survival`. Built specifically for genetic data, `michiganCoxSurv` allows the user to filter on info score (imputation quality metric) and minor allele frequency from the reference panel used for imputation using `RefPanelAF` as the input arguement for `maf.filter`. Users are also provided with the sample minor allele frequency in the `sangerCoxSurv` output. 

Samples can be selected for analyses by providing a vector of `sample.ids`. The output from Sanger imputation server returns the samples as `SAMP1, ..., SAMPN`, where `N` is the total number of samples. The sample order corresponds to the sample order in the vcf file you used for imputation. Note, sample order can also be found in the `.fam` file if genotyping data were initially in `.bed`, `.bim` and `.fam` (PLINK) format prior to conversion to VCF. If no sample list is specified all samples are included in the analyses.


```r
vcf.file <- system.file(package="gwasurvivr",
                        "extdata", 
                        "michigan.chr14.dose.vcf.gz")
pheno.fl <- system.file(package="gwasurvivr",
                        "extdata", 
                        "simulated_pheno.txt")
pheno.file <- read_delim(pheno.fl, delim=" ", col_names=TRUE)
pheno.file %>%
    head()
```

```
## # A tibble: 6 x 8
##    ID_1 ID_2  event  time   age DrugTxYes sex    group       
##   <int> <chr> <int> <dbl> <dbl>     <int> <chr>  <chr>       
## 1     1 SAMP1     0 12     33.9         0 male   control     
## 2     2 SAMP2     1  7.61  58.7         1 male   experimental
## 3     3 SAMP3     0 12     39.4         0 female control     
## 4     4 SAMP4     0  4.3   38.8         0 male   control     
## 5     5 SAMP5     0 12     43.6         0 male   experimental
## 6     6 SAMP6     1  2.6   57.7         0 male   control
```

```r
# recode sex column and remove first column 
pheno.file <- pheno.file %>%
        mutate(SexFemale=if_else(sex=="female", 1L, 0L))

# select only experimental group sample.ids
sample.ids <- pheno.file %>%
        filter(group=="experimental") %$%
        ID_2 
head(sample.ids)
```

```
## [1] "SAMP2"  "SAMP5"  "SAMP7"  "SAMP9"  "SAMP11" "SAMP12"
```

In this example, we will select samples from the `experimental` group and will test survival only on these patients. The first column in the `pheno.file` are sample IDs (we will match on these). We include `age`, `DrugTxYes`, and `sex` in the survival model as covariates. 

We perform the analysis using the `experimental` group to demonstrate how one may want to prepare their data if not initially all samples are patients or cases (i.e. a case-control study and survival of cases is of interest). We also are showing how the IDs (`sample.ids`) need to be a vector of class `character`. The `chunk.size` refers to size of each data chunk read in and is defaulted to 10,000 rows, users can customize that to their needs. The larger the `chunk.size` the more memory (RAM) required to run the analysis. The recommended `chunk.size=500` and probably should not exceed `chunk.size=5000`. 

By default survival estimates and pvalues for the SNP adjusted for other covariates are outputted (`print.covs='only'`), however users can select `print.covs=all` to get the coefficient estimates for covariates included in the model. Depending on the number of covariates included this can add substantially to output file size. 
Next we run `michiganCoxSurv` with the default, `print.covs="only"`, load the results into R and provide descriptions of output by column. We will then run the analysis again using `print.covs="all"`. `verbose=TRUE` is used for these examples so the function display messages while running.

Use `?michiganCoxSurv` for argument specific documentation.


### Single SNP analysis

`print.covs="only"`


```r
michiganCoxSurv(vcf.file=vcf.file,
                covariate.file=pheno.file,
                id.column="ID_2",
                sample.ids=sample.ids,
                time.to.event="time",
                event="event",
                covariates=c("age", "SexFemale", "DrugTxYes"),
                inter.term=NULL,
                print.covs="only",
                out.file="michigan_only",
                info.filter=0.3,
                maf.filter=0.005,
                chunk.size=500,
                verbose=TRUE,
                clusterObj=NULL)
```



```
## Analysis started on 2018-05-02 at 14:00:02
```

```
## Covariates included in the models are: age, DrugTxYes, SexFemale
```

```
## 253 samples are included in the analysis
```

```
## Analyzing chunk 0-500
```

```
## Analyzing chunk 500-1000
```

```
## Analysis completed on 2018-05-02 at 14:00:07
```

```
## 240 SNPs were removed from the analysis for not meeting the threshold criteria.
```

```
## List of removed SNPs can be found in /var/folders/1w/bb5rrzjn4v9bvq_hlptzq4bc0000gn/T//Rtmp3LPLbf/michigan_only3e6b7063c12f.snps_removed
```

```
## 260 SNPs were analyzed in total
```

```
## The survival output can be found at /var/folders/1w/bb5rrzjn4v9bvq_hlptzq4bc0000gn/T//Rtmp3LPLbf/michigan_only3e6b7063c12f.coxph
```

Here we load the data and glimpse the first few values in each column outputted from the SNP*interaction survival analyses using `print.covs="only"`.


```r
michigan_only <- read_tsv("michigan_only.coxph")
```




```r
michigan_only %>% 
    head() %>%
    glimpse()
```

```
## Observations: 6
## Variables: 21
## $ RSID          <chr> "rs201487625", "rs28881575", "rs77206417", "rs13...
## $ TYPED         <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
## $ CHR           <int> 14, 14, 14, 14, 14, 14
## $ POS           <int> 19258506, 19264217, 19264589, 19264800, 19264875...
## $ REF           <chr> "A", "A", "C", "C", "G", "T"
## $ ALT           <chr> "G", "G", "T", "T", "A", "A"
## $ AF            <dbl> 0.7533720, 0.9868190, 0.0463043, 0.0105328, 0.31...
## $ MAF           <dbl> 0.2466280, 0.0131810, 0.0463043, 0.0105328, 0.31...
## $ SAMP_FREQ_ALT <dbl> 0.7861, 0.9822, 0.0363, 0.0159, 0.3101, 0.0199
## $ SAMP_MAF      <dbl> 0.2139, 0.0178, 0.0363, 0.0159, 0.3101, 0.0199
## $ INFO          <dbl> 0.925998, 0.647238, 0.487405, 0.378574, 0.583201...
## $ ER2           <dbl> NA, NA, NA, NA, NA, NA
## $ PVALUE        <dbl> 0.5221423, 0.3859810, 0.7960395, 0.5338763, 0.83...
## $ HR            <dbl> 0.8965579, 2.0522850, 1.1445003, 1.5619517, 1.03...
## $ HR_lowerCI    <dbl> 0.6417414, 0.4039371, 0.4112794, 0.3832609, 0.71...
## $ HR_upperCI    <dbl> 1.252554, 10.427054, 3.184893, 6.365620, 1.51236...
## $ COEF          <dbl> -0.10919241, 0.71895381, 0.13496810, 0.44593611,...
## $ SE.COEF       <dbl> 0.1706007, 0.8293112, 0.5221687, 0.7168242, 0.19...
## $ Z             <dbl> -0.6400466, 0.8669288, 0.2584761, 0.6220997, 0.2...
## $ N             <int> 253, 253, 253, 253, 253, 253
## $ NEVENT        <int> 100, 100, 100, 100, 100, 100
```

### SNP with covariate interaction
A SNP*covariate interaction can be implemented using the `inter.term` argument. In this example, we will use `DrugTxYes` from the covariate file as the covariate we want to test for interaction with the SNP. 

`print.covs="only"`


```r
michiganCoxSurv(vcf.file=vcf.file,
                covariate.file=pheno.file,
                id.column="ID_2",
                sample.ids=sample.ids,
                time.to.event="time",
                event="event",
                covariates=c("age", "SexFemale", "DrugTxYes"),
                inter.term="DrugTxYes",
                print.covs="only",
                out.file="michigan_intx_only",
                info.filter=0.3,
                maf.filter=0.005,
                chunk.size=500,
                verbose=FALSE,
                clusterObj=NULL)
```




Here we load the data and glimpse the first few values in each column outputted from the SNP*interaction survival analyses using `print.covs="only"`


```r
michigan_intx_only <- read_tsv("michigan_intx_only.coxph")
```




```r
michigan_intx_only %>% 
    head() %>%
    glimpse()
```

```
## Observations: 6
## Variables: 21
## $ RSID          <chr> "rs201487625", "rs28881575", "rs77206417", "rs13...
## $ TYPED         <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
## $ CHR           <int> 14, 14, 14, 14, 14, 14
## $ POS           <int> 19258506, 19264217, 19264589, 19264800, 19264875...
## $ REF           <chr> "A", "A", "C", "C", "G", "T"
## $ ALT           <chr> "G", "G", "T", "T", "A", "A"
## $ AF            <dbl> 0.7533720, 0.9868190, 0.0463043, 0.0105328, 0.31...
## $ MAF           <dbl> 0.2466280, 0.0131810, 0.0463043, 0.0105328, 0.31...
## $ SAMP_FREQ_ALT <dbl> 0.7861, 0.9822, 0.0363, 0.0159, 0.3101, 0.0199
## $ SAMP_MAF      <dbl> 0.2139, 0.0178, 0.0363, 0.0159, 0.3101, 0.0199
## $ INFO          <dbl> 0.925998, 0.647238, 0.487405, 0.378574, 0.583201...
## $ ER2           <dbl> NA, NA, NA, NA, NA, NA
## $ PVALUE        <dbl> 0.45928864, 0.02343794, 0.29679339, 0.31659650, ...
## $ HR            <dbl> 1.331283e+00, 4.124811e+01, 3.413024e-01, 7.0027...
## $ HR_lowerCI    <dbl> 6.239217e-01, 1.653024e+00, 4.530005e-02, 2.3043...
## $ HR_upperCI    <dbl> 2.840604e+00, 1.029269e+03, 2.571462e+00, 2.1281...
## $ COEF          <dbl> 0.2861432, 3.7196052, -1.0749863, 11.1566440, -0...
## $ SE.COEF       <dbl> 0.3866702, 1.6413260, 1.0303372, 11.1401955, 0.4...
## $ Z             <dbl> 0.7400187, 2.2662196, -1.0433345, 1.0014765, -0....
## $ N             <int> 253, 253, 253, 253, 253, 253
## $ NEVENT        <int> 100, 100, 100, 100, 100, 100
```


## Sanger Imputation Server
[Sanger Imputation Server](https://imputation.sanger.ac.uk/) pre-phases typed genotypes using either SHAPEIT or EAGLE, imputes genotypes using PBWT algorithm and outputs a `.vcf.gz` file for each chromosome. These `.vcf.gz` files are used as input for `gwasurvivr`. The function, `sangerCoxSurv`  uses a modification of cox proportional hazard regression from the R library `survival`. Built specifically for genetic data, `sangerCoxSurv` allows the user to filter on info score (imputation quality metric) and minor allele frequency from the reference panel used for imputation using `RefPanelAF` as the input arguement for `maf.filter`. Users are also provided with the sample minor allele frequency in the `sangerCoxSurv` output. 

Samples can be selected for analyses by providing a vector of `sample.ids`. The output from Sanger imputation server returns the samples as `SAMP1, ..., SAMPN`, where `N` is the total number of samples. The sample order corresponds to the sample order in the vcf file you used for imputation. Note, sample order can also be found in the `.fam` file if genotyping data were initially in `.bed`, `.bim` and `.fam` (PLINK) format prior to conversion to VCF. If no sample list is specified all samples are included in the analyses.


```r
vcf.file <- system.file(package="gwasurvivr",
                        "extdata", 
                        "sanger.pbwt_reference_impute.vcf.gz")
pheno.fl <- system.file(package="gwasurvivr",
                        "extdata", 
                        "simulated_pheno.txt")
pheno.file <- read_delim(pheno.fl,
                         delim=" ",
                         col_names=TRUE)
```

```
## Parsed with column specification:
## cols(
##   ID_1 = col_integer(),
##   ID_2 = col_character(),
##   event = col_integer(),
##   time = col_double(),
##   age = col_double(),
##   DrugTxYes = col_integer(),
##   sex = col_character(),
##   group = col_character()
## )
```

```r
pheno.file %>%
    head()
```

```
## # A tibble: 6 x 8
##    ID_1 ID_2  event  time   age DrugTxYes sex    group       
##   <int> <chr> <int> <dbl> <dbl>     <int> <chr>  <chr>       
## 1     1 SAMP1     0 12     33.9         0 male   control     
## 2     2 SAMP2     1  7.61  58.7         1 male   experimental
## 3     3 SAMP3     0 12     39.4         0 female control     
## 4     4 SAMP4     0  4.3   38.8         0 male   control     
## 5     5 SAMP5     0 12     43.6         0 male   experimental
## 6     6 SAMP6     1  2.6   57.7         0 male   control
```

In this example, we will select samples from the `experimental` group and will test survival only on these patients. The first column in the `pheno.file` are sample IDs (we will match on these). We include `age`, `DrugTxYes`, and `sex` in the survival model as covariates. 


```r
# create new sex variable where Females are 1 and males are 0
pheno.file <- pheno.file %>%
        mutate(SexFemale=if_else(sex=="female", 1L, 0L))

# select only experimental group sample.ids
sample.ids <- pheno.file %>%
        filter(group=="experimental") %$%
        ID_2 
sample.ids %>%
    head()
```

```
## [1] "SAMP2"  "SAMP5"  "SAMP7"  "SAMP9"  "SAMP11" "SAMP12"
```

We perform the analysis using the `experimental` group to demonstrate how one may want to prepare their data if not initially all samples are patients or cases (i.e. a case-control study and survival of cases is of interest). We also are showing how the IDs (`sample.ids`) need to be a vector of class `character`. The `chunk.size` refers to size of each data chunk read in and is defaulted to 10,000 rows, users can customize that to their needs. The larger the `chunk.size` the more memory (RAM) required to run the analysis. The recommended `chunk.size=500` and probably should not exceed `chunk.size=5000`. 

By default survival estimates and pvalues for the SNP adjusted for other covariates are outputted (`print.covs='only'`), however users can select `print.covs=all` to get the coefficient estimates for covariates included in the model. Depending on the number of covariates included this can add substantially to output file size. 
Next we run `sangerCoxSurv` with the default, `print.covs="only"`, load the results into R and provide descriptions of output by column. We will then run the analysis again using `print.covs="all"`. `verbose=TRUE` is used for these examples so the function display messages while running.
>>>>>>> 241dd9d56fcd1de0d18a0202696e6f658748d977


For R >= 3.5:  
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("gwasurvivr", version = "devel")
```

Alternatively:  
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("suchestoncampbelllab/gwasurvivr")

```

For R >= 3.4 and R < 3.5:  
```
source("https://bioconductor.org/biocLite.R")
biocLite("gwasurvivr")
```

# How to use package
Please refer to the [vignette](http://bioconductor.org/packages/devel/bioc/vignettes/gwasurvivr/inst/doc/gwasurvivr_Introduction.html) for a detailed description on how to use gwasurvivr functions for survival analysis (Cox proportional hazard model).



