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

# Installation
This package is currently available on [Bioconductor devel branch](https://bioconductor.org/packages/devel/bioc/html/gwasurvivr.html) or by using `devtools` library for R >= 3.4 and going to the Sucheston Campbell Lab GitHub repository (this page).  If using R 3.5, use `BiocManager` to install the package, if using R >= 3.4, `BiocInstaller` or `biocLite` can be used.


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



