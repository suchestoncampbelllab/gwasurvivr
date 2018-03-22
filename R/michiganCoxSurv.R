#' Fit cox survival to all variants in a .vcf.gz file from Sanger imputation server
#' 
#' Performs survival analysis using Cox proportional hazard models on imputed genetic data stored in compressed VCF files 
#' 
#' @param vcf.file character(1) path to VCF file.
#' @param vcf.file character(1) path to corresponding info file.
#' @param pheno.file matrix(1) comprising phenotype data. 
#' @param time character(1) string that matches time column name in pheno.file
#' @param event character(1) string that matches event column name in pheno.file
#' @param covariates character vector with matching column names in pheno.file of covariates of interest
#' @param sample.ids character vector with sample ids to include in analysis
#' @param chunk.size integer(1) number of variants to process per thread
#' @param info.filter integer(1) of imputation quality score filter (i.e. 0.7 will filter info > 0.7)
#' @param maf.filter integer(1) filter out minor allele frequency below threshold (i.e. 0.005 will filter MAF > 0.005)
#' @param output.name character(1) string with output name
#' 
#' @return 
#' Saves text file directly to disk that contains survival analysis results
#' 
#' @examples 
#' vcf.file <- system.file(package="gwasurvivr","extdata", "sanger.pbwt_reference_impute.vcf.gz")
#' pheno.fl <- system.file(package="gwasurvivr", "extdata", "simulated_pheno.txt")
#' pheno.file <- read.table(pheno.fl, sep=" ", header=TRUE, stringsAsFactors = FALSE)
#' library(tidyverse)
#' library(magrittr)
#' pheno.file <- pheno.file %>% mutate(SexFemale=case_when(sex=="female"~1L, sex=="male"~0L)) %>% dplyr::select(-ID_1)
#' sample.ids <- pheno.file %>% filter(group=="experimental") %$% ID_2 
#' michiganCoxSurv(vcf.file=vcf.file, pheno.file=pheno.file, time.to.event="time", event="event", covariates=c("age", "SexFemale", "bmiOVWT"), sample.ids=sample.ids, output.name="sanger_example", chunk.size=10000, info.filter=0.7, maf.filter=0.005, verbose=TRUE)
#' @importFrom survival Surv coxph.fit
#' @import parallel
#' @import VariantAnnotation
#' 
#' @export

michiganCoxSurv <- function(vcf.file,
                            pheno.file,
                            time.to.event,
                            event,
                            covariates,
                            sample.ids,
                            output.name,
                            chunk.size,
                            info.filter,
                            maf.filter,
                            verbose=TRUE){
        
}