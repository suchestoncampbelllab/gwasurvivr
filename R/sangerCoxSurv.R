#' Fit cox survival to all variants in a .vcf.gz file from Sanger imputation server
#' 
#' Performs survival analysis using Cox proportional hazard models on imputed genetic data stored in compressed VCF files 
#' 
#' @param vcf.file character(1) path to VCF file.
#' @param pheno.file matrix(1) comprising phenotype data. 
#' @param time.to.event character(1) string that matches time column name in pheno.file
#' @param event character(1) string that matches event column name in pheno.file
#' @param covariates character vector with matching column names in pheno.file of covariates of interest
#' @param sample.ids character vector with sample ids to include in analysis
#' @param chunk.size integer(1) number of variants to process per thread
#' @param info.filter integer(1) of imputation quality score filter (i.e. 0.7 will filter info > 0.7)
#' @param maf.filter integer(1) filter out minor allele frequency below threshold (i.e. 0.005 will filter MAF > 0.005)
#' @param out.file character(1) string with output name
#' @param verbose logical(1) for messages that describe which part of the analysis is currently being run
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
#' pheno.file <- pheno.file %>%  
#'                     mutate(SexFemale=case_when(
#'                                       sex=="female"~1L,
#'                                       sex=="male"~0L)
#'                           ) %>% 
#'                     select(-ID_1)
#' sample.ids <- pheno.file %>%
#'                     filter(group=="experimental") %$%
#'                     ID_2 
#' sangerCoxSurv(vcf.file=vcf.file,
#'               pheno.file=pheno.file,
#'               time.to.event="time",
#'               event="event",
#'               covariates=c("age", "SexFemale", "bmiOVWT"),
#'               sample.ids=sample.ids,
#'               out.file="sanger_example",
#'               chunk.size=10000,
#'               info.filter=0.7,
#'               maf.filter=0.005,
#'               verbose=TRUE)
#' 
#' @importFrom survival Surv coxph.fit
#' @importFrom utils write.table
#' @importFrom matrixStats rowMeans2
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom stats pnorm
#' @import parallel
#' @import VariantAnnotation
#' 
#' @export

sangerCoxSurv <- function(vcf.file,
                          covariate.file,
                          id.column,
                          sample.ids=NULL, 
                          time.to.event, 
                          event,
                          covariates,
                          inter.term=NULL,
                          print.covs="only",
                          out.file,
                          maf.filter=0.05,
                          info.filter=NULL,
                          chunk.size=5000,
                          verbose=TRUE,
                          clusterObj=NULL
){
    ################################################
    # #### Phenotype data wrangling ################
 
    cox.params <- coxPheno(covariate.file, covariates, id.column,inter.term, time.to.event, event, verbose)
    ################################################
    
    ####################################################
    ########## Cluster object ########################
    # create cluster object depending on user pref or OS type,
    # also create option to input number of cores
    if(!is.null(clusterObj)){
        cl <- clusterObj
    }else if(.Platform$OS.type == "unix") {
        cl <- makeForkCluster(getOption("gwasurvivr.cores", detectCores()))
    } else {
        cl <- makeCluster(getOption("gwasurvivr.cores", detectCores()))
    }
    on.exit(stopCluster(cl), add=TRUE)
    #################################################################
    
    #### open VCF file connection ######
    vcf <- VcfFile(vcf.file, yieldSize=chunk.size)
    open(vcf)
    
    ##################################################
    ####### read first chunk #########################
    chunk.start <- 0
    if(verbose) message("Analyzing chunk ", chunk.start, "-", chunk.start+chunk.size)    
    
    data <- readVcf(vcf, param=ScanVcfParam(geno="DS", info=c("RefPanelAF", "TYPED", "INFO")))
    out.list <- coxVcfSanger(data, cox.params, cl, inter.term, print.covs)
    write.table(
        out.list$res,
        paste0(out.file, ".coxph"),
        append = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )
    write.table(
        out.list$dropped.snps,
        paste0(out.file, ".snps_removed"),
        append = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )
    
    
    chunk.start <- chunk.size
    snps_removed <- nrow(out.list$dropped.snps)
    
    ##### Start repeat loop #####
    # # get genotype probabilities by chunks
    # # apply the survival function and save output
    
    repeat{ 
        # read in just dosage data from Vcf file
        if(verbose) message("Analyzing chunk ", chunk.start, "-", chunk.start+chunk.size)    
        
        data <- readVcf(vcf, param=ScanVcfParam(geno="DS", info=c("RefPanelAF", "TYPED", "INFO")))
        
        if(nrow(data)==0){
            break
        }
        
        out.list <- coxVcfSanger(data, cox.params, cl, inter.term, print.covs)
        write.table(
            out.list$res,
            paste0(out.file, ".coxph"),
            append = TRUE,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = "\t"
        )
        write.table(
            out.list$dropped.snps,
            paste0(out.file, ".snps_removed"),
            append = TRUE,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = "\t"
        )
        
        
        chunk.start <- chunk.start+chunk.size
        snps_removed <- snps_removed+nrow(out.list$dropped.snps)
        
        
    }
    close(vcf)
    if(verbose) message("Analysis completed on ", format(Sys.time(), "%Y-%m-%d"), " at ", format(Sys.time(), "%H:%M:%S"))
    if(verbose) message(snps_removed, " SNPs were removed from the analysis for not meeting the threshold criteria.")
    if(verbose) message("List of removed SNPs can be found in ", paste0(out.file, ".snps_removed"))
}
