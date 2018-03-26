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
#' sangerCoxSurv(vcf.file=vcf.file, pheno.file=pheno.file, time.to.event="time", event="event", covariates=c("age", "SexFemale", "bmiOVWT"), sample.ids=sample.ids, output.name="sanger_example", chunk.size=10000, info.filter=0.7, maf.filter=0.005, verbose=TRUE)
#' @importFrom survival Surv coxph.fit
#' @import parallel
#' @import VariantAnnotation
#' 
#' @export

sangerCoxSurv <- function(vcf.file,
                          pheno.file,
                          time.to.event, 
                          event,
                          covariates,
                          sample.ids,
                          output.name,
                          chunk.size,
                          info.filter,
                          maf.filter,
                          verbose=TRUE
                       
){
    
    if(verbose) message("Analysis started on ", format(Sys.time(), "%Y-%m-%d"), " at ", format(Sys.time(), "%H:%M:%S"))
   
    # subset phenotype file for sample ids
    pheno.file <- pheno.file[match(sample.ids, pheno.file[[1]]), ]
    if(verbose) message("Analysis running for ", nrow(pheno.file), " samples.")
    
    # covariates are defined in pheno.file
    ok.covs <- colnames(pheno.file)[colnames(pheno.file) %in% covariates]
    if(verbose) message("Covariates included in the models are: ", paste(ok.covs, collapse=", "))
    if(verbose) message("If your covariates of interest are not included in the model\nplease stop the analysis and make sure user defined covariates\nmatch the column names in the pheno.file")
    
    

    vcf <- VcfFile(vcf.file, yieldSize=chunk.size)
    open(vcf)
    
    chunk.start <- 0
    snps_removed <- 0
    
    ## Save header for the cox.surv output
    write.table(t(c("RSID",
                    "TYPED",
                    "CHR",
                    "POS",
                    "REF",
                    "ALT",
                    "RefPanelAF",
                    "SAMP_MAF",
                    "INFO",
                    "COEF",
                    "SE.COEF",
                    "HR",
                    "HR_lowerCI",
                    "HR_upperCI",
                    "Z",
                    "PVALUE",
                    "N",
                    "NEVENT"
                    )), 
    paste0(output.name, ".coxph"),
    append = FALSE, 
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep="\t")
    
    # for a single machine
    cl <- makeForkCluster(getOption("gwasurvivr.cores", detectCores()))
    on.exit(stopCluster(cl), add=TRUE)
    
    # get genotype probabilities by chunks
    # apply the survival function and save output
    pheno.file <- pheno.file[,-1]
    pheno.file <- as.matrix(pheno.file[,c(time.to.event, event, covariates)])
    params <- coxParam(pheno.file, time.to.event, event, covariates, sample.ids)
    
    N <- nrow(pheno.file)
    NEVENTS <- sum(pheno.file[,event]==1)
    
    repeat{ 
        # read in just dosage data from Vcf file
        data <- readVcf(vcf, param=ScanVcfParam(geno="DS", info=c("RefPanelAF", "TYPED", "INFO")))
        
        if(nrow(data)==0){
            break
        }
        
        # read dosage data from collapsed vcf, subset for defined ids
        genotype <- geno(data)$DS[, match(sample.ids, colnames(data)) , drop=F]
        
        
        # grab info, REFPAN_AF, TYPED/IMPUTED, INFO
        # calculates sample MAF
        snp.ids <- rownames(data)
        snp.ranges <- data.frame(SummarizedExperiment::rowRanges(data))
        snp.ranges <- snp.ranges[,c("seqnames", "start", "REF", "ALT")]
        snp.meta <- data.frame(info(data))[,c("RefPanelAF", "TYPED", "INFO")]
        #rowRanges(data)$SAMP_MAF <- round(matrixStats::rowMeans2(genotype)*0.5, 4)
        samp.maf <- round(matrixStats::rowMeans2(genotype)*0.5, 4)

        snp.info <- cbind(RSID=snp.ids,
                           snp.ranges,
                           snp.meta,
                           SAMP_MAF=samp.maf)
        
        #### filtering by MAF and INFO Score #####
        # info > 0.7
        # maf > 0.005 & maf < 0.995

        snp.maf.filt <- snp.info[snp.info$INFO > info.filter,]
        idx <- snp.maf.filt$RefPanelAF > maf.filter & snp.maf.filt$RefPanelAF < (1-maf.filter)
        snp.maf.filt <- snp.maf.filt[idx,]
        
        ## record removed SNPs by filtering
        snp_maf_removed <- rownames(genotype)[-idx]
        snps_removed <- snps_removed + length(snp_maf_removed)
        
        write.table(snp_maf_removed, 
                    file= paste0(output.name, ".MAF_INFO_removed"),
                    append = TRUE, 
                    row.names = FALSE,
                    col.names = FALSE,
                    quote = FALSE,
                    sep="\t")
        
        # clean data

        snp.maf.filt$ALT <- sapply(snp.maf.filt$ALT, as.character)
        snp.maf.filt$RefPanelAF <- sapply(snp.maf.filt$RefPanelAF, as.numeric)
        colnames(snp.maf.filt) <- c("RSID", "CHR", "POS", "REF", "ALT", "RefPanelAF", "TYPED", "INFO", "SAMP_MAF")
        snp.maf.filt <- snp.maf.filt[,c("RSID", "TYPED", "CHR", "POS", "REF", "ALT", "RefPanelAF", "SAMP_MAF", "INFO")]
        
        # now we need to filter the genotype file
        genotype <- genotype[idx,]
        
        # message user
        if(verbose) message("Analyzing chunk ", chunk.start, "-", chunk.start+chunk.size)
        
        # apply survival function
        snp.out <- t(parApply(cl=cl, X=genotype, MARGIN=1, FUN=survFit, params))
        colnames(snp.out) <- c("COEF", "SE.COEF")
        
        
        # calculate statistics
        Z <- snp.out[,1]/snp.out[,2]
        PVALUE <- 2*pnorm(abs(Z), lower.tail=F)
        HR <- exp(snp.out[,1])
        HR_lowerCI <- exp(snp.out[,1] - 1.96*snp.out[,2])
        HR_upperCI <- exp(snp.out[,1] + 1.96*snp.out[,2])
        cox.out <- cbind(snp.maf.filt,
                         snp.out,
                         HR, 
                         HR_lowerCI,
                         HR_upperCI,
                         Z,
                         PVALUE,
                         N,
                         NEVENTS)
        
        write.table(cox.out,
                    paste0(output.name, ".coxph"),
                    append = TRUE, 
                    row.names = FALSE,
                    col.names = FALSE,
                    quote = FALSE,
                    sep="\t")
        
        chunk.start <- chunk.start+chunk.size
        
    }
    close(vcf)
    if(verbose) message("Analysis completed on ", format(Sys.time(), "%Y-%m-%d"), " at ", format(Sys.time(), "%H:%M:%S"))
    if(verbose) message(length(snp_maf_removed), " SNPs were removed from the analysis for not meeting the threshold criteria.")
    if(verbose) message("List of removed SNPs can be found in ", paste0(output.name, ".MAF_INFO_removed"))
}




