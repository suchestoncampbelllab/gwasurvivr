coxVcf <- function(x, data, cox.params, cl){

    covariates <- x$covariates
    maf.filter <- x$maf.filter
    # info.filter <- x$info.filter
    inter.term <- x$inter.term
    print.covs <- x$print.covs
    ####### Get genotype data ############ 
    # read dosage data from collapsed vcf, subset for defined ids
    genotypes <- geno(data)$DS[, cox.params$ids, drop=FALSE]
    ########################################
    
    ##################################################
    ######## collect the snp info into df ############
    # grab info, REFPAN_AF, TYPED/IMPUTED, INFO
    # calculates sample MAF
    snp.ids <- rownames(data)
    snp.ranges <- data.frame(rowRanges(data))[,c("seqnames",
                                                 "start",
                                                 "REF", 
                                                 "ALT")]
    snp.ranges$ALT <- sapply(snp.ranges$ALT, as.character)
    
    snp.ranges <- addSnpRangesVectors(x, snp.ranges)
    
    snp.meta <- data.frame(info(data))
    
    snp.meta <- addSnpMetaVectors(x, snp.meta)
    
    samp.exp_alt <- round(rowMeans2(genotypes)*0.5, 4)
    samp.maf <- ifelse(samp.exp_alt > 0.5, 1-samp.exp_alt, samp.exp_alt)
    snp <- cbind(RSID=snp.ids,
                 snp.ranges,
                 snp.meta,
                 SAMP_FREQ_ALT=samp.exp_alt,
                 SAMP_MAF=samp.maf)
    ##################################################
    
    ##################################################
    ##### SNP Checks #################################
    # remove snps with SD less than 1e-4
    # to put this in perspective:
    # a sample size of 100 000 000 with only 1 person being 1 and rest 0,
    # has an SD = 1e-4
    # x <- c(rep(0, 1e8),1)
    # sd(x)
    snp.keep <- rowSds(genotypes) > 1e-4
    if(!all(snp.keep)){
        genotypes <- genotypes[snp.keep,]
        snp.drop <- snp[!snp.keep,]
        snp <- snp[snp.keep,]
    }else{
        snp.drop <- data.frame()
    }
    
    
    empty.geno <- tryCatch(
        {

            snpRef <- getSnpRef(x, snp)
            x.filter <- getFilter(x)
            threshold_name <- getThresholdName(x)
            
            # Further filter by user defined thresholds
            if(!is.null(maf.filter)){
                
                ok.maf <- snpRef>maf.filter & snpRef<(1-maf.filter)
                snp.drop <- base::rbind(snp.drop,snp[!ok.maf,])
                snp <- snp[ok.maf,]
                if(all(!ok.maf)) stop("None of the SNPs pass the MAF threshold")
                genotypes <- genotypes[ok.maf,]
            }
            
            if(!is.null(x.filter)){
                ok.info <- getOkInfo(x, snp, x.filter)
                snp.drop <- base::rbind(snp.drop,snp[!ok.info,])
                snp <- snp[ok.info,]
                if(all(!ok.info)) {
                    stop(paste0("None of the SNPs pass the ",threshold_name, " threshold"))
                }
                genotypes <- genotypes[ok.info,]
            }
            #############################################################
            ########## clean and save dropped and kept SNP info #########
            # rearrange columns for snp info
            # snp.cols <- c("RSID",
            #               "CHR", 
            #               "POS", 
            #               "REF",
            #               "ALT",
            #               "RefPanelAF", 
            #               "TYPED", 
            #               "INFO",
            #               "SAMP_FREQ_ALT",
            #               "SAMP_MAF")
            colnames(snp) <- x$snp.cols
            colnames(snp.drop) <- x$snp.cols
            # snp.ord <- c("RSID",
            #              "TYPED", 
            #              "CHR",
            #              "POS",
            #              "REF",
            #              "ALT", 
            #              "RefPanelAF",
            #              "SAMP_FREQ_ALT",
            #              "SAMP_MAF",
            #              "INFO")
            snp <- snp[, x$snp.ord]
            snp.drop <- snp.drop[, x$snp.ord]
            ###########################################################
            
            ###########################################################
            ############### fit models in parallel ####################
            cox.out <- getGenotypesCoxOut(inter.term, genotypes, cl, cox.params,
                                 print.covs)
            #############################################
            mod.out <- list(dropped.snps=snp.drop)
            mod.out$res <- coxExtract(cox.out,
                                         snp, 
                                         print.covs)
            return(mod.out)
        },
        error=function(err) err
    )
    if(inherits(empty.geno, "error")){
        mod.out <- list(dropped.snps=snp.drop)
        mod.out$res <- NULL
        return(mod.out)
    } 
}
