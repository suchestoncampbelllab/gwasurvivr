runOnChunks <- function(genoData, chunk.size, verbose, 
                        cox.params, flip.dosage, exclude.snps, maf.filter, inter.term,
                        cl, print.covs, out.file, snp.cols, snp.ord,
                        funProcessSNPGenotypes) {

  # number of snps in segment
  snp.start <- 1
  snp.end <- nsnp(genoData)
  # number of dropped snps
  snp.drop.n <-0
  snp.n <- 0
  
  # get genotypes for certain chunk size
  nsnp.seg <- snp.end - snp.start + 1
  nchunks <- ceiling(nsnp.seg/chunk.size)
  
  for(i in seq_len(nchunks)){
    
    if(verbose) message("Analyzing part ", i, "/", nchunks, "...")
    
    # set up chunks
    next.chunk <- (i-1)*chunk.size
    next.chunk.start <- snp.start + next.chunk
    snp.chunk <- ifelse(next.chunk.start + chunk.size > snp.end,
                        snp.end - next.chunk.start + 1,
                        chunk.size)
    chunk.idx <- (next.chunk+1):(next.chunk+snp.chunk)
    
    # get genotypes for chunk
    genotypes <- getGenotype(genoData,
                             snp=c(next.chunk.start, snp.chunk),
                             scan=c(1,-1),
                             drop=FALSE)
    
    # get the snp info file
    snp <- getAnnotation(getSnpAnnotation(genoData))[chunk.idx,]
    # snp.cols <- c("snpID",
    #               "TYPED",
    #               "RSID",
    #               "POS",
    #               "A0",
    #               "A1",
    #               "CHR")
    # snp.ord <- c("RSID",
    #              "TYPED",
    #              "CHR",
    #              "POS",
    #              "A0",
    #              "A1")

    colnames(snp) <- snp.cols
    snp <- snp[, snp.ord]
    
    # grab sample file data
    scanAnn <- getAnnotation(getScanAnnotation(genoData))
    
    # assign rsIDs (pasted with imputation status) as rows 
    # and sample ID as columns to genotype file
    # dimnames(genotypes) <- list(paste(snp$TYPED, snp$RSID, sep=";"), 
    #                             scanAnn$ID_2)
    # 
    # # Subset genotypes by given samples
    # 
    # if(is.null(exclude.snps)){
    #   genotypes <- genotypes[,cox.params$ids]
    # } else {
    #   genotypes <- genotypes[!snp$RSID %in% exclude.snps,cox.params$ids]
    #   snp <- snp[! snp$RSID %in% exclude.snps, ]
    #   if(verbose) message(length(exclude.snps), " SNPs are excluded based on the exclusion list provided by the user")
    # }
    # flip dosage
    listSNPGenotype <- funProcessSNPGenotypes(
      snp = snp, genotypes = genotypes, scanAnn = scanAnn, 
      exclude.snps = exclude.snps, cox.params = cox.params, verbose = verbose)
    
    snp <- listSNPGenotype$snp
    genotypes <- listSNPGenotype$genotypes
    
    if(flip.dosage) genotypes <- 2 - genotypes
    ########################################################################
    
    ###############################################################
    ##### SNP filtering ###########################################
    
    ### Check snps for MAF = 0  ###
    # remove snps with SD less than 1e-4
    # to put this in perspective:
    # a sample size of 100 000 000 with only 1 person being 1 and rest 0,
    # has an SD = 1e-4
    # x <- c(rep(0, 1e8),1)
    # sd(x)
    ok.snp <- rowSds(genotypes, na.rm = TRUE) > 1e-4
    snp <- snp[ok.snp, ]
    genotypes <- genotypes[ok.snp, ]
    snp.drop <- snp[!ok.snp, ]
    
    # calculate MAF
    snp$exp_freq_A1 <- round(rowMeans2(genotypes, na.rm = TRUE)*0.5,4)
    snp$SAMP_MAF <- ifelse(snp$exp_freq_A1 > 0.5,
                           1-snp$exp_freq_A1,
                           snp$exp_freq_A1
    )
    
    # Further filter by user defined thresholds
    if (!is.null(maf.filter)) {
      ok.snp <- snp$SAMP_MAF > maf.filter
      genotypes <- genotypes[ok.snp,]
      snp <- snp[ok.snp,]
      
      if(nrow(snp.drop) > 0){
        snp.drop$exp_freq_A1 <- 1
        snp.drop$SAMP_MAF <- 0
        snp.drop <- rbind(snp.drop, snp[!ok.snp,])
      } else {
        snp.drop <- snp[!ok.snp,]
      }
    }
    
    if (nrow(snp.drop) > 0) {
      write.table(
        snp.drop,
        paste0(out.file, ".snps_removed"),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE,
        append = TRUE )
      snp.drop.n <- snp.drop.n+nrow(snp.drop)
    }
    
    if (nrow(genotypes) > 0) {
      # fit models in parallel
      cox.out <- getGenotypesCoxOut(inter.term, genotypes, cl, cox.params,
                                    print.covs)
      
      res <- coxExtract(cox.out,
                        snp,
                        print.covs)
      
      write.table(
        res,
        file = paste0(out.file, ".coxph"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        append = TRUE
      )
      
      snp.n <- nrow(genotypes) + snp.n
      
    }
    
    
    
    
    return(list(snp.drop.n = snp.drop.n, 
                snp.n = snp.n))
    
  }
}

