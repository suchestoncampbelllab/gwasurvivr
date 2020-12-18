runOnChunks <- function(x, genoData, cox.params, cl,
                        snp.cols, snp.ord) {

  # number of snps in segment
  snp.start <- 1
  snp.end <- nsnp(genoData)
  # number of dropped snps
  snp.drop.n <-0
  snp.n <- 0
  
  # get genotypes for certain chunk size
  nsnp.seg <- snp.end - snp.start + 1
  nchunks <- ceiling(nsnp.seg/x$chunk.size)
  
  for(i in seq_len(nchunks)){
    
    if(x$verbose) message("Analyzing part ", i, "/", nchunks, "...")
    
    # set up chunks
    next.chunk <- (i-1)*x$chunk.size
    next.chunk.start <- snp.start + next.chunk
    snp.chunk <- ifelse(next.chunk.start + x$chunk.size > snp.end,
                        snp.end - next.chunk.start + 1,
                        x$chunk.size)
    chunk.idx <- (next.chunk+1):(next.chunk+snp.chunk)
    
    # get genotypes for chunk
    genotypes <- getGenotype(genoData,
                             snp=c(next.chunk.start, snp.chunk),
                             scan=c(1,-1),
                             drop=FALSE)
    
    # get the snp info file
    snp <- getAnnotation(getSnpAnnotation(genoData))[chunk.idx,]


    colnames(snp) <- snp.cols
    snp <- snp[, snp.ord]
    
    # grab sample file data
    scanAnn <- getAnnotation(getScanAnnotation(genoData))
    

    listSNPGenotype <- processSNPGenotypes(x,
      snp = snp, genotypes = genotypes, scanAnn = scanAnn, 
      exclude.snps = x$exclude.snps, cox.params = cox.params, verbose = x$verbose)
    
    snp <- listSNPGenotype$snp
    genotypes <- listSNPGenotype$genotypes
    
    if(x$flip.dosage) genotypes <- 2 - genotypes
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
    if (!is.null(x$maf.filter)) {
      ok.snp <- snp$SAMP_MAF > x$maf.filter
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
        paste0(x$out.file, ".snps_removed"),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE,
        append = TRUE )
      snp.drop.n <- snp.drop.n+nrow(snp.drop)
    }
    
    if (nrow(genotypes) > 0) {
      # fit models in parallel
      cox.out <- getGenotypesCoxOut(x$inter.term, genotypes, cl, cox.params,
                                    x$print.covs)
      
      res <- coxExtract(cox.out,
                        snp,
                        x$print.covs)
      
      write.table(
        res,
        file = paste0(x$out.file, ".coxph"),
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

