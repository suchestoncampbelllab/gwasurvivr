loadProcessWrite.MichiganSangerCoxSurv<- function(x, cl, cox.params) {
  
  ################################################
  #### open VCF file connection ##################
  vcf <- VcfFile(x$vcf.file, yieldSize=x$chunk.size)
  open(vcf)
  
  ################################################
  ####### read first chunk #######################
  chunk.start <- 0
  # Initialize so an empty first chunk still returns valid counts (#25).
  snps_removed <- 0
  snps_analyzed <- 0

  
  ################################################
  ##### Start repeat loop ########################
  # get genotype probabilities by chunks
  # apply the survival function and save output
  repeat{ 
    # read in just dosage data from Vcf file
    if(x$verbose) message("Analyzing chunk ",
                          chunk.start,
                          "-",
                          chunk.start+x$chunk.size)    
    
    data <- readVcf(vcf,
                    param=ScanVcfParam(geno=x$geno.field,
                                       info=x$ScanVcfParamInfo))
    
    if(nrow(data)==0) break
    
    
    out.list <- coxVcf(x, data, cox.params, cl)
    
    if(chunk.start == 0) {
      
      write.table(
        out.list$res,
        paste0(x$out.file, ".coxph"),
        append = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE,
        sep = "\t"
      )
      write.table(
        out.list$dropped.snps,
        paste0(x$out.file, ".snps_removed"),
        append = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        sep = "\t"
      )
      
      chunk.start <- x$chunk.size
      snps_removed <- nrow(out.list$dropped.snps)
      snps_analyzed <- nrow(out.list$res)
      
    } else {
      
      
      write.table(
        out.list$res,
        paste0(x$out.file, ".coxph"),
        append = TRUE,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        sep = "\t"
      )
      write.table(
        out.list$dropped.snps,
        paste0(x$out.file, ".snps_removed"),
        append = TRUE,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        sep = "\t"
      )
      chunk.start <- chunk.start + x$chunk.size
      snps_removed <- snps_removed + nrow(out.list$dropped.snps)
      snps_analyzed <-  snps_analyzed + nrow(out.list$res)
    }
  }
  
  ################################################
  close(vcf)

  # Named list: coxSurv() / closing_messages() read these by name. An unnamed
  # list left `snps_analyzed` NULL -> blank "SNPs were analyzed" message (#25).
  return(list(snps_removed = snps_removed, snps_analyzed = snps_analyzed))

}

coxVcf <- function(x, data, cox.params, cl){
  
  covariates <- x$covariates
  maf.filter <- x$maf.filter
  # info.filter <- x$info.filter
  inter.term <- x$inter.term
  print.covs <- x$print.covs
  ####### Get genotype data ############
  # Convert the requested VCF genotype field (DS/GT/HDS) to ALT dosage, then
  # subset for defined ids.
  genotypes <- genotypeFieldToDosage(geno(data)[[x$geno.field]],
                                     x$geno.field, verbose = x$verbose)
  genotypes <- genotypes[, cox.params$ids, drop=FALSE]
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
  #
  # Every filter below is NA-safe (NA -> drop) and every subset uses
  # drop=FALSE, so each variant in the chunk lands in exactly one of
  # .coxph / .snps_removed (#42). Variants with any missing dosage in the
  # analysed samples are removed: the Cox fit cannot use them.
  colnames(snp) <- x$snp.cols
  keep <- !rowAnyNAs(genotypes)
  keep[keep] <- rowSds(genotypes[keep, , drop=FALSE]) > 1e-4

  # Further filter by user defined thresholds
  if(!is.null(maf.filter)){
    snpRef <- getSnpRef(x, snp)
    keep <- keep & !is.na(snpRef) &
      snpRef > maf.filter & snpRef < (1-maf.filter)
  }
  x.filter <- getFilter(x)
  if(!is.null(x.filter)){
    ok.info <- getOkInfo(x, snp, x.filter)
    keep <- keep & !is.na(ok.info) & ok.info
  }

  #############################################################
  ########## clean and save dropped and kept SNP info #########
  snp <- snp[, x$snp.ord, drop=FALSE]
  # Effect allele is unambiguous for VCF dosages: DS/GT/HDS count the ALT
  # allele, so the hazard ratio is per 1 copy of ALT. Emit it explicitly so
  # users never have to infer which allele the HR refers to.
  snp$EFFECT_ALLELE <- snp$ALT
  snp$OTHER_ALLELE  <- snp$REF

  snp.drop <- snp[!keep, , drop=FALSE]
  snp <- snp[keep, , drop=FALSE]
  genotypes <- genotypes[keep, , drop=FALSE]

  mod.out <- list(dropped.snps=snp.drop, res=NULL)
  if(nrow(genotypes) == 0) return(mod.out)

  ###########################################################
  ############### fit models in parallel ####################
  # Errors are NOT swallowed here: previously any failure silently discarded
  # every variant in the chunk from both output files (#42).
  cox.out <- tryCatch(
    getGenotypesCoxOut(inter.term, genotypes, cl, cox.params, print.covs),
    error=function(err) {
      stop("Cox fit failed for variants ", snp$RSID[1], " .. ",
           snp$RSID[nrow(snp)], ": ", conditionMessage(err), call.=FALSE)
    })
  mod.out$res <- coxExtract(cox.out, snp, print.covs)
  mod.out
}
