createSangerCoxSurv <- function(vcf.file,
                                  covariate.file,
                                  id.column,
                                  sample.ids,
                                  time.to.event,
                                  event,
                                  covariates,
                                  inter.term,
                                  print.covs,
                                  out.file,
                                  maf.filter,
                                info.filter,
                                  chunk.size,
                                  verbose,
                                  clusterObj){
  
  cox_surv <- list(vcf.file = vcf.file,
                   covariate.file = covariate.file,
                   id.column = id.column,
                   sample.ids=sample.ids,
                   time.to.event = time.to.event,
                   event = event,
                   covariates = covariates,
                   inter.term=inter.term,
                   print.covs=print.covs,
                   out.file = out.file,
                   maf.filter=maf.filter,
                   info.filter=info.filter,
                   chunk.size=chunk.size,
                   verbose=verbose,
                   clusterObj=clusterObj,
                   ScanVcfParamInfo = c("RefPanelAF", "TYPED", "INFO"),
                   snp.cols = c("RSID",
                                 "CHR", 
                                 "POS", 
                                 "REF",
                                 "ALT",
                                 "RefPanelAF", 
                                 "TYPED", 
                                 "INFO",
                                 "SAMP_FREQ_ALT",
                                 "SAMP_MAF"),
                   snp.ord = c("RSID",
                                "TYPED", 
                                "CHR",
                                "POS",
                                "REF",
                                "ALT", 
                                "RefPanelAF",
                                "SAMP_FREQ_ALT",
                                "SAMP_MAF",
                                "INFO")
                   )
  
  class(cox_surv) <- "SangerCoxSurv"
  
  return(cox_surv)
}


loadProcessWrite.SangerCoxSurv <- function(x,
                                             cl,
                                             cox.params) {
  
  
  #### open VCF file connection ##################
  vcf <- VcfFile(x$vcf.file, yieldSize=x$chunk.size)
  open(vcf)
  
  ################################################
  ####### read first chunk #######################
  chunk.start <- 0
  
  
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
                    param=ScanVcfParam(geno="DS", 
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
      
      chunk.start <- chunk.start+x$chunk.size
      snps_removed <- snps_removed+nrow(out.list$dropped.snps)
      snps_analyzed <-  snps_analyzed+nrow(out.list$res)
    }
    
  }
  ################################################
  
  close(vcf)
  
  return(list(snps_removed, snps_analyzed))
  
}


addSnpRangesVectors.SangerCoxSurv <- function(x, snp.ranges) {
  return(snp.ranges)
}

addSnpMetaVectors.SangerCoxSurv <- function(x, snp.meta){
  snp.meta$RefPanelAF <- sapply(snp.meta$RefPanelAF, as.numeric)
  return(snp.meta)
}

getSnpRef.SangerCoxSurv <- function(x, snp){
  snp$RefPanelAF
}
getFilter.SangerCoxSurv <- function(x){
  x$info.filter
}
getThresholdName.SangerCoxSurv <- function(x){
  return("info")
}
getOkInfo.SangerCoxSurv <- function(x, snp, x.filter){
  snp$INFO >= x.filter
}