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
  
  class(cox_surv) <- c("SangerCoxSurv", "MichiganSangerCoxSurv")
  
  return(cox_surv)
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