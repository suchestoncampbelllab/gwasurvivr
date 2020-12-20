
createMichiganCoxSurv <- function(vcf.file,
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
                                  r2.filter,
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
                   r2.filter=r2.filter,
                   chunk.size=chunk.size,
                   verbose=verbose,
                   clusterObj=clusterObj,
                   ScanVcfParamInfo = c("AF", "MAF", "R2", "ER2"),
                   snp.cols = c("RSID",
                                 "CHR",
                                 "POS", 
                                 "REF", 
                                 "ALT", 
                                 "AF", 
                                 "MAF", 
                                 "R2",
                                 "ER2",
                                 "TYPED",
                                 "SAMP_FREQ_ALT",
                                 "SAMP_MAF"),
                   snp.ord = c("RSID",
                                "TYPED",
                                "CHR", 
                                "POS",
                                "REF", 
                                "ALT", 
                                "AF",
                                "MAF",
                                "SAMP_FREQ_ALT",
                                "SAMP_MAF",
                                "R2",
                                "ER2"))
  
  class(cox_surv) <- c("MichiganCoxSurv", "MichiganSangerCoxSurv")
  
  return(cox_surv)
}

addSnpMetaVectors.MichiganCoxSurv <- function(x, snp.meta){
  snp.meta$TYPED <- NA
  snp.meta$TYPED[!is.na(snp.meta$ER2)] <- TRUE
  snp.meta$TYPED[is.na(snp.meta$ER2)] <- FALSE
  return(snp.meta)
}

addSnpRangesVectors.MichiganCoxSurv <- function(x, snp.ranges) {
  snp.ranges$REF <- sapply(snp.ranges$REF, as.character)
  return(snp.ranges)
}

getSnpRef.MichiganCoxSurv <- function(x, snp){
  snp$MAF
}
getFilter.MichiganCoxSurv <- function(x, snp){
  x$r2.filter
}
getThresholdName.MichiganCoxSurv <- function(x){
  return("R2")
}
getOkInfo.MichiganCoxSurv <- function(x, snp, x.filter){
  !is.na(snp$R2 >= x.filter)
}