
createMichiganCoxSurvModified <- function(vcf.file,
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
                   chunk.size=chunk.size,
                   verbose=verbose,
                   clusterObj=clusterObj,
                   ScanVcfParamInfo = c("PR"),
                   snp.cols = c("RSID",
                                 "CHR",
                                 "POS", 
                                 "REF", 
                                 "ALT", 
                                 "PR",
                                 "SAMP_MAF"),
                   snp.ord = c("RSID",
                                "CHR", 
                                "POS",
                                "REF", 
                                "ALT", 
                                "PR",
                                "SAMP_MAF"))
  
  class(cox_surv) <- c("MichiganCoxSurvModified", "MichiganSangerCoxSurvModified")
  
  return(cox_surv)
}

addSnpMetaVectors.MichiganCoxSurvModified <- function(x, snp.meta){
  snp.meta$TYPED <- TRUE
  return(snp.meta)
}

addSnpMetaVectors <- function(x, snp.meta){
    snp.meta$TYPED <- TRUE
    return(snp.meta)
}


addSnpRangesVectors.MichiganCoxSurvModified <- function(x, snp.ranges) {
  snp.ranges$REF <- sapply(snp.ranges$REF, as.character)
  return(snp.ranges)
}

getSnpRef.MichiganCoxSurvModified <- function(x, snp){
  snp$MAF
}
getFilter.MichiganCoxSurvModified <- function(x, snp){
  x$r2.filter
}

getThresholdName.MichiganCoxSurvModified <- function(x){
  return("R2")
}

getOkInfo.MichiganCoxSurvModified <- function(x, snp, x.filter){
  !is.na(snp$R2 >= x.filter)
}