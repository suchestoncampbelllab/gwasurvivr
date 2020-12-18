createGdsCoxSurv <- function(gdsfile,
                             covariate.file,
                             id.column,
                             sample.ids,
                             time.to.event,
                             event,
                             covariates,
                             inter.term,
                             print.covs,
                             out.file,
                             chunk.size,
                             maf.filter,
                             exclude.snps,
                             flip.dosage,
                             verbose,
                             clusterObj){
  
  cox_surv <- list(gdsfile = gdsfile,
                   covariate.file = covariate.file,
                   id.column = id.column,
                   sample.ids = sample.ids,
                   time.to.event = time.to.event, 
                   event = event,
                   covariates = covariates,
                   inter.term = inter.term,
                   print.covs = print.covs,
                   out.file = out.file,
                   chunk.size = chunk.size,
                   maf.filter = maf.filter,
                   exclude.snps = exclude.snps,
                   flip.dosage = flip.dosage,
                   verbose = verbose,
                   clusterObj = clusterObj,
                   columnHeadings = c("RSID", "TYPED", "CHR", "POS", "A0","A1", "exp_freq_A1", 
                                      "SAMP_MAF"),
                   snp.df = data.frame(matrix(data = rep(NA, 16), ncol = 8 )),
                   snp.cols = c("snpID","TYPED","RSID","POS","A0","A1", "CHR"),
                   snp.ord = c("RSID","TYPED", "CHR","POS","A0","A1"))
  
  class(cox_surv) <- c("GdsCoxSurv", "PlinkGdsImpute2CoxSurv")
  
  return(cox_surv)
}

createSnpSpike.GdsCoxSurv <- function(x, cox.params){
  rbind(c(rnorm(nrow(cox.params$pheno.file)-3), rep(NA, 3)),
        c(rnorm(nrow(cox.params$pheno.file)-4), rep(NA, 4)))
}


processSNPGenotypes.GdsCoxSurv <- function(x, snp, genotypes, scanAnn, 
                                           exclude.snps = NULL, 
                                           cox.params, verbose) {
  # assign rsIDs (pasted with imputation status) as rows
  # and sample ID as columns to genotype file
  dimnames(genotypes) <- list(paste(snp$snp, snp$rsID, sep=";"),
                              scanAnn$ID_2)
  
  # Subset genotypes by given samples
  genotypes <- genotypes[,cox.params$ids]
  
  return(list(snp = snp, genotypes = genotypes))
}



getGenoData.GdsCoxSurv <- function(x){
  
  # read genotype
  ## need to add if statement about dimensions
  # set default "snp,scan" -- 
  # in GWASTools documentation say it needs to be in this orientation
  gds <- GdsGenotypeReader(x$gdsfile, genotypeDim="scan,snp")
  # close gds file on exit of the function
  # on.exit(close(gds), add=TRUE)
  # aux files
  snpfile <- replaceFileExt(file.path = x$gdsfile, ext = ".snp.rdata")
  scanfile <- replaceFileExt(file.path = x$gdsfile, ext = ".scan.rdata")
  # read in snp data
  snpAnnot <- getobj(snpfile)
  # read scan
  scanAnnot <- getobj(scanfile)
  # put into GenotypeData coding 
  genoData <- GenotypeData(gds,
                           snpAnnot=snpAnnot,
                           scanAnnot=scanAnnot)
  
  return(genoData)
}