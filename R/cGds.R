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
                   clusterObj = clusterObj)
  
  class(cox_surv) <- "GdsCoxSurv"
  
  return(cox_surv)
}


loadProcessWrite.GdsCoxSurv <- function(x,
                                            cl,
                                            cox.params) {
  
  writeFileHeadings(
    cols = c("RSID", "TYPED", "CHR", "POS", "A0","A1", "exp_freq_A1", 
             "SAMP_MAF"), 
    out.file = x$out.file,
    inter.term = x$inter.term,
    snp.df = data.frame(matrix(data = rep(NA, 16), ncol = 8 ) ), 
    snp.spike = rbind(c(rnorm(nrow(cox.params$pheno.file)-3), rep(NA, 3)),
                      c(rnorm(nrow(cox.params$pheno.file)-4), rep(NA, 4))),
    print.covs = x$print.covs,
    cox.params = cox.params)
  
  ############################################################################
  ##### Load Genotype data ###################################################
  
  genoData <- getGenoData(x, x$gdsfile)
  
  ############################################################################
  ##### Genotype data wrangling ##############################################
  
  results <- runOnChunks(x, genoData, x$chunk.size, x$verbose, 
                         cox.params, x$flip.dosage, x$exclude.snps, 
                         x$maf.filter, x$inter.term,
                         cl, x$print.covs, x$out.file, 
                         snp.cols = c("snpID","TYPED","RSID","POS","A0","A1", "CHR"),
                         snp.ord = c("RSID","TYPED", "CHR","POS","A0","A1"))
  
  return(list(snps_removed = results$snp.drop.n, 
              snps_analyzed = results$snp.n))
  
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



getGenoData.GdsCoxSurv <- function(x, gdsfile){
  
  # read genotype
  ## need to add if statement about dimensions
  # set default "snp,scan" -- 
  # in GWASTools documentation say it needs to be in this orientation
  gds <- GdsGenotypeReader(gdsfile, genotypeDim="scan,snp")
  # close gds file on exit of the function
  # on.exit(close(gds), add=TRUE)
  # aux files
  snpfile <- replaceFileExt(file.path = gdsfile, ext = ".snp.rdata")
  scanfile <- replaceFileExt(file.path = gdsfile, ext = ".scan.rdata")
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