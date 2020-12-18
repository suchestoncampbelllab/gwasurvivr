createImpute2CoxSurv <- function(impute.file,
                                 sample.file,
                                 chr,
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
                                 clusterObj,
                                 keepGDS){
  
  cox_surv <- list(impute.file = impute.file,
                   sample.file = sample.file,
                   chr = chr,
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
                   keepGDS = keepGDS,
                   columnHeadings = c("RSID", "TYPED", "CHR", "POS", "A0","A1", "exp_freq_A1",
                                      "SAMP_MAF"),
                   snp.df = data.frame(matrix(data = rep(NA, 16), ncol = 8 )),
                   snp.cols = c("snpID","TYPED","RSID", "POS","A0","A1","CHR"),
                   snp.ord = c("RSID","TYPED","CHR","POS","A0","A1"))
  
  class(cox_surv) <- c("Impute2CoxSurv", "PlinkGdsImpute2CoxSurv")
  
  return(cox_surv)
}

createSnpSpike.Impute2CoxSurv <- function(x, cox.params){
  rbind(c(rnorm(nrow(cox.params$pheno.file)-3), rep(NA, 3)),
        c(rnorm(nrow(cox.params$pheno.file)-4), rep(NA, 4)))
}


processSNPGenotypes.Impute2CoxSurv <- function(x, snp, genotypes, scanAnn, 
                                               exclude.snps, 
                                               cox.params, verbose) {
  # assign rsIDs (pasted with imputation status) as rows 
  # and sample ID as columns to genotype file
  dimnames(genotypes) <- list(paste(snp$TYPED, snp$RSID, sep=";"), 
                              scanAnn$ID_2)
  
  # Subset genotypes by given samples
  if(is.null(exclude.snps)){
    genotypes <- genotypes[,cox.params$ids]
  } else {
    genotypes <- genotypes[!snp$RSID %in% exclude.snps,cox.params$ids]
    snp <- snp[! snp$RSID %in% exclude.snps, ]
    if(verbose) message(length(exclude.snps), " SNPs are excluded based on the exclusion list provided by the user")
  }
  
  return(list(snp = snp, genotypes = genotypes))
}


getGenoData.Impute2CoxSurv <- function(x){
  
  if (x$keepGDS){
    gdsfile <- replaceFileExt(file.path = x$impute.file, ext = ".gds")
    snpfile <- replaceFileExt(file.path = x$impute.file, ext = ".snp.rdata")
    scanfile <- replaceFileExt(file.path = x$impute.file, ext = ".scan.rdata")
  } else {
    gdsfile <- tempfile(pattern="", fileext = ".gds")
    snpfile <- tempfile(pattern="", fileext = ".snp.rdata")
    scanfile <- tempfile(pattern="", fileext = ".scan.rdata")
    on.exit(unlink(c(gdsfile, snpfile, scanfile), recursive = TRUE), add=TRUE)
  }
  
  comp_time <- system.time(
    imputedDosageFile(input.files=c(x$impute.file, x$sample.file),
                      filename=gdsfile,
                      chromosome=as.numeric(x$chr),
                      input.type="IMPUTE2",
                      input.dosage=FALSE,
                      output.type = "dosage",
                      file.type="gds",
                      snp.annot.filename = snpfile,
                      scan.annot.filename = scanfile,
                      verbose=TRUE)
  )
  
  messageCompressionTime(comp_time)
  
  # read genotype
  ## need to add if statement about dimensions
  # set default "snp,scan" -- 
  # in GWASTools documentation say it needs to be in this orientation
  gds <- GdsGenotypeReader(gdsfile, genotypeDim="scan,snp")
  # close gds file on exit of the function
  # on.exit(close(gds), add=TRUE)
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