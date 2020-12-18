createPlinkCoxSurv <- function(b.file,
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
  
  cox_surv <- list(b.file = b.file,
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
                   flip.dosage =flip.dosage,
                   verbose = verbose,
                   clusterObj = clusterObj,
                   keepGDS = keepGDS,
                   columnHeadings = c("RSID","CHR", "POS","A0","A1", "exp_freq_A1","SAMP_MAF"),
                   snp.df = data.frame(t(rep(NA, 7))),
                   snp.cols = c("snpID","RSID","CHR","POS","A0","A1"),
                   snp.ord = c("RSID","CHR","POS","A0","A1"))
  
  class(cox_surv) <- c("PlinkCoxSurv", "PlinkGdsImpute2CoxSurv")
  
  return(cox_surv)
}

createSnpSpike.PlinkCoxSurv <- function(x, cox.params){
  rbind(rnorm(nrow(cox.params$pheno.file)),
        rnorm(nrow(cox.params$pheno.file)))
}


processSNPGenotypes.PlinkCoxSurv <- function(x, snp, genotypes, scanAnn, 
                                             exclude.snps = NULL, 
                                             cox.params, verbose){
  
    # assign rsIDs (pasted with imputation status) as rows
    # and sample ID as columns to genotype file
    dimnames(genotypes) <- list(snp$RSID,
                                scanAnn$scanID)
    
    # Subset genotypes by given samples
    blankSNPs <- snp$A0 == "0" & snp$A1 == "0"
    genotypes <- genotypes[!blankSNPs,cox.params$ids]
    snp <- snp[!blankSNPs,]
    
    return(list(snp = snp, genotypes = genotypes))
}


getGenoData.PlinkCoxSurv <- function(x) {
  
  bim.file <- replaceFileExt(file.path =  x$b.file, ext = ".bim")
  fam.file <- replaceFileExt(file.path =  x$b.file, ext = ".fam")
  
  if (x$keepGDS){
    gdsfile <- replaceFileExt(file.path =  x$b.file, ext = ".gds")
    on.exit(unlink(gdsfile, recursive = TRUE), add=TRUE)
  } else {
    gdsfile <- tempfile(pattern="", fileext = ".gds")
  }
  
  
  comp_time <- system.time(snpgdsBED2GDS( x$b.file, 
                                         fam.file,
                                         bim.file,
                                         gdsfile,
                                         cvt.chr="int",
                                         cvt.snpid="int",
                                         verbose=TRUE)
  )
  
  messageCompressionTime(comp_time)
  
  # read genotype
  ## need to add if statement about dimensions
  # set default "snp,scan" -- 
  # in GWASTools documentation say it needs to be in this orientation
  gds <- GdsGenotypeReader(gdsfile,
                           YchromCode=24L, 
                           XYchromCode=25L)
  
  
  scanID <- getScanID(gds)
  scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID,
                                                  stringsAsFactors=FALSE))
  snpID <- getSnpID(gds)
  chromosome <- getChromosome(gds)
  position <- getPosition(gds)
  alleleA <- getAlleleA(gds)
  alleleB <- getAlleleB(gds)
  rsID <- getVariable(gds, "snp.rs.id")
  # requires snpID
  snpAnnot <- SnpAnnotationDataFrame(data.frame(snpID,
                                                rsID,
                                                chromosome,
                                                position,
                                                alleleA,
                                                alleleB,
                                                stringsAsFactors=FALSE),
                                     YchromCode=24L,
                                     XYchromCode=25L)
  genoData <- GenotypeData(gds,
                           scanAnnot=scanAnnot, 
                           snpAnnot=snpAnnot)
  
  return(genoData)
}