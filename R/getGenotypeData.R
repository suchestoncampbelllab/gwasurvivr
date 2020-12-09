getGdsGenotypeData <- function(gdsfile){
  
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

getPlinkGenoData <- function(gdsfile, keepGDS, b.file) {
  
  bim.file <- replaceFileExt(file.path = b.file, ext = ".bim")
  fam.file <- replaceFileExt(file.path = b.file, ext = ".fam")
  
  if (keepGDS){
    gdsfile <- replaceFileExt(file.path = b.file, ext = ".gds")
    on.exit(unlink(gdsfile, recursive = TRUE), add=TRUE)
  } else {
    gdsfile <- tempfile(pattern="", fileext = ".gds")
  }
  
  
  comp_time <- system.time(snpgdsBED2GDS(b.file, 
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


getImpute2GenoData <- function(gdsfile, snpfile, scanfile, impute.file,
                               keepGDS, sample.file, chr){
  
  if (keepGDS){
    gdsfile <- replaceFileExt(file.path = impute.file, ext = ".gds")
    snpfile <- replaceFileExt(file.path = impute.file, ext = ".snp.rdata")
    scanfile <- replaceFileExt(file.path = impute.file, ext = ".scan.rdata")
  } else {
    gdsfile <- tempfile(pattern="", fileext = ".gds")
    snpfile <- tempfile(pattern="", fileext = ".snp.rdata")
    scanfile <- tempfile(pattern="", fileext = ".scan.rdata")
    on.exit(unlink(c(gdsfile, snpfile, scanfile), recursive = TRUE), add=TRUE)
  }
  
  comp_time <- system.time(
    imputedDosageFile(input.files=c(impute.file, sample.file),
                      filename=gdsfile,
                      chromosome=as.numeric(chr),
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