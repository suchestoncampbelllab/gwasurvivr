impute2ProcessSNPGenotypes <- function(snp, genotypes, scanAnn, exclude.snps, 
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


plinkProcessSNPGenotypes <- function(snp, genotypes, scanAnn, exclude.snps = NULL, 
                                       cox.params, verbose) {
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


gdsProcessSNPGenotypes <- function(snp, genotypes, scanAnn, exclude.snps = NULL, 
                                     cox.params, verbose) {
  # assign rsIDs (pasted with imputation status) as rows
  # and sample ID as columns to genotype file
  dimnames(genotypes) <- list(paste(snp$snp, snp$rsID, sep=";"),
                              scanAnn$ID_2)
  
  # Subset genotypes by given samples
  genotypes <- genotypes[,cox.params$ids]
  
  return(list(snp = snp, genotypes = genotypes))
}
