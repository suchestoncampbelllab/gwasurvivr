loadProcessWrite.PlinkGdsImpute2CoxSurv<- function(x, cl, cox.params) {
  
  writeFileHeadings(x, cox.params = cox.params)
  
  genoData <- getGenoData(x)
  
  results <- runOnChunks(x, genoData, cox.params, cl)
  
  return(list(snps_removed = results$snp.drop.n, 
              snps_analyzed = results$snp.n))
}

writeFileHeadings <- function(x, cox.params){
  
  # set up columns for output
  write.table( t(x$columnHeadings),
               paste0(x$out.file, ".snps_removed"),
               row.names = FALSE,
               col.names=FALSE,
               sep="\t",
               quote = FALSE,
               append = FALSE)
  
  colnames(x$snp.df) <- x$columnHeadings
  rownames(x$snp.df) <- NULL
  
  cox.out <- getSnpSpikeCoxOut(
    x$inter.term, 
    snp.spike = createSnpSpike(x, cox.params),
    cox.params, 
    x$print.covs)
  
  res.cols <- colnames(coxExtract(cox.out,
                                  x$snp.df,
                                  print.covs=x$print.covs))
  
  write.table(t(res.cols),
              paste0(x$out.file, ".coxph"),
              row.names = FALSE,
              col.names=FALSE,
              sep="\t",
              quote = FALSE,
              append = FALSE)
  
}