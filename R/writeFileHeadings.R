writeFileHeadings <- function(cols, out.file, inter.term, snp.df, snp.spike,
                              print.covs, cox.params){
  
  # set up columns for output
  write.table( t(cols),
               paste0(out.file, ".snps_removed"),
               row.names = FALSE,
               col.names=FALSE,
               sep="\t",
               quote = FALSE,
               append = FALSE)

  colnames(snp.df) <- cols
  rownames(snp.df) <- NULL
  
  cox.out <- getSnpSpikeCoxOut(
    inter.term, 
    snp.spike = rbind(rnorm(nrow(cox.params$pheno.file)),
                      rnorm(nrow(cox.params$pheno.file))),
    cox.params, 
    print.covs)
  
  res.cols <- colnames(coxExtract(cox.out,
                                  snp.df,
                                  print.covs=print.covs))
  
  write.table(t(res.cols),
              paste0(out.file, ".coxph"),
              row.names = FALSE,
              col.names=FALSE,
              sep="\t",
              quote = FALSE,
              append = FALSE)
  
}