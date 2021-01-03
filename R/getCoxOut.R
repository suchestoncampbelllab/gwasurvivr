getGenotypesCoxOut <- function(inter.term, genotypes, cl, cox.params,
                        print.covs) {
  
  if(is.null(inter.term)){
    if(is.matrix(genotypes)){
      cox.out <- t(parApply(cl=cl,
                            X=genotypes, 
                            MARGIN=1, 
                            FUN=survFit, 
                            cox.params=cox.params,
                            print.covs=print.covs))
    } else {
      cox.out <- survFit(genotypes,
                         cox.params=cox.params,
                         print.covs=print.covs) 
    }
  } else if(inter.term %in% covariates) {
    if(is.matrix(genotypes)){
      cox.out <- t(parApply(cl=cl,
                            X=genotypes,
                            MARGIN=1,
                            FUN=survFitInt, 
                            cox.params=cox.params, 
                            cov.interaction=inter.term,
                            print.covs=print.covs))
    } else {
      cox.out <- survFitInt(genotypes,
                            cox.params=cox.params,
                            cov.interaction=inter.term, 
                            print.covs=print.covs)
    }
  }
  
  return(cox.out)
        
}


getSnpSpikeCoxOut <- function(inter.term, snp.spike, cox.params, print.covs){
  
  if(is.null(inter.term)){
    cox.out <- t(apply(snp.spike, 1, survFit,
                       cox.params=cox.params,
                       print.covs=print.covs))
  } else {
    cox.out <- t(apply(snp.spike,
                       1,
                       survFitInt,
                       cox.params=cox.params,
                       cov.interaction=inter.term, 
                       print.covs=print.covs) )
  }
  
}

