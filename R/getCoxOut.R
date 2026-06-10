getGenotypesCoxOut <- function(inter.term,
                               genotypes,
                               cl,
                               cox.params,
                               print.covs) {

        if (is.null(inter.term)) {
            if (is.matrix(genotypes)) {
                cox.out <- t(
                    parApply(
                        cl = cl,
                        X = genotypes,
                        MARGIN = 1,
                        FUN = survFit,
                        cox.params = cox.params,
                        print.covs = print.covs
                    )
                )
            } else {
                cox.out <- survFit(genotypes,
                                   cox.params = cox.params,
                                   print.covs = print.covs)
            }
        } else if (inter.term %in% colnames(cox.params$pheno.file)) {
            # The interaction covariate's name lives in the cox.params pheno
            # matrix columns. The previous code referenced an undefined
            # `covariates` symbol here, which errored on every interaction run
            # and produced empty output (#32).
            if (is.matrix(genotypes)) {
                cox.out <- t(
                    parApply(
                        cl = cl,
                        X = genotypes,
                        MARGIN = 1,
                        FUN = survFitInt,
                        cox.params = cox.params,
                        cov.interaction = inter.term,
                        print.covs = print.covs
                    )
                )
            } else {
                cox.out <- survFitInt(
                    genotypes,
                    cox.params = cox.params,
                    cov.interaction = inter.term,
                    print.covs = print.covs
                )
            }
        } else {
            stop("Interaction term '", inter.term, "' is not among the model ",
                 "covariates (", paste(colnames(cox.params$pheno.file),
                                       collapse = ", "),
                 "). It must be included in `covariates`.")
        }

    return(cox.out)

}


getSnpSpikeCoxOut <-
    function(inter.term,
             snp.spike,
             cox.params,
             print.covs) {
        if (is.null(inter.term)) {
            cox.out <- t(apply(
                snp.spike,
                1,
                survFit,
                cox.params = cox.params,
                print.covs = print.covs
            ))
        } else {
            cox.out <- t(
                apply(
                    snp.spike,
                    1,
                    survFitInt,
                    cox.params = cox.params,
                    cov.interaction = inter.term,
                    print.covs = print.covs
                )
            )
        }
        
    }
