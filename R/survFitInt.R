survFitInt <-
    function(SNP,
             cox.params,
             cov.interaction,
             print.covs = "only"
             ) {
        
    ## creating model matrix
    X <- cbind(INTER.TERM=cox.params$pheno.file[,cov.interaction]*SNP, SNP, cox.params$pheno.file)
    
    ## run fit with pre-defined parameters including INIT
    fit <- coxph.fit(X,
                     cox.params$Y,
                     cox.params$STRATA,
                     cox.params$OFFSET,
                     c(0,cox.params$INIT), 
                     cox.params$CONTROL,
                     cox.params$WEIGHTS,
                     cox.params$METHOD, 
                     cox.params$ROWNAMES)
    
    ## extract statistics
    if(print.covs=="only") {
        coef <- fit$coefficients[1]
        serr <- sqrt(diag(fit$var)[1])
        res <- cbind(coef, serr)
        return(res)
        
    } else if(print.covs=="some"){
        coef <- fit$coefficients
        serr <- sqrt(diag(fit$var))
        res <- cbind(coef, serr)
        drop <- colnames(cox.params$pheno.file)[!colnames(cox.params$pheno.file) %in% cov.interaction]
        res <- res[!rownames(res) %in% drop,]
        res.names <- dimnames(res)
        res <- c(res)
        names(res) <- c(paste(toupper(res.names[[2]][1]), res.names[[1]], sep="_"),
                        paste(toupper(res.names[[2]][2]), res.names[[1]], sep="_"))
        return(res)
        
    } else if(print.covs=="all"){
        coef <- fit$coefficients
        serr <- sqrt(diag(fit$var))
        res <- cbind(coef, serr)
        res.names <- dimnames(res)
        res <- c(res)
        names(res) <- c(paste(toupper(res.names[[2]][1]), res.names[[1]], sep="_"),
                        paste(toupper(res.names[[2]][2]), res.names[[1]], sep="_"))
        return(res)

    }
}
