survFitInt <- function(SNP,
                       cox.params,
                       cov.interaction,
                       print.covs) 
    {
    ## creating model matrix
    X <- cbind(INTER.TERM=cox.params$pheno.file[,cov.interaction]*SNP,
               SNP,
               cox.params$pheno.file)
    
    ## remove NA samples in genotype data
    X <- X[!is.na(SNP),]
    Y <- cox.params$Y[!is.na(SNP)]
    
    ## run fit with pre-defined parameters including INIT
    fit <- coxph.fit(X,
                     Y,
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
        n.sample <- length(Y)
        n.event <- sum(!grepl("[+]", as.character(Y)))
        return(list(res=res, n.sample=n.sample, n.event=n.event))
        
    } else if(print.covs=="some"){
        coef <- fit$coefficients
        serr <- sqrt(diag(fit$var))
        res <- cbind(coef, serr)
        drop <- colnames(cox.params$pheno.file)[!colnames(cox.params$pheno.file)
                                                %in% cov.interaction]
        res <- res[!rownames(res) %in% drop,]
        res.names <- dimnames(res)
        res <- c(res)
        names(res) <- c(paste(toupper(res.names[[2]][1]),
                              res.names[[1]],
                              sep="_"),
                        paste(toupper(res.names[[2]][2]),
                              res.names[[1]], 
                              sep="_"))
        n.sample <- length(Y)
        n.event <- sum(!grepl("[+]", as.character(Y)))
        return(list(res=res, n.sample=n.sample, n.event=n.event))
        
    } else if(print.covs=="all"){
        coef <- fit$coefficients
        serr <- sqrt(diag(fit$var))
        res <- cbind(coef, serr)
        res.names <- dimnames(res)
        res <- c(res)
        names(res) <- c(paste(toupper(res.names[[2]][1]),
                              res.names[[1]],
                              sep="_"),
                        paste(toupper(res.names[[2]][2]),
                              res.names[[1]],
                              sep="_"))
        n.sample <- length(Y)
        n.event <- sum(!grepl("[+]", as.character(Y)))
        return(list(res=res, n.sample=n.sample, n.event=n.event))
    }
}
