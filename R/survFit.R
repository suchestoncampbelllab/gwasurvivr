survFit <- function(SNP, cox.params, print.covs){
        
    if(is.null(cox.params$pheno.file)){
        X <- matrix(SNP, ncol = 1)
    }else{
        ## creating model matrix
        X <- cbind(SNP, cox.params$pheno.file)
    }
        
        ## run fit with pre-defined parameters including INIT
        fit <- coxph.fit(X,
                         cox.params$Y,
                         cox.params$STRATA,
                         cox.params$OFFSET,
                         cox.params$INIT, 
                         cox.params$CONTROL,
                         cox.params$WEIGHTS,
                         cox.params$METHOD, 
                         cox.params$ROWNAMES)
        
        
        ## extract statistics
        if(print.covs=="only") {
            coef <- fit$coefficients[1]
            serr <- sqrt(diag(fit$var)[1])
            res <- cbind(coef=coef, serr=serr)
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
