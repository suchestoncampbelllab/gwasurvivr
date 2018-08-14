survFit <- function(SNP, cox.params, print.covs){
        
    if(is.null(cox.params$pheno.file)){
        SNP <- SNP[!is.na(SNP)]
        X <- matrix(SNP, ncol = 1)
        Y <- cox.params$Y[!is.na(SNP)]
        ROWNAMES <- cox.params$ROWNAMES[!is.na(SNP)]
    }else{
        ## creating model matrix
        X <- cbind(SNP, cox.params$pheno.file)
        X <- X[!is.na(SNP),]
        Y <- cox.params$Y[!is.na(SNP)]
        ROWNAMES <- cox.params$ROWNAMES[!is.na(SNP)]
    }
        
        ## run fit with pre-defined parameters including INIT
        fit <- coxph.fit(X,
                         Y,
                         cox.params$STRATA,
                         cox.params$OFFSET,
                         cox.params$INIT, 
                         cox.params$CONTROL,
                         cox.params$WEIGHTS,
                         cox.params$METHOD, 
                         ROWNAMES)
        
        
        ## extract statistics
        if(print.covs=="only") {
            coef <- fit$coefficients[1]
            serr <- sqrt(diag(fit$var)[1])
            n.sample <- nrow(X)
            n.event <- sum(!grepl("[+]", as.character(Y)))
            res <- cbind(coef=coef, serr=serr, n.sample=n.sample, n.event=n.event)
            return(res)
            
        } else if(print.covs=="all"){
                
            coef <- fit$coefficients
            serr <- sqrt(diag(fit$var))
            res <- cbind(coef, serr)
            res.names <- dimnames(res)
            n.sample <- nrow(X)
            n.event <- sum(!grepl("[+]", as.character(Y)))
            res <- c(res, n.sample, n.event)
            names(res) <- c(paste(toupper(res.names[[2]][1]),
                                  res.names[[1]],
                                  sep="_"),
                            paste(toupper(res.names[[2]][2]),
                                  res.names[[1]],
                                  sep="_"),
                            "N", 
                            "N.EVENT"
                            )
            return(res)

        }
}
