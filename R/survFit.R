#' @export 

### define survFit
.survFit <- function(input.genotype, params){
        
        ## creating model matrix
        X <- cbind(input.genotype, params$pheno.file[,covariates])
        
        ## run fit with pre-defined parameters including INIT
        fit <- coxph.fit(X,
                         params$Y,
                         params$STRATA,
                         params$OFFSET,
                         params$INIT, 
                         params$CONTROL,
                         params$WEIGHTS,
                         params$METHOD, 
                         params$ROWNAMES)
        
        ## extract statistics
        coef <- fit$coefficients[1]
        serr <- sqrt(diag(fit$var)[1])
        cbind(coef, serr)
}

