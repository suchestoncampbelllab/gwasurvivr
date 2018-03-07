#' @export 

# additional helper functions for survival models
.coxParam <- function(pheno.file, time.to.event, event, covariates, sample.ids){
        ### build arguments for coxph.fit ###
        Y <- Surv(time=pheno.file[,time.to.event], event=pheno.file[,event])
        rownames(Y) <- as.character(seq_len(nrow(Y)))
        STRATA <- NULL
        CONTROL <- structure(
                list(eps = 1e-09,
                     toler.chol = 1.81898940354586e-12,
                     iter.max = 20L, # potentially select a more optimal max
                     toler.inf = 3.16227766016838e-05,
                     outer.max = 10L, 
                     timefix = TRUE),
                .Names = c(
                        "eps",
                        "toler.chol",
                        "iter.max", 
                        "toler.inf", 
                        "outer.max", 
                        "timefix"
                )
        )
        
        OFFSET <- NULL
        WEIGHTS <- NULL
        METHOD <- "efron"
        ROWNAMES <- sample.ids
        
        # define INIT
        INIT <- NULL
        init.fit <- coxph.fit(pheno.file[,covariates], 
                              Y, 
                              STRATA,
                              OFFSET,
                              INIT, 
                              CONTROL,
                              WEIGHTS,
                              METHOD, 
                              ROWNAMES)
        
        INIT <- c(0,  init.fit$coefficients)
        
        params <- list(Y=Y, STRATA=STRATA, OFFSET=OFFSET, CONTROL=CONTROL, WEIGHTS=WEIGHTS, METHOD=METHOD, ROWNAMES=ROWNAMES, INIT=INIT)
        return(params)
}
