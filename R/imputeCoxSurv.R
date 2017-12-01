#'
#'
#'
#'
#'
#'
#'
#'
#'

imputeCoxSurv <- function(se, time, event, covariates){
        library(survival)
        library(SummarizedExperiment)
        colData(se) <- colData(se)[,c(time, event, covariates)]
        
        # maybe add operability with different assays
        formula <- paste0("Surv(time=",
                          time,
                          ", event=",
                          event,
                          ") ~ input.genotype + ",
                          paste(covariates, collapse=" + "))
        
        
        survFit <- function(input.genotype){
                env <- new.env()
                stat <- tryCatch(coxph(as.formula(formula), data=colData(se)),
                                 warning=function(warn) NA,
                                 error=function(err) NA
                )
                if(any(is.na(stat))){
                        rep(NA, 7)
                }else{
                        m <- summary(stat)
                        c(m$coef[1,], n=m$n, nevents=m$nevent)
                }
        }
        
        sres <- t( apply( assay(se), 1, function(i) survFit(i)  ))
        
        # change survival result matrix column names
        colnames(sres) <- c("coef", "exp.coef","se.coef","z", "p.value", "n", "nevents")
        
        # remove rownames so we don't have a duplicated rownames issue
        rownames(sres) <- NULL
        
        sres <- data.frame(sres)
        
        sres <- cbind(sres,
                      lb=exp(sres$coef-sres$se.coef),
                      ub=exp(sres$coef+sres$se.coef))
        
        
        sres <- sres[,c("coef", "se.coef", "exp.coef", "lb", "ub", "z", "p.value", "n", "nevents")]
        
        # add survival results into summarizedexperiment object
        mcols(rowRanges(se)) <- cbind(mcols(rowRanges(se)), data.frame(sres))
        
        metadata(se) <- list(formula=formula)
        ## we also may want to keep the formula info
        return(se)
}