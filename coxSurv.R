## Survival Function
library(survival)


coxSurv <- function(se, ncore=1){
        library(survival)
        .env <- environment()

        # maybe add operability with different assays
        formula <- paste0("Surv(time=",
                      names(colData(se))[2],
                      ", event=",
                      names(colData(se))[3],
                      ") ~ input.genotype + ",
                      paste(names(colData(se)[4:length(names(colData(se)))]), collapse=" + "))
        
        survFit <- function(input.genotype){
                stat <- tryCatch(coxph(as.formula(formula), data=data.frame(colData(se))),
                                 warning=function(warn)NA,
                                 error=function(err)NA
                                 )
                if(any(is.na(stat))){
                        rep(NA, 6)
                }else{
                        m <- summary(stat)
                        c(m$coef[1,], n=m$n)
                }
        }
        
        #parallel when ncore>1
        if(ncore==1){
                sres <- t( apply( assay(se), 1, function(i) survFit(i)  ))
        }else if(ncore>1){
                cl <- makeCluster(ncore, type="SOCK")
                clusterExport(cl, c(colData(se),"formula", "survFit"), envir=.env)
                clusterEvalQ(cl, library(survival))
                sres <- t(parApply(cl, assay(se), 1, function(i) survFit(i)) )
                stopCluster(cl)
        }else{
                stop("ncore should be >= 1")
        }
        
        # change survival result matrix column names
        colnames(sres) <- c("coef", "exp.coef","se.coef","z", "p-value", "n")
        # remove rownames so we don't have a duplicated rownames issue
        rownames(sres) <- NULL
        
        # add survival results into summarizedexperiment object
        mcols(rowRanges(se)) <- cbind(mcols(rowRanges(se)), data.frame(sres))
        
        metadata(se) <- list(formula=formula)
        ## we also may want to keep the formula info
        
        return(se)
        
        
}

library(microbenchmark)
microbenchmark(se.surv <- coxSurv(se), times=1)
