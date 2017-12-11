#' Fit cox survival to all variants in a SummarizedExperiment object
#'
#' Performs survival analysis using Cox proportional hazard models on imputed genetic data stored in SummarizedExperiment object
#'
#' @param se SummarizedExperiment object as generated with `readImputeGds` and `addCov` functions 
#' or includes 1 assay containing allele dosages; rowData with SNP information;
#' colData including time, event and covariates.
#' @param time character(1) string that matches time column name in pheno.file
#' @param event character(1) string that matches event column name in pheno.file
#' @param covariates character vector with matching column names in pheno.file of covariates of interest
#' 
#' @return
#' Generates se where survival results are kept in rowData.
#'
#'
#'

imputeCoxSurv <- function(se, time, event, covariates){
    
    message("Analysis started on ", format(Sys.time(), "%d/%b/%y"), " at ", format(Sys.time(), "%H:%M:%S"))
    
    pheno.file <- as.matrix(colData(se)[,c(time, event, covariates)])
    
    # define Ns
    n.sample <- nrow(pheno.file)
    n.event <- sum(pheno.file[,event])
    message("Analysis running for ", n.sample, " samples.")
        
    # covariates are defined in pheno.file
    ok.covs <- colnames(pheno.file)[colnames(pheno.file) %in% covariates]
    message("Covariates included in the models are: ", paste(ok.covs, collapse=", "))
    message("If your covariates of interest are not included in the model\nplease stop the analysis and make sure user defined covariates\nmatch the column names in the colData(se)")
        
    ### Check fo snps with MAF = 0  ###  
    indx <- sort(unique(c(which(matrixStats::rowSds(assay(se,1)) == 0))))
    
    ### which() may produce integer(0), have a check
    if(length(indx) != 0){
        
        message(length(indx), " SNPs were removed from the analysis for not meeting the threshold criteria.")
        
        ## Remove MAF=0 snps
        se <- se[-indx,]
        
        # save list of snps that were removed
        metadata(se) <- list(removedSNPs=rowRanges(se)$rsid[indx])
    }
    
    message("Running survival models...")
        
    ### build arguments for coxph.fit ###
    Y <- Surv(time=colData(se)[[time]], event=colData(se)[[event]])
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
    ROWNAMES <- rownames(pheno.file)
        
    # define INIT
    INIT <- NULL
    init.fit <- coxph.fit(pheno.file[,covariates], 
                          Y, STRATA, OFFSET, INIT, 
                          CONTROL, WEIGHTS, METHOD, ROWNAMES)
    
    INIT <- c(0,  init.fit$coefficients)
        
    ### define survFit
    survFit <- function(input.genotype){
        
        ## only thing that's changing
        X <- cbind(input.genotype,pheno.file[,covariates])
        
        ## run fit with pre-defined parameters including INIT
        fit <- coxph.fit(X, Y, STRATA, OFFSET, INIT, 
                         CONTROL, WEIGHTS, METHOD, ROWNAMES)
            
        ## extract statistics
        coef <- fit$coefficients[1]
        serr <- sqrt(diag(fit$var)[1])
        cbind(coef, serr)
    }
    
    ## detect cores for a single machine
    cl <- makeForkCluster(detectCores())
    
    ## apply survival function
    snp.out <- t(parApply(cl=cl, X=assay(se,1), MARGIN=1, FUN=survFit))
    # snp.out <- parRapply(cl=cl, x=assay(se,1), FUN=survFit)
    
    z <- snp.out[,1]/snp.out[,2]
    pval <- 2*pnorm(abs(z), lower.tail=F)
    hr <- exp(snp.out[,1])
    lowerCI <-exp(snp.out[,1]-1.96*snp.out[,2])
    upperCI <-exp(snp.out[,1]+1.96*snp.out[,2])
    
    sres <- cbind(snp.out, hr, lowerCI, upperCI, z, pval, n.sample, n.event)
    colnames(sres) <- c("coef", "se.coef", "exp.coef", "lb", "ub", "z", "p.value", "n", "nevents")
    rownames(sres) <- NULL # remove rownames so we don't have a duplicated rownames issue
    
    ## add survival results into summarizedexperiment object
    mcols(rowRanges(se)) <- cbind(mcols(rowRanges(se)), data.frame(sres)[,,drop=T])
        
    return(se)
    message("Analysis completed on ", format(Sys.time(), "%d/%b/%y"), " at ", format(Sys.time(), "%H:%M:%S"))
}
