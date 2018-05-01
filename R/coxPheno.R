coxPheno <- function(pheno.file, covariates, id.column, inter.term, time.to.event, event, sample.ids, verbose){
    #### Phenotype data wrangling #####
    ## id column shold be provided!
    if (missing(id.column) ) {
        stop("Name of the ID column is not provided")
    }
    
    ### SUBSET BY SAMPLE IDS IF GIVEN
    # user can provide null for sample.ids if not wishing to subset samples
    if(!is.null(sample.ids)){
        # only keep samples given with sample.ids argument
        pheno.file <- pheno.file[pheno.file[[id.column]] %in% sample.ids,]
    }
    

    if(!is.null(covariates)){
        # covariates are defined in pheno.file
        if(!is.null(inter.term)) {
            if(!inter.term %in% colnames(pheno.file)) stop("inter.term term is missing.")
            covariates <- base::unique(covariates, inter.term)
        } 
        if(!is.null(inter.term)) covariates <- base::unique(covariates, inter.term)
        ok.covs <- colnames(pheno.file)[colnames(pheno.file) %in% covariates]
        if (verbose) message("Covariates included in the models are: ", paste(ok.covs, collapse=", "))
        if(!is.null(inter.term) & verbose) message("Models will include interaction term: SNP*",inter.term)
        
        ### drop NAs
        pheno.file <- pheno.file[,c(id.column, time.to.event, event, ok.covs)]
        pheno.file <- pheno.file[complete.cases(pheno.file),]
        ids <- pheno.file[[id.column]]
        
        ## covariates should be numeric!
        pheno.file <- as.matrix(pheno.file[,c(time.to.event, event)])
        
        if (!is.numeric(pheno.file) ) {
            stop("Provided covariates must be numeric!\ne.g. categorical variables should be recoded as indicator or dummy variables.")
        }
    } else {
        
        ### drop NAs
        pheno.file <- pheno.file[,c(id.column, time.to.event, event)]
        pheno.file <- pheno.file[complete.cases(pheno.file),]
        ids <- pheno.file[[id.column]]
        
        ## time-event should be numeric!
        pheno.file <- as.matrix(pheno.file[,c(time.to.event, event)])
        
        
        if (!is.numeric(pheno.file) ) {
            stop("Time and event columns must be numeric!")
        }
    }
    

    # build coxph.fit parameters
    cox.params <- coxParam(pheno.file, time.to.event, event, covariates, ids, verbose)
    return(cox.params)
}
