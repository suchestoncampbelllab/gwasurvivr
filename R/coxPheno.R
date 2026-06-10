coxPheno <- function(pheno.file,
                     covariates,
                     id.column,
                     inter.term,
                     time.to.event,
                     event,
                     sample.ids,
                     verbose,
                     start.time = NULL){
    #### Phenotype data wrangling #####
    ## id column shold be provided!
    if (missing(id.column) ) {
        stop("Name of the ID column is not provided")
    }

    ## Coerce to a plain data.frame. tibbles / grouped_df (dplyr) keep different
    ## `[` semantics and a grouped tibble made `pheno.file[, covariates]` throw
    ## "subscript out of bounds" (#39). as.data.frame() also drops grouping.
    pheno.file <- as.data.frame(pheno.file, stringsAsFactors = FALSE)

    if (!id.column %in% colnames(pheno.file)) {
        stop("The id.column '", id.column, "' was not found in the covariate ",
             "file columns: ", paste(colnames(pheno.file), collapse = ", "))
    }

    ### SUBSET BY SAMPLE IDS IF GIVEN
    # user can provide null for sample.ids if not wishing to subset samples
    if(!is.null(sample.ids)){
        # only keep samples given with sample.ids argument
        pheno.file <- pheno.file[pheno.file[[id.column]] %in% sample.ids,]
        if(nrow(pheno.file) < 1) {
            stop("None of the provided sample.ids match values in the '",
                 id.column, "' column of the covariate file. Check that ",
                 "sample.ids and id.column refer to the same identifiers.")}
    }


    if(!is.null(covariates)){
        # covariates are defined in pheno.file
        if(!is.null(inter.term)) {
            if(!inter.term %in% colnames(pheno.file)) {
                stop("inter.term term is missing.")
            }
            covariates <- unique(covariates, inter.term)
        } 
        if(!is.null(inter.term)) covariates <- unique(covariates, inter.term)
        ok.covs <- colnames(pheno.file)[colnames(pheno.file) %in% covariates]
        if (verbose) message("Covariates included in the models are: ",
                             paste(ok.covs, collapse=", "))
        if(!is.null(inter.term) & verbose) {
            message("Models will include interaction term: SNP*", inter.term)
        }
        
        ### drop NAs (start.time is included when left-truncation is requested)
        pheno.file <- pheno.file[,c(id.column, start.time, time.to.event,
                                    event, ok.covs)]
        pheno.file <- pheno.file[complete.cases(pheno.file),]
        ids <- pheno.file[[id.column]]

        ## covariates should be numeric!
        pheno.file <- as.matrix(pheno.file[,c(start.time, time.to.event,
                                              event, ok.covs)])
        
        if (!is.numeric(pheno.file) ) {
            stop("Provided covariates must be numeric!\n",
                 "e.g. categorical variables should be recoded as indicator",
                 " or dummy variables.")
        }
    } else {
        
        ### drop NAs (start.time is included when left-truncation is requested)
        pheno.file <- pheno.file[,c(id.column, start.time, time.to.event, event)]
        pheno.file <- pheno.file[complete.cases(pheno.file),]
        ids <- pheno.file[[id.column]]

        ## time-event should be numeric!
        pheno.file <- as.matrix(pheno.file[,c(start.time, time.to.event, event)])
        
        
        if (!is.numeric(pheno.file) ) {
            stop("Time and event columns must be numeric!")
        }
    }
    

    # build coxph.fit parameters
    cox.params <- coxParam(pheno.file,
                           time.to.event,
                           event,
                           covariates,
                           ids,
                           verbose,
                           start.time = start.time)
    return(cox.params)
}
