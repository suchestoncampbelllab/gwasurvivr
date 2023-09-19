coxParam_intime <-
    function(pheno.file,
             time.start,
             time.stop,
             event,
             covariates,
             sample.ids,
             verbose) {
        #########################################
        ##### get the Ns ####
        # define Ns
        n.sample <- nrow(pheno.file)
        n.event <- sum(as.numeric(pheno.file[, event]))
        if (verbose)
            message(n.sample, " samples are included in the analysis")
        #########################################
        ##### build arguments for coxph.fit #####
        Y <-
            Surv(time = pheno.file[, time.start],
                 time2 = pheno.file[, time.stop],
                 event = pheno.file[, event])
        rownames(Y) <- as.character(seq_len(nrow(Y)))
        STRATA <- NULL
        CONTROL <- structure(
            list(
                eps = 1e-09,
                toler.chol = 1.81898940354586e-12,
                iter.max = 20L,
                # potentially select a more optimal max
                toler.inf = 3.16227766016838e-05,
                outer.max = 10L,
                timefix = TRUE
            ),
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
        #########################################
        
        #########################################
        ##### initialization ########
        #based on number of covariates e.g. 0,1 or more
        if (is.null(covariates)) {
            INIT <- NULL
            pheno.file <- matrix(pheno.file[, covariates], ncol = 1)
            pheno.file <- NULL
            
            
        } else if (length(covariates) == 1L) {
            pheno.file <- matrix(pheno.file[, covariates], ncol = 1)
            colnames(pheno.file) <- covariates
            INIT <- NULL
            init.fit <- coxph.fit(pheno.file,
                                  Y,
                                  STRATA,
                                  OFFSET,
                                  INIT,
                                  CONTROL,
                                  WEIGHTS,
                                  METHOD,
                                  ROWNAMES)
            
            INIT <- c(0,  init.fit$coefficients)
            
        } else {
            pheno.file <- pheno.file[, covariates]
            INIT <- NULL
            init.fit <- coxph.fit(pheno.file,
                                  Y,
                                  STRATA,
                                  OFFSET,
                                  INIT,
                                  CONTROL,
                                  WEIGHTS,
                                  METHOD,
                                  ROWNAMES)
            
            INIT <- c(0,  init.fit$coefficients)
        }
        
        #########################################
        cox.params <-
            list(
                pheno.file = pheno.file,
                Y = Y,
                STRATA = STRATA,
                OFFSET = OFFSET,
                CONTROL = CONTROL,
                WEIGHTS = WEIGHTS,
                METHOD = METHOD,
                ROWNAMES = ROWNAMES,
                INIT = INIT,
                n.sample = n.sample,
                n.event = n.event,
                ids = sample.ids
            )
        return(cox.params)
    }
