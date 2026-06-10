coxParam <-
    function(pheno.file,
             time.to.event,
             event,
             covariates,
             sample.ids,
             verbose,
             start.time = NULL) {
        #########################################
        # Resolve the survival fitter once; it is passed through cox.params so
        # it is available on every parallel worker. Counting-process /
        # left-truncated (start, stop] data requires agreg.fit; standard
        # right-censored data uses coxph.fit. Both share the same signature.
        coxph_fit <- if (is.null(start.time)) {
            .resolve_coxph_fit()
        } else {
            .resolve_agreg_fit()
        }
        #########################################
        ##### get the Ns ####
        # define Ns
        n.sample <- nrow(pheno.file)
        n.event <- sum(as.numeric(pheno.file[, event]))
        if (verbose)
            message(n.sample, " samples are included in the analysis")
        #########################################
        ##### build arguments for coxph.fit #####
        # If start.time is supplied, use left-truncated / counting-process
        # (start, stop] intervals; otherwise standard right-censored survival.
        if (is.null(start.time)) {
            Y <- Surv(time = pheno.file[, time.to.event],
                      event = pheno.file[, event])
        } else {
            Y <- Surv(time = pheno.file[, start.time],
                      time2 = pheno.file[, time.to.event],
                      event = pheno.file[, event])
        }
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
            # No covariates: no warm-start fit, and pheno.file carries no
            # covariate columns. n.sample (computed above) is the source of
            # truth for sample count -- see createSnpSpike.* (#6).
            INIT <- NULL
            pheno.file <- NULL


        } else if (length(covariates) == 1L) {
            pheno.file <- matrix(pheno.file[, covariates], ncol = 1)
            colnames(pheno.file) <- covariates
            INIT <- NULL
            init.fit <- coxph_fit(pheno.file,
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
            init.fit <- coxph_fit(pheno.file,
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
                ids = sample.ids,
                coxph_fit = coxph_fit
            )
        return(cox.params)
    }
