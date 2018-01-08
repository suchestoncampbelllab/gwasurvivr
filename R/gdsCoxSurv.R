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


gdsCoxSurv <- function(gdsfile, 
                       scanfile, 
                       snpfile, 
                       infofile, 
                       covfile, 
                       sample.ids, 
                       time, 
                       event,
                       covariates,
                       effect.allele,
                       verbose=TRUE){
        # read genotype
        gds <- GdsGenotypeReader(gdsfile)
        # close gds file on exit of the function
        on.exit(close(gds))
        # read in snp data
        snpAnnot <- getobj(snpfile)
        # read scan
        scanAnnot <- getobj(scanfile)
        # put into GenotypeData coding 
        genoData <- GenotypeData(gds,
                                 snpAnnot=snpAnnot,
                                 scanAnnot=scanAnnot)
        # store genotype, sample info and, and snp info
        genotypes <- getGenotype(genoData)
        # grab map data and rename some columns
        snp <- getAnnotation(getSnpAnnotation(genoData)) 
        colnames(snp)[colnames(snp)=="snpID"] <- "snp.index"
        colnames(snp)[colnames(snp)=="rsID"] <- "rsid"
        colnames(snp)[colnames(snp)=="snp"] <- "snpid"
        # grab sample file data
        scanAnn <- getAnnotation(getScanAnnotation(genoData))
        # read in info table
        infofile <- read.table(infofile,
                               header = FALSE,
                               stringsAsFactors = FALSE)
        colnames(infofile) <- c("snpid",
                                "rsid",
                                "position",
                                "exp_freq_a1",
                                "info",
                                "certainty",
                                "type",
                                "info_type0", 
                                "concord_type0",
                                "r2_type0")
        # select columns of interest
        infofile <- infofile[,c("snpid",
                                "rsid",
                                "position",
                                "exp_freq_a1",
                                "info",
                                "certainty")]
        
        infofile <- infofile[infofile$rsid %in% snp$rsid,]
        
        # add snp.index so we can avoid some duplicate warnings
        #        infofile[duplicated(infofile$rsid) | duplicated(infofile$rsid, fromLast=TRUE),]
        
        
        infofile$snp.index <- seq_len(nrow(infofile))
        
        
        
        # merge snp file with info file
        snp <- merge(infofile,
                     snp,
                     by=c("snp.index", "snpid", "rsid", "position"),
                     all.y=TRUE)
        
        # change order into what we want final outside to be
        colnames(snp)[colnames(snp)=="chromosome"] <- "chr"
        colnames(snp)[colnames(snp)=="alleleA"] <- "allele0"
        colnames(snp)[colnames(snp)=="alleleB"] <- "allele1"
        
        snp <- snp[,c("snp.index",
                      "chr",
                      "position",
                      "snpid",
                      "rsid",
                      "allele0",
                      "allele1", 
                      "exp_freq_a1", 
                      "info",
                      "certainty")]
        
        dimnames(genotypes) <- list(snp$rsid, scanAnn$ID_2)
        
        # flip dosage
        genotypes <- 2 - genotypes
        
        # add covariates to scan file
        # covfile <- read.table(covfile,
        #                       header=TRUE,
        #                       stringsAsFactors=FALSE,
        #                       sep="\t")
        
        colnames(covfile)[1] <- "ID_2" 
        
        # only keep samples with complete data
        covfile <- covfile[complete.cases(covfile),]
        covfile <- covfile[covfile$ID_2 %in% sample.ids,]
        # subset genotype data for patients of interest
        genotypes <- genotypes[,colnames(genotypes) %in% covfile[[1]]]
        
        scanAnn <- scanAnn[,c("ID_2", "missing", "sex")]
        colnames(scanAnn)[colnames(scanAnn)=="sex"] <- "sexFromSample"
        
        
        scanAnn <- merge(scanAnn,
                         covfile,
                         by="ID_2",
                         all.x=TRUE)
        
        # fix order
        scanAnn <- scanAnn[match(colnames(genotypes), scanAnn$ID_2),]
        
        ## STARTING ANALYSIS PORTION
        if (verbose) message("Analysis started on ", format(Sys.time(), "%Y-%m-%d"), " at ", format(Sys.time(), "%H:%M:%S"))
        
        pheno.file <- as.matrix(scanAnn[,c(time, event, covariates)])
        
        # define Ns
        n.sample <- nrow(pheno.file)
        n.event <- sum(pheno.file[,event])
        if (verbose) message("Analysis running for ", n.sample, " samples.")
        
        # covariates are defined in pheno.file
        ok.covs <- colnames(pheno.file)[colnames(pheno.file) %in% covariates]
        if (verbose) message("Covariates included in the models are: ", paste(ok.covs, collapse=", "))
        if (verbose) message("If your covariates of interest are not included in the model\nplease stop the analysis and make sure user defined covariates\nmatch the column names in the colData(se)")
        
        ### Check of snps with MAF = 0  ###  
        indx <- sort(unique(c(which(matrixStats::rowSds(genotypes) == 0))))
        
        ### which() may produce integer(0), have a check
        if(length(indx) != 0){
                message(length(indx), " SNPs were removed from the analysis for not meeting the threshold criteria.")
                ## Remove MAF=0 snps
                genotypes <- genotypes[-indx,]
                # remove from snp list too that we will merge back later
                snp <- snp[-indx,]
                # save list of snps that were removed
                rm.snps <- list(removedSNPs=snp$rsid[indx])
        }
        
        if (verbose) message("Running survival models...")
        
        ### build arguments for coxph.fit ###
        Y <- Surv(time=scanAnn[[time]], event=scanAnn[[event]])
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
                              Y, 
                              STRATA,
                              OFFSET,
                              INIT, 
                              CONTROL,
                              WEIGHTS,
                              METHOD, 
                              ROWNAMES)
        
        INIT <- c(0,  init.fit$coefficients)
        
        ### define survFit
        survFit <- function(input.genotype){
                
                ## creating model matrix
                X <- cbind(input.genotype, pheno.file[,covariates])
                
                ## run fit with pre-defined parameters including INIT
                fit <- coxph.fit(X,
                                 Y,
                                 STRATA,
                                 OFFSET,
                                 INIT, 
                                 CONTROL,
                                 WEIGHTS,
                                 METHOD, 
                                 ROWNAMES)
                
                ## extract statistics
                coef <- fit$coefficients[1]
                serr <- sqrt(diag(fit$var)[1])
                cbind(coef, serr)
        }
        
        ## detect cores for a single machine
        cl <- makeForkCluster(detectCores())
        on.exit(stopCluster(cl))
        ## apply survival function
        snp.out <- t(parApply(cl=cl, X=genotypes, MARGIN=1, FUN=survFit))
        # snp.out <- parRapply(cl=cl, x=assay(se,1), FUN=survFit)
        
        z <- snp.out[,1]/snp.out[,2]
        pval <- 2*pnorm(abs(z), lower.tail=F)
        hr <- exp(snp.out[,1])
        lowerCI <- exp(snp.out[,1]-1.96*snp.out[,2])
        upperCI <- exp(snp.out[,1]+1.96*snp.out[,2])
        
        sres <- cbind(snp.out, hr, lowerCI, upperCI, z, pval, n.sample, n.event)
        colnames(sres) <- c("coef", "se.coef", "exp.coef", "lb", "ub", "z", "p.value", "n", "nevents")
        rownames(sres) <- NULL # remove rownames so we don't have a duplicated rownames issue
        
        res <- cbind(snp,sres)
        
        return(res)
        message("Analysis completed on ", format(Sys.time(), "%Y/%m/%d"), " at ", format(Sys.time(), "%H:%M:%S"))
}


