#' Fit cox survival to all variants from a standard IMPUTE2 output after genotype imputation
#'
#' Performs survival analysis using Cox proportional hazard models on imputed genetic data from IMPUTE2 output
#'
#' @param impute.file character(1) of IMPUTE2 file 
#' @param sample.file character(1) of sample file affiliated with IMPUTE2 file
#' @param chromosome numeric(1) denoting chromosome number
#' @param infofile character(1) of info file affiliated with IMPUTE2 file
#' @param covfile data.frame(1) or matrix(1) comprising phenotype information; first column is required to be sample IDs
#' @param sample.ids character vector of sample IDs to keep in survival analysis
#' @param time character(1) of column name in covfile that represents the time interval of interest in the analysis
#' @param event character(1) of column name in covfile that represents the event of interest to be included in the analysis
#' @param covariates character vector with exact names of columns in covfile to include in analysis
#' @param outfile character(1) of output file name (do not include extension) 
#' @param flip.dosage logical(1) to flip which allele the dosage was calculated on, default=TRUE
#' @param verbose logical(1) for messages that describe which part of the analysis is currently being run
#' 
#' @return
#' Saves text file directly to disk that contains survival analysis results
#' 
#' @importFrom survival Surv coxph.fit
#' @import parallel
#' @import GWASTools
#' @import dplyr
#' 
#' @export

gdsCoxSurv <- function(impute.file,
                       sample.file,
                       chromosome,
                       infofile,
                       covfile, 
                       sample.ids, 
                       time, 
                       event,
                       covariates,
                       outfile,
                       flip.dosage,
                       verbose=TRUE){
        
        gdsfile <- paste0(outfile, ".gds")
        snpfile <- paste0(outfile, ".snp.rdata")
        scanfile <- paste0(outfile, ".scan.rdata")
        imputedDosageFile(input.files=c(impute.file, sample.file),
                          filename=gdsfile,
                          chromosome=as.numeric(chromosome),
                          input.type="IMPUTE2",
                          input.dosage=F,
                          file.type="gds",
                          snp.annot.filename = snpfile,
                          scan.annot.filename = scanfile)
        
        # read genotype
        gds <- GdsGenotypeReader(gdsfile, genotypeDim = "snp,scan")
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
                               header = TRUE,
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
        snp <- snp %>% left_join(infofile) 

        #         merge(infofile,
        #              snp,
        #              by=c("snp.index", "snpid", "rsid", "position"),
        #              all.y=TRUE)
        # 
        # # fix snp order
        # infofile$snpid_rsid <- paste(infofile$snpid, infofile$rsid, sep=";")
        # snp$snpid_rsid <- paste(snp$snpid, snp$rsid, sep=";")
        # 
        # snp <- snp[match(infofile$snpid_rsid, snp$snpid_rsid),]
        
        # remove extra cols
        #infofile$snpid_rsid <- NULL
        #snp$snpid_rsid <- NULL
        
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
        if(flip.dosage) genotypes <- 2 - genotypes
        
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
        
        colnames(scanAnn)[colnames(scanAnn)=="sex"] <- "sex.sample"
        scanAnn <- scanAnn[,c("ID_2", "missing", "sex.sample")]
        
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
        pval <- 2*pnorm(abs(z), lower.tail=FALSE)
        hr <- exp(snp.out[,1])
        lowerCI <- exp(snp.out[,1]-1.96*snp.out[,2])
        upperCI <- exp(snp.out[,1]+1.96*snp.out[,2])
        
        sres <- cbind(snp.out, hr, lowerCI, upperCI, z, pval, n.sample, n.event)
        colnames(sres) <- c("coef", "se.coef", "exp.coef", "lb", "ub", "z", "p.value", "n", "nevents")
        rownames(sres) <- NULL # remove rownames so we don't have a duplicated rownames issue
        
        res <- cbind(snp,sres)
        write.table(res, file=paste0(outfile, ".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
        
        if (verbose) message("Analysis completed on ", format(Sys.time(), "%Y-%m-%d"), " at ", format(Sys.time(), "%H:%M:%S"))
        return(res)
}






