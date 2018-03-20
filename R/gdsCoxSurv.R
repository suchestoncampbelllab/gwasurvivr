#' Fit cox survival to all variants from a standard IMPUTE2 output after genotype imputation
#'
#' Performs survival analysis using Cox proportional hazard models on imputed genetic data from IMPUTE2 output
#'
#' @param impute.file character(1) of IMPUTE2 file 
#' @param sample.file character(1) of sample file affiliated with IMPUTE2 file
#' @param chromosome numeric(1) denoting chromosome number
#' @param covfile data.frame(1) or matrix(1) comprising phenotype information; first column is required to be sample IDs
#' @param sample.ids character(1) vector of sample IDs to keep in survival analysis
#' @param time.to.event character(1) of column name in covfile that represents the time interval of interest in the analysis
#' @param event character(1) of column name in covfile that represents the event of interest to be included in the analysis
#' @param covariates character(1) vector with exact names of columns in covfile to include in analysis
#' @param outfile character(1) of output file name (do not include extension) 
#' @param maf.filter numeric(1) to filter minor allele frequency (MAF)
#' @param flip.dosage logical(1) to flip which allele the dosage was calculated on, default=TRUE
#' @param verbose logical(1) for messages that describe which part of the analysis is currently being run
#'
#' @return
#' Saves text file directly to disk that contains survival analysis results
#'  
#' @importFrom survival Surv coxph.fit
#' @import parallel
#' @import GWASTools
#' @import matrixStats
#' 
#' @export

gdsCoxSurv <- function(impute.file,
                       sample.file,
                       chromosome,
                       covfile, 
                       sample.ids, 
                       time.to.event, 
                       event,
                       covariates,
                       outfile,
                       maf.filter=NULL,
                       flip.dosage=TRUE,
                       verbose=TRUE
                       ){
        
        
        gdsfile <- tempfile(pattern="", fileext = ".gds")
        snpfile <- tempfile(pattern="", fileext = ".snp.rdata")
        scanfile <- tempfile(pattern="", fileext = ".scan.rdata")
        
        GWASTools::imputedDosageFile(input.files=c(impute.file, sample.file),
                                     filename=gdsfile,
                                     chromosome=as.numeric(chromosome),
                                     input.type="IMPUTE2",
                                     input.dosage=FALSE,
                                     file.type="gds",
                                     snp.annot.filename = snpfile,
                                     scan.annot.filename = scanfile)
        
        # read genotype
        ## need to add if statement about dimensions
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

        # grab sample file data
        scanAnn <- getAnnotation(getScanAnnotation(genoData))
        
        # assign rsIDs (pasted with imputation status) as rows 
        # and sample ID as columns to genotype file
        dimnames(genotypes) <- list(paste(snp$snp, snp$rsID, sep=";"), 
                                    scanAnn$ID_2)
        # flip dosage
        if(flip.dosage) genotypes <- 2 - genotypes
        
        # user needs to input covfile 
        colnames(covfile)[1] <- "ID_2" 
        
        # only keep samples given with sample.ids argument
        covfile <- covfile[covfile$ID_2 %in% sample.ids,]
        
        # subset genotype data for patients of interest
        genotypes <- genotypes[,covfile[[1]]]
        
        # calculate MAF
        snp$exp_freq_A1 <- round(1-matrixStats::rowMeans2(genotypes)*0.5,3)
        
        # calculate info score
        obs.mean <- matrixStats::rowMeans2(genotypes)
        obs.var <- matrixStats::rowVars(genotypes)
        p <- obs.mean/2
        p_all <- 2*p*(1-p)
        info.score <- round(obs.var/p_all,3)
        info.score[info.score>1] <- 1
        snp$info <- info.score
        
        # rearrange columns
        colnames(snp) <- c("snpID", "snpid", "rsid", "position", 
                           "A0", "A1", "chr", "exp_freq_A1", "info")
        
        snp <- snp[,c("chr",
                      "position",
                      "snpid",
                      "rsid",
                      "A0",
                      "A1", 
                      "exp_freq_A1", 
                      "info")]
        
        ### Check snps for MAF = 0  ###
        # if MAF threshold is set, subset by the given value,
        # otherwise check for SD = 0 and remove
        
        if(!is.null(maf.filter)){
                maf.idx <- snp$exp_freq_A1<maf.filter & snp$exp_freq_A1>(1-maf.filter)
                rm.snps <- snp[maf.idx, c("snpid", "rsid", "exp_freq_A1", "info")]
                snp <- snp[!maf.idx,]
                genotypes <- genotypes[!maf.idx,]
                
                if(verbose) message(sum(maf.idx), " SNPs were removed from the analysis for not meeting the given MAF < ", maf.filter)
                
        } else {
            indx <- sort(unique(c(which(matrixStats::rowSds(genotypes) == 0))))

            ### which() may produce integer(0), have a check
            if(length(indx) != 0){
                    ## Remove MAF=0 snps
                    genotypes <- genotypes[-indx,]
                    # remove from snp list too that we will merge back later
                    snp <- snp[-indx,]
                    # save list of snps that were removed
                    rm.snps <- snp[-indx, c("snpid", "rsid", "exp_freq_A1", "info")]
                    if(verbose) message(length(indx), " SNPs were removed from the analysis for sd = 0")
                    }
        }
        
        write.table(rm.snps, paste0(outfile, ".snps_removed"),
                    row.names = FALSE, sep="\t", quote = FALSE)
        if(verbose) message("List of removed SNPs are saved to ", paste0(outfile, ".snps_removed"))
        
        ## STARTING ANALYSIS PORTION
        
        pheno.file <- as.matrix(covfile[,c(time.to.event, event, covariates)])
        
        # covariates are defined in pheno.file
        ok.covs <- colnames(pheno.file)[colnames(pheno.file) %in% covariates]
        if (verbose) message("Covariates included in the models are: ", paste(ok.covs, collapse=", "))
        
        if (!typeof(pheno.file) %in% c("numeric", "double", "integer") ) {
            stop("Provided covariates must be numeric! e.g. categorical variables should be recoded as indicator or dummy variables.")
        }
        
        # define Ns
        n.sample <- nrow(pheno.file)
        n.event <- sum(as.numeric(pheno.file[,event]))
        if (verbose) message(n.sample, " samples are included in the analysis")
        
        # build coxph.fit parameters
        params <- coxParam(pheno.file, time.to.event, event, covariates, sample.ids)
        
        if (verbose) message("Analysis started on ", format(Sys.time(), "%Y-%m-%d"), " at ", format(Sys.time(), "%H:%M:%S"))
        
        # create cluster, also create option to input number of cores
        cl <- makeForkCluster(getOption("gwasurvivr.cores", detectCores()))
        on.exit(stopCluster(cl))
    
        snp.out <- t(parApply(cl=cl, X=genotypes, MARGIN=1, FUN=survFit, params))
    
        
        # calculate z-score
        z <- snp.out[,1]/snp.out[,2]
        # calculate p-value
        pval <- 2*pnorm(abs(z), lower.tail=FALSE)
        # calculate hazard ratio
        hr <- exp(snp.out[,1])
        # confidence interval HR
        lowerCI <- exp(snp.out[,1]-1.96*snp.out[,2])
        upperCI <- exp(snp.out[,1]+1.96*snp.out[,2])
        
        # putting everything back together
        sres <- cbind(snp.out, hr, lowerCI, upperCI, z, pval, n.sample, n.event)
        colnames(sres) <- c("coef", "se.coef", "exp.coef", "lb", "ub", "z", "p.value", "n", "nevents")
        rownames(sres) <- NULL # remove rownames so we don't have a duplicated rownames issue
        
        res <- data.frame(cbind(snp,sres))
        write.table(res, file=paste0(outfile, ".surv.results"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
        
        if (verbose) message("Analysis completed on ", format(Sys.time(), "%Y-%m-%d"), " at ", format(Sys.time(), "%H:%M:%S"))
        # on.exit(unlink(tempdir(), recursive = TRUE))
}




