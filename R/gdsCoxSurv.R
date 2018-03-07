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
#' @param maf.filer numeric(1) to filter minor allele frequency (MAF)
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
        
        gdsfile <- paste0(outfile, ".gds")
        snpfile <- paste0(outfile, ".snp.rdata")
        scanfile <- paste0(outfile, ".scan.rdata")
        
        # see if files exist already ... if not convert to GDS ... still need to test if this works if files dont exist
        if(!file.exists(gdsfile) | !file.exists(snpfile) | !file.exists(scanfile)){
                GWASTools::imputedDosageFile(input.files=c(impute.file, sample.file),
                                             filename=tempfile(pattern=outfile, fileext = ".gds"),
                                             chromosome=as.numeric(chromosome),
                                             input.type="IMPUTE2",
                                             input.dosage=FALSE,
                                             file.type="gds",
                                             snp.annot.filename = tempfile(pattern=outfile, fileext = ".snp.rdata"),
                                             scan.annot.filename = tempfile(pattern=outfile, fileext = ".scan.rdata"))
                # read in files from temp directory
                gdsfile <- paste0(tempdir(), "/", dir(path=tempdir(), pattern=".gds"))
                snpfile <- paste0(tempdir(), "/", dir(path=tempdir(), pattern=".snp.rdata"))
                scanfile <- paste0(tempdir(), "/", dir(path=tempdir(), pattern=".scan.rdata"))
        }
        
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
        colnames(snp)[colnames(snp)=="snpID"] <- "snp.index"
        colnames(snp)[colnames(snp)=="rsID"] <- "rsid"
        colnames(snp)[colnames(snp)=="snp"] <- "snpid"
        # grab sample file data
        scanAnn <- getAnnotation(getScanAnnotation(genoData))
        
        # change order into what we want final outside to be
        colnames(snp)[colnames(snp)=="chromosome"] <- "chr"
        colnames(snp)[colnames(snp)=="alleleA"] <- "A0"
        colnames(snp)[colnames(snp)=="alleleB"] <- "A1"
        
        
        dimnames(genotypes) <- list(paste(snp$snpid, snp$rsid, sep=";"), scanAnn$ID_2)
        
        # flip dosage
        if(flip.dosage) genotypes <- 2 - genotypes
        
        
        # add covariates to scan file
        # user needs to input covfile 
        # covfile <- read.table(covfile,
        #                       header=TRUE,
        #                       stringsAsFactors=FALSE,
        #                       sep="\t")
        # 
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
        snp <- snp[,c("chr",
                      "position",
                      "snpid",
                      "rsid",
                      "A0",
                      "A1", 
                      "exp_freq_A1", 
                      "info")]
        
        if(!is.null(maf.filter)){
                maf.idx <- snp$exp_freq_A1>maf.filter & snp$exp_freq_A1<(1-maf.filter)
                rm.maf <- snp[!maf.idx, c("snpid", "rsid", "exp_freq_A1")]
                snp <- snp[maf.idx,]
                genotypes <- genotypes[maf.idx,]
        }
        
        ## STARTING ANALYSIS PORTION
        if (verbose) message("Analysis started on ", format(Sys.time(), "%Y-%m-%d"), " at ", format(Sys.time(), "%H:%M:%S"))
        
        pheno.file <- as.matrix(scanAnn[,c(time.to.event, event, covariates)])
        
        # define Ns
        n.sample <- nrow(pheno.file)
        n.event <- sum(as.numeric(pheno.file[,event]))
        if (verbose) message("Analysis running for ", n.sample, " samples.")
        
        # covariates are defined in pheno.file
        ok.covs <- colnames(pheno.file)[colnames(pheno.file) %in% covariates]
        if (verbose) message("Covariates included in the models are: ", paste(ok.covs, collapse=", "))

        ### Check of snps with MAF = 0  ###  
        indx <- sort(unique(c(which(matrixStats::rowSds(genotypes) == 0))))
        
        ### which() may produce integer(0), have a check
        if(length(indx) != 0){
                if(verbose) message(length(indx), " SNPs were removed from the analysis for not meeting the threshold criteria.")
                ## Remove MAF=0 snps
                genotypes <- genotypes[-indx,]
                # remove from snp list too that we will merge back later
                snp <- snp[-indx,]
                # save list of snps that were removed
                rm.snps <- list(removedSNPs=snp$rsid[indx])
        }
        
        if (verbose) message("Running survival models...")
        
        # build coxph.fit parameters
        params <- .coxParam(pheno.file, time.to.event, event, covariates, sample.ids)
        
        snp.out <- t(apply(genotypes, 1, .survFit, params)) 

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
        write.table(res, file=paste0(outfile, ".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
        
        if (verbose) message("Analysis completed on ", format(Sys.time(), "%Y-%m-%d"), " at ", format(Sys.time(), "%H:%M:%S"))
        #return(res)
}




