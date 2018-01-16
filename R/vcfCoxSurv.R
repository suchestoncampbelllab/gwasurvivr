#' Fit cox survival to all variants in a VCF file
#' 
#' Performs survival analysis using Cox proportional hazard models on imputed genetic data stored in compressed VCF files 
#' 
#' @param vcf.file character(1) path to VCF file.
#' @param vcf.file character(1) path to corresponding info file.
#' @param chunk.size integer(1) number of variants to process per thread
#' @param pheno.file numeric matrix(1) comprising phenotype data. 
#' @param time character(1) string that matches time column name in pheno.file
#' @param event character(1) string that matches event column name in pheno.file
#' @param covariates character vector with matching column names in pheno.file of covariates of interest
#' @param sample.ids character vector with sample ids to include in analysis
#' @param output.name character(1) string with output name
#' 
#' @return 
#' Saves text file directly to disk that contains survival analysis results
#' 
#' @examples 
#' set.seed(222)
#' vcf.file <- system.file(package="SurvivR", "extdata", "chr21.25000005-25500000.vcf.gz")
#' fl <- system.file(package="SurvivR", "extdata", "pheno_file.txt") 
#' pheno.file <- read.table(fl, sep="\t", header = T)
#' time <- "intxsurv_1Y"
#' event <- "dead_1Y" 
#' covariates <- c("distatD", "age") 
#' sample.ids = sample(rownames(pheno.file), size=190)
#' pheno.file$distatD <- as.integer(pheno.file$distatD) - 1L
#' pheno.file <- as.matrix(pheno.file)
#' output.name <- tempfile()
#' vcfCoxSurv(vcf.file, chunk.size, pheno.file, time, event, covariates, sample.ids, output.name)
#'  
#' @importFrom survival Surv coxph.fit
#' @import parallel
#' @import VariantAnnotation
#' 
#' @export

vcfCoxSurv <- function(vcf.file, # character, path to vcf file
                       info.file, # character, path to corresponding info file
                       MAF, # double, defining the MAF threshold. Any SNP with lower MAF will be excluded from analysis.
                       Rsq, # double, defining the Rsq threshold.Any SNP with lower Rsq will be excluded from analysis.
                       chunk.size, # integer, defines the size of the chunk
                       pheno.file, # this needs to be disussed either a matrix or file path to a file with specific format
                       time, # character, column defining time
                       event, # character, column defining event
                       covariates, # character vector, columns defining covariates
                       sample.ids, # character vector, list of samples that will be analyzed, could also be a file path?
                       output.name # character, name of the output file
                       
){
    
    message("Analysis started on ", format(Sys.time(), "%Y-%m-%d"), " at ", format(Sys.time(), "%H:%M:%S"))
   
    
    # subset phenotype file for sample ids
    pheno.file <- pheno.file[sample.ids, ]
    message("Analysis running for ", nrow(pheno.file), " samples.")
    
    # covariates are defined in pheno.file
    ok.covs <- colnames(pheno.file)[colnames(pheno.file) %in% covariates]
    message("Covariates included in the models are: ", paste(ok.covs, collapse=", "))
    message("If your covariates of interest are not included in the model\nplease stop the analysis and make sure user defined covariates\nmatch the column names in the pheno.file")
    
    
    ### define arguments for the survfit
    
    ### building arguments for coxph.fit ###
    Y <- Surv(time=pheno.file[,time], event=pheno.file[,event])
    rownames(Y) <- as.character(seq_len(nrow(Y)))
    STRATA <- NULL
    CONTROL <- structure(
        list(
            eps = 1e-09,
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
                          Y, STRATA, OFFSET, INIT, CONTROL, WEIGHTS, METHOD, ROWNAMES
    )
    myINIT <- c(0,  init.fit$coefficients)
    
    ### define survFit
    survFit <- function(input.genotype){
        
        ## only thing that's changing
        X <- cbind(input.genotype,pheno.file[,covariates])
        ## run fit
        fit <- coxph.fit(
            X, Y, STRATA, OFFSET, myINIT, CONTROL, WEIGHTS, METHOD, ROWNAMES
        )
        
        ## extract statistics
        coef <- fit$coefficients[1]
        se <- sqrt(diag(fit$var)[1])
        iter <- fit$iter
        
        cbind(coef, se, iter)
    }
    
    vcf <- VcfFile(vcf.file, yieldSize=chunk.size)
    open(vcf)
    
    chunk.start <- 0
    # chunk_end <- chunk.size
    snps_removed <- 0
    
    ## Save header for the cox.surv output
    #
    # InputName rsid Chr Pos EA NonEA CoefValue HR SE LowerCI UpperCI 
    # Waldpv LRTpv ModLRTpv EAF MAF Infoscore
    
    
    write.table(t(c("snp",
                    "genotyped",
                    "ref_0_allele", 
                    "alt_1_allele", 
                    "alt_allele_frq", 
                    "maf",
                    "hr",
                    "lowerCI",
                    "upperCI",
                    "p.value",
                    "z",
                    "coef",
                    "se.coef",
                    "iter",
                    "AvgCall",
                    "Rsq",
                    "LooRsq",
                    "EmpR",
                    "EmpRsq",
                    "Dose0",
                    "Dose1"    
                    )), 
    paste0(output.name, ".coxph"),
    append = F, 
    row.names = F,
    col.names = F,
    quote = F,
    sep="\t")
    
    # for a single machine
    cl <- makeForkCluster(detectCores())
    
    # get genotype probabilities by chunks
    # apply the survival function and save output
    
    repeat{ 
        # read in just dosage data from Vcf file
        data <- readVcf(vcf, param=ScanVcfParam(geno="DS"))
        
        if(nrow(data)==0){
            break
        }
        
        # read info file
        snp.info <- fread(info.file, 
                          skip=chunk.start, 
                          nrows = chunk.size,
                          na.strings = "-")
        
        # read dosage data from collapsed vcf, subset for defined ids
        genotype <- geno(data)$DS[, sample.ids, drop=F]
        
        ## Check for sd of the snps, remove snps that doesn't meet
        ## MAF and Rsq thresholds
        indx <- sort(unique(c(which(matrixStats::rowSds(genotype) == 0),
                              which(snp.info[["MAF"]] < MAF & snp.info[["Rsq"]] < Rsq))))
        
        if(length(indx) != 0){          
        ## Save list of snps that have a MAF=0
        snps_maf0 <- rownames(genotype)[indx]
        snps_removed <- snps_removed + length(snps_maf0)
        write.table(snps_maf0, 
                    file= paste0(output.name, ".MAF0snps"),
                    append = TRUE, 
                    row.names = FALSE,
                    col.names = FALSE,
                    quote = FALSE,
                    sep="\t")
        
        ## Remove MAF=0 snps
        genotype <- genotype[-indx,]
        snp.info <- snp.info[-indx,]
        }
        
        # message user
        message("Analyzing chunk ", chunk.start, "-", chunk.start+chunk.size)
        
        # apply survival function
        snp.out <- t(parApply(cl=cl, X=genotype, MARGIN=1, FUN=survFit))
        # snp.out <- survFit(test.input.genotype)
        
        z <- snp.out[,1]/snp.out[,2]
        pval <- 2*pnorm(abs(z), lower.tail=F)
        hr <- exp(snp.out[,1])
        lowerCI <-exp(snp.out[,1]-1.96*snp.out[,2])
        upperCI <-exp(snp.out[,1]+1.96*snp.out[,2])
        cox.out <- cbind(snp.info[,c(1, #"SNP", 
                                     8, #"Genotyped",
                                     2, #"REF.0.", 
                                     3, #"ALT.1.", 
                                     4, #"ALT_Frq", 
                                     5  #"MAF"
                                     ), with=FALSE],
                         hr, lowerCI, upperCI, pval, z, 
                         snp.out, 
                         snp.info[,-c(1:5,8), with=F])
        
        write.table(cox.out, 
                    paste0(output.name, ".coxph"),
                    append = TRUE, 
                    row.names = FALSE,
                    col.names = FALSE,
                    quote = FALSE,
                    sep="\t")
        
        chunk.start <- chunk.start+chunk.size
        
    }
    close(vcf)
    message("Analysis completed on ", format(Sys.time(), "%Y-%m-%d"), " at ", format(Sys.time(), "%H:%M:%S"))
    message(length(snps_maf0), " SNPs were removed from the analysis for not meeting the threshold criteria.")
    message("List of removed SNPs can be found in ", paste0(output.name, ".MAF0snps"))
}




