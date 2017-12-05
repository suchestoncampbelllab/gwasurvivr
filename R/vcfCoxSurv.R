#' Fit cox survival to all variants in a VCF file
#' 
#' Performs survival analysis using Cox proportional hazard models on imputed genetic data stored in compressed VCF files 
#' 
#' @param vcf.file character(1) path to VCF file.
#' @param chunk.size integer(1) number of variants to process per thread
#' @param pheno.file matrix(1) comprising phenotype data
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
#' @importFrom parallel detectCores
#' @importFrom parallel parApply
#' @export

vcfCoxSurv <- function(vcf.file, # character, path to vcf file
                       chunk.size, # integer, defines the size of the chunk
                       pheno.file, # this needs to be disussed either a matrix or file path to a file with specific format
                       time, # character, column defining time
                       event, # character, column defining event
                       covariates, # character vector, columns defining covariates
                       sample.ids, # character vector, list of samples that will be analyzed, could also be a file path?
                       output.name # character, name of the output file
){
    
    # subset phenotype file for sample ids
    pheno.file <- pheno.file[sample.ids, ]
    
    ### define survFit
    survFit <- function(input.genotype){
        
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
        
        # maybe change init?
        INIT <- NULL
        
        ## only thing that's changing
        X <- cbind(input.genotype,pheno.file[,covariates])
        ## run fit
        fit <- coxph.fit(
            X, Y, STRATA, OFFSET, INIT, CONTROL, WEIGHTS, METHOD, ROWNAMES
        )
        
        ## extract statistics
        coef <- fit$coefficients[1]
        se <- sqrt(diag(fit$var)[1])
        cbind(coef, se)
    }
    
    vcf <- VcfFile(vcf.file, yieldSize=chunk.size)
    open(vcf)
    
    chunk_start <- 0
    chunk_end <- chunk.size
    
    write.table(t(c("snp",
                    "coef",
                    # "exp.coef",
                    "se.coef",
                    "z",
                    "p.value"
                    # "lower.CI95",
                    # "upper.CI95",
                    # "n",
                    # "n.event"
    )), 
    paste0(output.name, ".coxph"),
    append = F, 
    row.names = F,
    col.names = F,
    quote = F,
    sep="\t")
    
    # get genotype probabilities by chunks
    # apply the survival function and save output
    
    # for a single machine
    cl <- makeForkCluster(detectCores())
    
    repeat{ 
        # read in just dosage data from Vcf file
        data <- readVcf(vcf, param=ScanVcfParam(geno="DS"))
        
        if(nrow(data)==0){
            break
        }
        # read dosage data from collapsed vcf, subset for defined ids
        genotype <- geno(data)$DS[, sample.ids, drop=F]
        
        ## Add step to check for sd of the snps with:
        indx <- which(matrixStats::rowSds(genotype) == 0)
        
        ## Save list of snps that have a MAF=0
        snps_maf0 <- rownames(genotype)[indx]
        write.table(snps_maf0, 
                    file= paste0(output.name, ".MAF0snps"),
                    append = T, 
                    row.names = T,
                    col.names = F,
                    quote = F,
                    sep="\t")
        
        ## Remove MAF=0 snps
        genotype <- genotype[-indx,]
        
        # message user
        message("Analyzing chunk ", chunk_start, "-", chunk_end)
        
        # apply survival function
        snp.out <- t(parApply(cl=cl, X=genotype, MARGIN=1, FUN=survFit))
        
        
        z <- snp.out[,1]/snp.out[,2]
        pval <- 2*pnorm(abs(z), lower.tail=F)
        snp.out <- cbind(snp.out,z,pval)
        
        write.table(snp.out, 
                    paste0(output.name, ".coxph"),
                    append = T, 
                    row.names = T,
                    col.names = F,
                    quote = F,
                    sep="\t")
        
        chunk_start <- chunk_start+chunk.size
        chunk_end <- chunk_end+chunk.size
        
    }
    
    close(vcf)
}




