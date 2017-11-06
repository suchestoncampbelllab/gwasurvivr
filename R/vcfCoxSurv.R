library(VariantAnnotation)
library(survival)
library(data.table)
library(dplyr)
library(tidyr)
library(microbenchmark)
library(parallel)

pdata <-fread("/projects/rpci/lsuchest/lsuchest/DBMT_PhenoData/DBMTpheno_EA_long_20171023.txt")


## write function
vcf.file="./chr21_sub/chr21.25000000-26000000.dose.vcf.recode.vcf.gz"
chunk.size=10000
time="intxsurv_1Y"
event="dead_1Y"
covariates=c("distatD", "age")
pheno.file = pdata %>% 
                mutate(sample.ids=paste0("SAMP", 1:nrow(.))) %>%
                dplyr::select(sample.ids, intxsurv_1Y, age, distatD, intxsurv_1Y, dead_1Y)
sample.ids = paste0("SAMP", sample(1:1000, size=200))
output.name="test_survivR/chr21"
info.file <- "chr21.info"


<<<<<<< HEAD:R/vcfCoxSurv.R
vcfCoxSurv <- function(vcf.file, chunk.size, info.file, pheno.file, time, event, covariates, sample.ids, output.name){
=======

vcfCoxSurv <- function(vcf.file,
                       chunk.size,
                       pheno.file,
                       time,
                       event, 
                       covariates,
                       sample.ids,
                       output.name){

>>>>>>> 5ce982623072949f0f454d06eee4fe7971e04272:vcfCoxSurv.R
        
        # subset phenotype file for sample ids
        pheno.file <- pheno.file[match(sample.ids, pheno.file$sample.ids), ]
        
        # define the formula for CoxPH
        formula <- paste0("Surv(time=",
                          time,
                          ", event=",
                          event,
                          ") ~ input.genotype + ",
                          paste(covariates, collapse=" + "))
        
        # assign survival function with the defined formula

        survFit <- function(input.genotype) {
                # make sure this function always returns something that always has same structure
                
                # following deals with snps that has no variance 
                # e.g. (have the same genotpye probability accross samples)
                
                # process data ahead of time to exclude rows that we know would fail
                # e.g. NA or no variability
                # filter data before hand before it enters function
                
                
                fit <-coxph(formula=as.formula(formula), data=pheno.file)
                res <- summary(fit)
                
                if(dim(res$coef)[1] == 1+length(covariates)){
                        out <- c(res$coef[1,], res$conf.int[1,-c(1:2)], n=res$n, n.events=res$nevent)
                }else{
                        out <- rep(NA, 9)
                }
                
                # fit <- tryCatch(coxph(as.formula(formula), data=pheno.file),
                #                  warning=function(warn) NA,
                #                  error=function(err) NA
                # )
                # 
                # if(anyNA(fit)){
                #         rep(NA, 9)
                # }else{
                #         m <- summary(fit)
                #         c(m$coef[1,], m$conf.int[1,-c(1:2)], n=m$n, nevents=m$nevent)
                # }


        }

        vcf <- VcfFile(vcf.file, yieldSize=chunk.size)
        open(vcf)
        
        chunk_start <- 0
        chunk_end <- chunk.size
        
        write.table(t(c("snp",
                        "coef",
                        "exp.coef",
                        "se.coef",
                        "z",
                        "p.value",
                        "lower.CI95",
                        "upper.CI95",
                        "n",
                        "n.event")), 
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
<<<<<<< HEAD:R/vcfCoxSurv.R
        
        # read in info file
        info <- read.table(info.file, header=T, na.strings="-", stringsAsFactors = F, sep="\t")
        info <- info %>%
                rename(REF="REF.0.",
                       ALT="ALT.1.")
=======
>>>>>>> 5ce982623072949f0f454d06eee4fe7971e04272:vcfCoxSurv.R
        
        microbenchmark(
        repeat{ 
                # read in just dosage data from Vcf file
                data <- readVcf(vcf, param=ScanVcfParam(geno="DS"))
                
                if(nrow(data)==0){
                        break
                }
                # read dosage data from collapsed vcf, subset for defined ids
<<<<<<< HEAD:R/vcfCoxSurv.R
                genotype <- geno(data)$DS[1:100, sample.ids, drop=F] 
                genotype <- genotype %>%
                        data.table(keep.rownames = T) %>%
                        rename(SNP=rn)
                
                write.table(genotype, "genotype.dosage", quote=F, row.names=F, sep="\t", col.names=T)
                
                info.data <- data.frame(cbind(rownames(data), data.frame(rowRanges(data)))) %>%
                        rename(SNP=rownames.data.) %>%
                        select(SNP, REF, ALT) %>%
                        mutate(SNP=as.character(SNP),
                               ALT=as.character(ALT))
                
                genotype.merged <- data.frame(info.data, data.table(genotype))
                
=======
                genotype <- geno(data)$DS[, sample.ids, drop=F]
>>>>>>> 5ce982623072949f0f454d06eee4fe7971e04272:vcfCoxSurv.R
                
                # message user
                message("Analyzing chunk ", chunk_start, "-", chunk_end)
                
                # apply survival function
                
                snp.out <- t(parApply(cl=cl, X=genotype, MARGIN=1, FUN=survFit))
                

                # change colnames to be more programming friendly
                colnames(snp.out) <- c("coef",
                                       "exp.coef",
                                       "se.coef",
                                       "z",
                                       "p.value",
                                       "lower.CI95",
                                       "upper.CI95",
                                       "n",
                                       "n.event")
                
                snp.out <- cbind(snp=rownames(snp.out), snp.out)
                
                write.table(snp.out, 
                            paste0(output.name, ".coxph"),
                            append = T, 
                            row.names = F,
                            col.names = F,
                            quote = F,
                            sep="\t")
                
                chunk_start <- chunk_start+chunk.size
                chunk_end <- chunk_end+chunk.size
                
        },
        times = 1)

        close(vcf)
}




