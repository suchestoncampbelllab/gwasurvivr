library(VariantAnnotation)
library(survival)
library(data.table)
library(dplyr); library(tidyr)

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



vcfCoxSurv <- function(vcf.file, chunk.size, pheno.file, time, event, 
                       covariates, sample.ids, output.name){

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
                fit <- coxph(formula=as.formula(formula), data=pheno.file)
                res <- summary(fit)
                
                # following deals with snps that has no variance 
                # e.g. (have the same genotpye probability accross samples)
                fit <- tryCatch(coxph(as.formula(formula), data=pheno.file),
                                 warning=function(warn) NA,
                                 error=function(err) NA
                )
                
                if(anyNA(fit)){
                        rep(NA, 9)
                }else{
                        m <- summary(fit)
                        c(m$coef[1,], m$conf.int[1,-c(1:2)] ,n=m$n, nevents=m$nevent)
                }
                                

        }
           
        # define vcf.file and chunks, open vcf file
        vcf <- VcfFile(vcf.file, yieldSize=chunk.size)
        open(vcf)
        chunk_start <- 0
        chunk_end <- chunk.size
        
        write.table(data.frame(t(c("coef", "exp.coef", "se.coef", "z", "p.value",
                                       "lower.CI95", "upper.CI95", "n", "n.event"))), 
                            paste0(output.name, ".coxph"),
                            append = F, 
                            row.names = F,
                            col.names = F,
                            quote = F,
                            sep="\t")
        
        # get genotype probabilities by chunks
        # apply the survival function and save output
        repeat{ 
                # read in just dosage data from Vcf file
                data <- readVcf(vcf, param=ScanVcfParam(geno="DS"))
                # read dosage data from collapsed vcf, subset for defined ids
                genotype <- geno(data)$DS[, sample.ids]

                # message user
                message("Analyzing chunk ", chunk_start, "-", chunk_end)

                # apply survival function
                snp.out <- t(apply(genotype, 1, survFit))
                # change colnames to be more programming friendly
                colnames(snp.out) <- c("coef", "exp.coef", "se.coef", "z", "p.value",
                                       "lower.CI95", "upper.CI95", "n","n.event")
                write.table(data.frame(snp.out), 
                            paste0(output.name, ".coxph"),
                            append = T, 
                            row.names = F,
                            col.names = F,
                            quote = F,
                            sep="\t")

                chunk_start <- chunk_start+chunk.size
                chunk_end <- chunk_end+chunk.size
                

                # if(chunk_end==10000){
                #         break
                # }



                if(nrow(data)==0){
                        break
                }
        }
        
        close(vcf.file)
}




