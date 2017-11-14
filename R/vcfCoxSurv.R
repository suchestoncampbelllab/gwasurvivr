vcfCoxSurv <- function(vcf.file, chunk.size, pheno.file, time, event, covariates, sample.ids, output.name){
        library(VariantAnnotation)
        library(survival)
        library(data.table)
        library(dplyr)
        library(tidyr)
        library(broom)
        library(microbenchmark)
        library(parallel)

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
                
                # number of NAs should be changed according to stats outputted
                if(sd(input.genotype)==0) {rep(NA, 8)} else { 
                    fit <-coxph(formula=as.formula(formula), data=pheno.file)
                    res <- summary(fit)
                    tidy(fit) %>% 
                        filter(!str_detect(term, paste0(covariates, collapse="|"))) %>%
                        mutate(n=res$n, n.event=res$nevent) %>%
                        select(-term) %>% as.numeric()
                }
        
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

        
        # read in info file
        # info <- read.table(info.file, header=T, na.strings="-", stringsAsFactors = F, sep="\t")
        # info <- info %>%
        #         rename(REF="REF.0.",
        #                ALT="ALT.1.")

        
        microbenchmark(
        repeat{ 
                # read in just dosage data from Vcf file
                data <- readVcf(vcf, param=ScanVcfParam(geno="DS"))
                
                if(nrow(data)==0){
                        break
                }
                # read dosage data from collapsed vcf, subset for defined ids
                genotype <- geno(data)$DS[, sample.ids, drop=F]
                
                ## Add step to check for sd of the snps with:
                # indx <- which(matrixStats::rowSds(genotype) == 0)
                
                ## List of snps tha have a MAF=0
                # snps_maf0 <- rownames(genotyoe)[indx]
                
                ## Remove MAF=0 snps
                # genotype <- genotype[-indx,]
                
                ## then save a file 
                
                # genotype <- genotype %>%
                #         data.table(keep.rownames = T) %>%
                #         rename(SNP=rn)
                
                # write.table(genotype, "genotype.dosage", quote=F, row.names=F, sep="\t", col.names=T)
                
                # info.data <- data.frame(cbind(rownames(data), 
                #                               data.frame(rowRanges(data)))) %>%
                #         rename(SNP=rownames.data.) %>%
                #         select(SNP, REF, ALT) %>%
                #         mutate(SNP=as.character(SNP),
                #                ALT=as.character(ALT)) 
                # 
                # genotype.merged <- data.frame(info.data, data.table(genotype))
                # 
                # genotype <- geno(data)$DS[, sample.ids, drop=F]
                
                # message user
                message("Analyzing chunk ", chunk_start, "-", chunk_end)
                
                # apply survival function
                snp.out <- t(parApply(cl=cl, X=genotype, MARGIN=1, FUN=survFit))
                
                
                
                
                # change colnames to be more programming friendly
                # colnames(snp.out) <- c("coef",
                #                        # "exp.coef",
                #                        "se.coef",
                #                        "z",
                #                        "p.value",
                #                        "lower.CI95",
                #                        "upper.CI95",
                #                        "n",
                #                        "n.event")
                
                # snp.out <- cbind(snp=rownames(snp.out), snp.out)
                
                write.table(snp.out, 
                            paste0(output.name, ".coxph"),
                            append = T, 
                            row.names = T,
                            col.names = F,
                            quote = F,
                            sep="\t")
                
                chunk_start <- chunk_start+chunk.size
                chunk_end <- chunk_end+chunk.size
                
        },
        times = 1)

        close(vcf)
}




