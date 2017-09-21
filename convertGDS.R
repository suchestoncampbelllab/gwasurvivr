library(GWASTools)
# so we can submit using R CMD BATCH
args <- commandArgs(TRUE)

# run GWAS tools to convert imputed data to GDS
# arguments: 1. chunk name, 2. chromosome number

convertImputeGds <- function(args){
        path <- "/projects/rpci/lsuchest/lsuchest/Rserve/ImputeData/var/db/gwas/imputed_data/BMT093013_forImpute/"
        gdsfile <- paste0(args[1], ".gds")
        snpfile <- paste0(args[1], ".snp.rdata")
        scanfile <- paste0(args[1], ".scan.rdata")
        imputedDosageFile(input.files=c(paste0(path, args[1], ".impute2"), paste0(path, args[1],".impute2_samples")),
                          filename=gdsfile,
                          chromosome=as.numeric(args[2]),
                          input.type="IMPUTE2",
                          input.dosage=F,
                          file.type="gds",
                          snp.annot.filename = snpfile,
                          scan.annot.filename = scanfile)
}

convertImputeGds(args)