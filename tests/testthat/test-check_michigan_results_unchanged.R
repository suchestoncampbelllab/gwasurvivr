local_edition(3) 
results_file_name <- "example_results"


vcf.file <- system.file(package="gwasurvivr",
                         "extdata",
                      "michigan.chr14.dose.vcf.gz")
pheno.fl <- system.file(package="gwasurvivr",
                      "extdata",
                   "simulated_pheno.txt")
pheno.file <- read.table(pheno.fl, 
                       sep=" ",
                       header=TRUE,
                       stringsAsFactors = FALSE)
pheno.file$SexFemale <- ifelse(pheno.file$sex=="female", 1L, 0L)
sample.ids <- pheno.file[pheno.file$group=="experimental",]$ID_2

michiganCoxSurv(vcf.file=vcf.file,
            covariate.file=pheno.file,
            id.column="ID_2",
            sample.ids=sample.ids,
            time.to.event="time",
            event="event",
            covariates=c("age", "SexFemale", "DrugTxYes"),
            inter.term=NULL,
            print.covs="only",
            out.file=results_file_name,
            r2.filter=0.3,
            maf.filter=0.005,
            chunk.size=50,
            verbose=FALSE,
            clusterObj=NULL)

if(file.exists(paste0(results_file_name, ".coxph"))){
  results_file <- read.table(paste0(results_file_name, ".coxph"), sep="\t", header=TRUE, stringsAsFactors = FALSE)
} else {
  results_file <- NULL
}

test_that("check michigan example results not changed", {
   expect_snapshot_output(results_file, cran = F)
})

if(file.exists(paste0(results_file_name, ".coxph"))){
  file.remove(paste0(results_file_name, ".coxph"))
}

if(file.exists(paste0(results_file_name, ".snps_removed"))){
  file.remove(paste0(results_file_name, ".snps_removed"))
}

