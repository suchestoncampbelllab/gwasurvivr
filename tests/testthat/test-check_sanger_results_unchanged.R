local_edition(3) 
results_file_name <- "example_results"


vcf.file <- system.file(package="gwasurvivr",
                      "extdata",
                      "sanger.pbwt_reference_impute.vcf.gz")
pheno.fl <- system.file(package="gwasurvivr",
                       "extdata",
                    "simulated_pheno.txt")
pheno.file <- read.table(pheno.fl,
                        sep=" ",
                        header=TRUE,
                        stringsAsFactors = FALSE)
pheno.file$SexFemale <- ifelse(pheno.file$sex=="female", 1L, 0L)
sample.ids <- pheno.file[pheno.file$group=="experimental",]$ID_2
sangerCoxSurv(vcf.file=vcf.file,
             covariate.file=pheno.file,
             id.column="ID_2",
             sample.ids=sample.ids,
             time.to.event="time",
             event="event",
             covariates=c("age", "SexFemale", "DrugTxYes"),
             inter.term=NULL,
             print.covs="only",
             out.file=results_file_name,
             info.filter=0.3,
             maf.filter=0.005,
             chunk.size=50,
             verbose=TRUE,
             clusterObj=NULL)

results_file <- read.table(paste0(results_file_name, ".coxph"), sep="\t", header=TRUE, stringsAsFactors = FALSE)

test_that("check sanger example results not changed", {
  skip_on_cran()
  expect_snapshot_output(results_file)
})

file.remove(paste0(results_file_name, ".coxph"))
file.remove(paste0(results_file_name, ".snps_removed"))

