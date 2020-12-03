local_edition(3) 
results_file_name <- "example_results"


impute.file <- system.file(package="gwasurvivr",
                           "extdata",
                           "impute_example.impute2.gz")
sample.file <- system.file(package="gwasurvivr",
                           "extdata", 
                           "impute_example.impute2_sample")
covariate.file <- system.file(package="gwasurvivr", 
                              "extdata",
                              "simulated_pheno.txt")
covariate.file <- read.table(covariate.file,
                             sep=" ",
                             header=TRUE,
                             stringsAsFactors = FALSE)
covariate.file$SexFemale <- ifelse(covariate.file$sex=="female", 1L, 0L)
sample.ids <- covariate.file[covariate.file$group=="experimental",]$ID_2

impute2CoxSurv(impute.file=impute.file,
              sample.file=sample.file,
              chr=14,
              covariate.file=covariate.file,
              id.column="ID_2",
              sample.ids=sample.ids,
              time.to.event="time",
              event="event",
              covariates=c("age", "SexFemale", "DrugTxYes"),
              inter.term=NULL,
              print.covs="only",
              out.file=results_file_name,
              chunk.size=50,
              maf.filter=0.005,
              exclude.snps=NULL,
              flip.dosage=TRUE,
              verbose=TRUE,
              clusterObj=NULL,
              keepGDS=FALSE)  

results_file <- read.table(paste0(results_file_name, ".coxph"), sep="\t", header=TRUE, stringsAsFactors = FALSE)

test_that("check impute2 results not changed from example", {
  skip_on_cran()
  expect_snapshot_output(results_file)
})

file.remove(paste0(results_file_name, ".coxph"))
file.remove(paste0(results_file_name, ".snps_removed"))
