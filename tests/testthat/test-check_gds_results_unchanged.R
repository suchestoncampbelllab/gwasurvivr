local_edition(3) 
results_file_name <- "example_results"

gdsfile <- system.file(package="gwasurvivr",
                       "extdata",
                       "gds_example.gds")
covariate.file <- system.file(package="gwasurvivr", 
                              "extdata",
                              "simulated_pheno.txt")
covariate.file <- read.table(covariate.file,
                             sep=" ",
                             header=TRUE,
                             stringsAsFactors = FALSE)
covariate.file$SexFemale <- ifelse(covariate.file$sex=="female", 1L, 0L)
sample.ids <- covariate.file[covariate.file$group=="experimental",]$ID_2
gdsCoxSurv(gdsfile=gdsfile,
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
           flip.dosage=TRUE,
           verbose=TRUE,
           clusterObj=NULL) 


if(file.exists(paste0(results_file_name, ".coxph"))){
  results_file <- read.table(paste0(results_file_name, ".coxph"), sep="\t", header=TRUE, stringsAsFactors = FALSE)
} else {
  results_file <- NULL
}

test_that("check gds example results not changed", {
  expect_snapshot_output(results_file, cran = F)
})

if(file.exists(paste0(results_file_name, ".coxph"))){
  file.remove(paste0(results_file_name, ".coxph"))
}

if(file.exists(paste0(results_file_name, ".snps_removed"))){
  file.remove(paste0(results_file_name, ".snps_removed"))
}