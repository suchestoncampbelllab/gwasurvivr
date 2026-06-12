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

if(file.exists(paste0(results_file_name, ".coxph"))){
  results_file <- read.table(paste0(results_file_name, ".coxph"), 
                             sep="\t", 
                             header=TRUE, 
                             stringsAsFactors = FALSE)
} else {
  results_file <- NULL
}

# Round numeric columns so the snapshot is stable across R / survival
# versions (exact low-order digits drift, especially for degenerate
# perfectly-separated fits). Numeric correctness is asserted independently,
# with tolerance, in test-stat-oracle.R.
if (!is.null(results_file)) {
  # Only round non-integer columns (statistics); leave integers such as POS,
  # CHR, N and N.EVENT exact.
  .frac <- vapply(results_file, function(col)
    is.numeric(col) && any(col %% 1 != 0, na.rm = TRUE), logical(1))
  results_file[.frac] <- lapply(results_file[.frac], signif, digits = 4)
}

test_that("check sanger example results not changed", {
  expect_snapshot_output(results_file, cran = F)
})

if(file.exists(paste0(results_file_name, ".coxph"))){
  file.remove(paste0(results_file_name, ".coxph"))
}

if(file.exists(paste0(results_file_name, ".snps_removed"))){
  file.remove(paste0(results_file_name, ".snps_removed"))
}
