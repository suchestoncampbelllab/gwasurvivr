# Regression tests for specific user-reported bugs.
# Each test names the GitHub issue it locks down.
local_edition(3)
options(gwasurvivr.cores = 2)

ext <- function(f) system.file(package = "gwasurvivr", "extdata", f)

pheno <- read.table(ext("simulated_pheno.txt"), sep = " ", header = TRUE,
                    stringsAsFactors = FALSE)
pheno$SexFemale <- ifelse(pheno$sex == "female", 1L, 0L)
exp.ids <- pheno[pheno$group == "experimental", ]$ID_2

vcf.file <- ext("michigan.chr14.dose.vcf.gz")
bed.file <- ext("plink_example.bed")

tmp <- function(name) file.path(tempdir(), name)
cleanup <- function(stub) {
    for (ext in c(".coxph", ".snps_removed")) {
        f <- paste0(stub, ext)
        if (file.exists(f)) file.remove(f)
    }
}

test_that("#6: analysis runs with no covariates (covariates = NULL)", {
    stub <- tmp("bug6")
    on.exit(cleanup(stub))
    expect_no_error(
        plinkCoxSurv(bed.file = bed.file, covariate.file = pheno,
                     id.column = "ID_2", sample.ids = exp.ids,
                     time.to.event = "time", event = "event",
                     covariates = NULL, inter.term = NULL, print.covs = "only",
                     out.file = stub, chunk.size = 50, maf.filter = 0.005,
                     exclude.snps = NULL, flip.dosage = TRUE, verbose = FALSE,
                     clusterObj = NULL)
    )
    out <- read.table(paste0(stub, ".coxph"), sep = "\t", header = TRUE)
    expect_gt(nrow(out), 0)
    # SNP-only model: no covariate columns, just the SNP effect statistics.
    expect_true(all(c("COEF", "SE.COEF", "HR", "PVALUE") %in% colnames(out)))
})

test_that("#25: 'SNPs analyzed' count is reported (not blank) for VCF path", {
    stub <- tmp("bug25")
    on.exit(cleanup(stub))
    msgs <- character()
    withCallingHandlers(
        michiganCoxSurv(vcf.file = vcf.file, covariate.file = pheno,
                        id.column = "ID_2", sample.ids = exp.ids,
                        time.to.event = "time", event = "event",
                        covariates = c("age", "SexFemale", "DrugTxYes"),
                        inter.term = NULL, print.covs = "only", out.file = stub,
                        maf.filter = 0.005, r2.filter = 0.3, chunk.size = 50,
                        verbose = TRUE, clusterObj = NULL),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        }
    )
    n.coxph <- nrow(read.table(paste0(stub, ".coxph"), sep = "\t", header = TRUE))
    analyzed <- grep("analyzed in total", msgs, value = TRUE)
    expect_length(analyzed, 1)
    # The reported number must match the actual number of rows written.
    reported <- as.integer(sub("\\D*(\\d+).*", "\\1", analyzed))
    expect_equal(reported, n.coxph)
})

test_that("#32: SNP*covariate interaction produces non-empty output", {
    stub <- tmp("bug32")
    on.exit(cleanup(stub))
    expect_no_error(
        michiganCoxSurv(vcf.file = vcf.file, covariate.file = pheno,
                        id.column = "ID_2", sample.ids = exp.ids,
                        time.to.event = "time", event = "event",
                        covariates = c("age", "SexFemale", "DrugTxYes"),
                        inter.term = "DrugTxYes", print.covs = "some",
                        out.file = stub, maf.filter = 0.005, r2.filter = 0.3,
                        chunk.size = 50, verbose = FALSE, clusterObj = NULL)
    )
    out <- read.table(paste0(stub, ".coxph"), sep = "\t", header = TRUE)
    expect_gt(nrow(out), 0)
    expect_true(any(grepl("INTER.TERM", colnames(out))))
})

test_that("#36: snps_removed lists the right SNPs with real (not bogus) MAF", {
    stub <- tmp("bug36")
    on.exit(cleanup(stub))
    # MAF filter 0.35 drops the SNP with SAMP_MAF 0.3056, keeps the 0.3889 SNPs.
    plinkCoxSurv(bed.file = bed.file, covariate.file = pheno, id.column = "ID_2",
                 sample.ids = exp.ids, time.to.event = "time", event = "event",
                 covariates = c("age", "SexFemale"), inter.term = NULL,
                 print.covs = "only", out.file = stub, chunk.size = 50,
                 maf.filter = 0.35, exclude.snps = NULL, flip.dosage = TRUE,
                 verbose = FALSE, clusterObj = NULL)
    removed <- read.table(paste0(stub, ".snps_removed"), sep = "\t",
                          header = TRUE, colClasses = "character")
    kept <- read.table(paste0(stub, ".coxph"), sep = "\t", header = TRUE)
    # Exactly the sub-threshold SNP is removed, with its REAL MAF (~0.3056).
    expect_true("rs34919020" %in% removed$RSID)
    expect_false("rs34919020" %in% kept$RSID)
    expect_equal(as.numeric(removed$SAMP_MAF[removed$RSID == "rs34919020"]),
                 0.3056, tolerance = 1e-4)
    # Kept SNPs all clear the threshold.
    expect_true(all(kept$SAMP_MAF > 0.35))
})

test_that("#39: a grouped tibble covariate.file is accepted", {
    skip_if_not_installed("dplyr")
    stub <- tmp("bug39")
    on.exit(cleanup(stub))
    grouped <- dplyr::group_by(dplyr::as_tibble(pheno), group)
    expect_no_error(
        plinkCoxSurv(bed.file = bed.file, covariate.file = grouped,
                     id.column = "ID_2", sample.ids = exp.ids,
                     time.to.event = "time", event = "event",
                     covariates = c("age", "SexFemale"), inter.term = NULL,
                     print.covs = "only", out.file = stub, chunk.size = 50,
                     maf.filter = 0.005, exclude.snps = NULL, flip.dosage = TRUE,
                     verbose = FALSE, clusterObj = NULL)
    )
})

test_that("coxPheno gives a clear error when sample.ids do not match id.column", {
    stub <- tmp("bugids")
    on.exit(cleanup(stub))
    expect_error(
        plinkCoxSurv(bed.file = bed.file, covariate.file = pheno,
                     id.column = "ID_2", sample.ids = c("nope1", "nope2"),
                     time.to.event = "time", event = "event",
                     covariates = c("age"), inter.term = NULL,
                     print.covs = "only", out.file = stub, chunk.size = 50,
                     maf.filter = 0.005, exclude.snps = NULL, flip.dosage = TRUE,
                     verbose = FALSE, clusterObj = NULL),
        "sample.ids"
    )
})
