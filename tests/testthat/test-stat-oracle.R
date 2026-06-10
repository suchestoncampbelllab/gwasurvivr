# Statistical validation oracle.
#
# gwasurvivr's speed comes from calling survival's internal coxph.fit() directly
# with a warm-start INIT. This test proves that fast path returns the SAME
# estimates as an independent, standard survival::coxph() formula fit -- the
# reference "oracle" -- for every SNP, across covariate / no-covariate /
# interaction models. Dosages are re-read independently from the source VCF so
# both the genotype extraction AND the Cox fit are validated.
#
# Tolerance: gwasurvivr (coxph.fit, eps=1e-9, Efron ties) and coxph() (same
# fitter, default eps=1e-9, Efron ties) should agree to fp precision. We use a
# relative tolerance of 1e-5 for COEF/SE/HR and 1e-4 for p-values (p-values
# compress tail differences).
local_edition(3)
suppressMessages({
    library(survival)
    library(VariantAnnotation)
})
options(gwasurvivr.cores = 2)

ext <- function(f) system.file(package = "gwasurvivr", "extdata", f)
vcf.file <- ext("michigan.chr14.dose.vcf.gz")

pheno <- read.table(ext("simulated_pheno.txt"), sep = " ", header = TRUE,
                    stringsAsFactors = FALSE)
pheno$SexFemale <- ifelse(pheno$sex == "female", 1L, 0L)
exp.ids <- pheno[pheno$group == "experimental", ]$ID_2

# Independent dosage matrix (ALT-allele dose), SNP x sample.
DS <- geno(readVcf(vcf.file, param = ScanVcfParam(geno = "DS")))$DS

tmp <- function(name) file.path(tempdir(), name)
cleanup <- function(stub) {
    for (e in c(".coxph", ".snps_removed")) {
        f <- paste0(stub, e)
        if (file.exists(f)) file.remove(f)
    }
}

# Build the analysis data.frame gwasurvivr uses: experimental samples, complete
# cases on the modelled columns, in the covariate-file row order.
analysis_frame <- function(covars) {
    cols <- c("ID_2", "time", "event", covars)
    df <- pheno[pheno$ID_2 %in% exp.ids, cols, drop = FALSE]
    df[stats::complete.cases(df), , drop = FALSE]
}

# Oracle fit for one SNP on a given frame, dropping NA-dosage samples to mirror
# survFit's per-SNP NA handling.
oracle_fit <- function(rsid, df, rhs) {
    df$SNP <- as.numeric(DS[rsid, df$ID_2])
    df <- df[!is.na(df$SNP), , drop = FALSE]
    fit <- coxph(as.formula(paste("Surv(time, event) ~", rhs)),
                 data = df, ties = "efron")
    list(coef = unname(coef(fit)["SNP"]),
         se = unname(sqrt(diag(vcov(fit))["SNP"])))
}

# Compare gwasurvivr's COEF/SE/HR/PVALUE for `column` against the oracle.
expect_matches_oracle <- function(gw, rsid, ofit, hr_col = "HR",
                                  p_col = "PVALUE", coef_col = "COEF",
                                  se_col = "SE.COEF") {
    row <- gw[gw$RSID == rsid, ]
    z <- ofit$coef / ofit$se
    p <- 2 * pnorm(abs(z), lower.tail = FALSE)
    expect_equal(row[[coef_col]], ofit$coef, tolerance = 1e-5)
    expect_equal(row[[se_col]], ofit$se, tolerance = 1e-5)
    expect_equal(row[[hr_col]], exp(ofit$coef), tolerance = 1e-5)
    expect_equal(row[[p_col]], p, tolerance = 1e-4)
}

test_that("Michigan/DS with covariates matches coxph() per SNP", {
    stub <- tmp("oracle_cov")
    on.exit(cleanup(stub))
    covars <- c("age", "SexFemale", "DrugTxYes")
    michiganCoxSurv(vcf.file = vcf.file, covariate.file = pheno,
                    id.column = "ID_2", sample.ids = exp.ids,
                    time.to.event = "time", event = "event",
                    covariates = covars, inter.term = NULL, print.covs = "only",
                    out.file = stub, maf.filter = 0.005, r2.filter = 0.3,
                    chunk.size = 50, verbose = FALSE, clusterObj = NULL)
    gw <- read.table(paste0(stub, ".coxph"), sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE)
    df <- analysis_frame(covars)
    rhs <- paste(c("SNP", covars), collapse = " + ")
    for (rsid in gw$RSID) {
        expect_matches_oracle(gw, rsid, oracle_fit(rsid, df, rhs))
    }
})

test_that("Michigan/DS with no covariates matches coxph() per SNP", {
    stub <- tmp("oracle_nocov")
    on.exit(cleanup(stub))
    michiganCoxSurv(vcf.file = vcf.file, covariate.file = pheno,
                    id.column = "ID_2", sample.ids = exp.ids,
                    time.to.event = "time", event = "event",
                    covariates = NULL, inter.term = NULL, print.covs = "only",
                    out.file = stub, maf.filter = 0.005, r2.filter = 0.3,
                    chunk.size = 50, verbose = FALSE, clusterObj = NULL)
    gw <- read.table(paste0(stub, ".coxph"), sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE)
    df <- analysis_frame(character(0))
    for (rsid in gw$RSID) {
        expect_matches_oracle(gw, rsid, oracle_fit(rsid, df, "SNP"))
    }
})

test_that("Michigan/GT (hard-call) matches coxph() on derived dosage per SNP", {
    skip_if_not(file.exists(ext("michigan.chr14.GT.vcf.gz")))
    stub <- tmp("oracle_gt")
    on.exit(cleanup(stub))
    gtvcf <- ext("michigan.chr14.GT.vcf.gz")
    covars <- c("age", "SexFemale", "DrugTxYes")
    michiganCoxSurv(vcf.file = gtvcf, covariate.file = pheno, id.column = "ID_2",
                    sample.ids = exp.ids, time.to.event = "time", event = "event",
                    covariates = covars, inter.term = NULL, print.covs = "only",
                    out.file = stub, maf.filter = 0.005, r2.filter = 0.3,
                    chunk.size = 50, genotype.field = "GT", verbose = FALSE,
                    clusterObj = NULL)
    gw <- read.table(paste0(stub, ".coxph"), sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE)
    # Independent hard-call -> dosage from the GT field.
    GT <- geno(readVcf(gtvcf, param = ScanVcfParam(geno = "GT")))$GT
    lut <- c("0/0" = 0, "0/1" = 1, "1/0" = 1, "1/1" = 2)
    df <- analysis_frame(covars)
    rhs <- paste(c("SNP", covars), collapse = " + ")
    for (rsid in gw$RSID) {
        df2 <- df
        df2$SNP <- unname(lut[GT[rsid, df2$ID_2]])
        df2 <- df2[!is.na(df2$SNP), , drop = FALSE]
        fit <- coxph(as.formula(paste("Surv(time, event) ~", rhs)),
                     data = df2, ties = "efron")
        row <- gw[gw$RSID == rsid, ]
        expect_equal(row[["COEF"]], unname(coef(fit)["SNP"]), tolerance = 1e-5)
        expect_equal(row[["SE.COEF"]],
                     unname(sqrt(diag(vcov(fit))["SNP"])), tolerance = 1e-5)
    }
})

test_that("Michigan/DS interval (left-truncated) matches coxph() per SNP", {
    stub <- tmp("oracle_interval")
    on.exit(cleanup(stub))
    covars <- c("age", "SexFemale", "DrugTxYes")
    # Synthetic entry time (start < stop) to exercise counting-process fitting.
    ph2 <- pheno
    ph2$start <- ph2$time / 2
    michiganCoxSurv(vcf.file = vcf.file, covariate.file = ph2, id.column = "ID_2",
                    sample.ids = exp.ids, time.to.event = "time", event = "event",
                    covariates = covars, inter.term = NULL, print.covs = "only",
                    out.file = stub, maf.filter = 0.005, r2.filter = 0.3,
                    chunk.size = 50, start.time = "start", verbose = FALSE,
                    clusterObj = NULL)
    gw <- read.table(paste0(stub, ".coxph"), sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE)
    df <- ph2[ph2$ID_2 %in% exp.ids, c("ID_2", "start", "time", "event", covars)]
    df <- df[stats::complete.cases(df), , drop = FALSE]
    rhs <- paste(c("SNP", covars), collapse = " + ")
    for (rsid in gw$RSID) {
        df2 <- df
        df2$SNP <- as.numeric(DS[rsid, df2$ID_2])
        df2 <- df2[!is.na(df2$SNP), , drop = FALSE]
        # Counting-process oracle: Surv(start, stop, event).
        fit <- coxph(as.formula(paste("Surv(start, time, event) ~", rhs)),
                     data = df2, ties = "efron")
        row <- gw[gw$RSID == rsid, ]
        expect_equal(row[["COEF"]], unname(coef(fit)["SNP"]), tolerance = 1e-5)
        expect_equal(row[["SE.COEF"]],
                     unname(sqrt(diag(vcov(fit))["SNP"])), tolerance = 1e-5)
    }
})

test_that("Michigan/DS SNP*covariate interaction matches coxph() per SNP", {
    stub <- tmp("oracle_int")
    on.exit(cleanup(stub))
    covars <- c("age", "SexFemale", "DrugTxYes")
    michiganCoxSurv(vcf.file = vcf.file, covariate.file = pheno,
                    id.column = "ID_2", sample.ids = exp.ids,
                    time.to.event = "time", event = "event",
                    covariates = covars, inter.term = "DrugTxYes",
                    print.covs = "some", out.file = stub, maf.filter = 0.005,
                    r2.filter = 0.3, chunk.size = 50, verbose = FALSE,
                    clusterObj = NULL)
    gw <- read.table(paste0(stub, ".coxph"), sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE)
    df <- analysis_frame(covars)
    rhs <- paste(c("SNP", covars, "SNP:DrugTxYes"), collapse = " + ")
    for (rsid in gw$RSID) {
        df2 <- df
        df2$SNP <- as.numeric(DS[rsid, df2$ID_2])
        df2 <- df2[!is.na(df2$SNP), , drop = FALSE]
        fit <- coxph(as.formula(paste("Surv(time, event) ~", rhs)),
                     data = df2, ties = "efron")
        # gwasurvivr names the interaction term INTER.TERM and SNP main effect SNP
        row <- gw[gw$RSID == rsid, ]
        expect_equal(row[["COEF_INTER.TERM"]],
                     unname(coef(fit)["SNP:DrugTxYes"]), tolerance = 1e-5)
        expect_equal(row[["COEF_SNP"]],
                     unname(coef(fit)["SNP"]), tolerance = 1e-5)
        expect_equal(row[["SE.COEF_INTER.TERM"]],
                     unname(sqrt(diag(vcov(fit))["SNP:DrugTxYes"])),
                     tolerance = 1e-5)
    }
})
