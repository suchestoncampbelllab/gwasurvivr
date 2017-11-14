
library(VariantAnnotation)
library(survival)
library(data.table)
library(dplyr)
library(tidyr)
library(broom)
library(microbenchmark)
library(parallel)

geno <- fread("~/GoogleDrive/Sucheston-Campbell Lab/survivR/data/genotype.dosage")
snps <- geno[[1]]
geno <- as.matrix(geno[,-1])
rownames(geno) <- snps

pheno <- read.table("~/GoogleDrive/Sucheston-Campbell Lab/survivR/data/pheno_file.txt", 
                    header=T, sep="\t", stringsAsFactors = T, row.names = "sample.ids")
pheno$distatD <- as.integer(pheno$distatD) -1


input.genotype <- geno[1,]

vcf.file="./chr21_sub/chr21.25000000-26000000.dose.vcf.recode.vcf.gz"
chunk.size=10000
time="intxsurv_1Y"
event="dead_1Y"
covariates=c("distatD", "age")
pheno.file = as.matrix(pheno)
sample.ids = paste0("SAMP", sample(1:1000, size=200))
output.name="test_survivR_chr21"



survFit <- function(input.genotype){}

# construct X, a matrix of covaritates and snp 
X <- cbind(
    input.genotype,
    pheno.file[,covariates]
)

# 

Y <- Surv(time=pheno.file[,time], event=pheno.file[,event])

rownames(Y) <- as.character(seq_len(nrow(Y)))

STRATA <- NULL

OFFSET <- rep(0, nrow(X))

INIT <- NULL

# check to see if control is interesting
# think of doing: 
# linear model -- you get a measure of goodness of fit after single iteration
## and first pass do filtering ... just a hand full of iterations
# and that will be enough info to determine if shuold continue on to do more expensive work
# sometimes statistic is figured out early iteration


# look into INIT 


# use predict as init 

INIT <- NULL


CONTROL <- structure(
    list(
        eps = 1e-09,
        toler.chol = 1.81898940354586e-12,
        iter.max = 20L,
        toler.inf = 3.16227766016838e-05,
        outer.max = 10L, 
        timefix = TRUE),
    .Names = c(
        "eps",
        "toler.chol",
        "iter.max", 
        "toler.inf", 
        "outer.max", 
        "timefix"
    )
)

WEIGHTS <- NULL
METHOD <- "efron"
ROWNAMES <- rownames(X)


## cox models find regression coefficents that together meaximize the log partial likelihood of the survival darta
# the inverse Hessian matrix (matrix of second-order partial derivatives of the log partial likelihood)
# evaluated at those estimated regression-coefficient values provides an estimate of the covariance-variance matrix among those regression coefficients
# therefore the stanrdard errors of the individual coefficients are the square root of the diaganol entries of the covariance matrix


# using coxph.fit with our created variables
fit <- coxph.fit(
    X, Y, STRATA, OFFSET, INIT, CONTROL, WEIGHTS, METHOD, ROWNAMES
)



summary(fit)
str(fit)
coef=fit$coefficients[1]
se=sqrt(diag(fit$var)[1])


# get summary data
(coef <- data.frame(
    coef=fit$coefficients[1],
    se=sqrt(diag(fit$var)[1])) %>%
    mutate(z=coef/se,
           pval=2*pnorm(abs(z), lower.tail = F),
           SNP=names(fit$coefficients[1])))


# compare with coxph function
surv.fun <- Surv(time = pheno$intxsurv_1Y, event=pheno$dead_1Y)
fit.org <- coxph(surv.fun ~ age + distatD, data=pheno)
summary(fit.org)


# system time for our coxph.fit
system.time({
    all <- apply(geno[1:5,], 1, function(x) {
        X <- cbind(x,X)
        fit <- coxph.fit(
            X, Y, STRATA, OFFSET, INIT, CONTROL, WEIGHTS, METHOD, ROWNAMES
        )
        (coef <- data.frame(
            coef=fit$coefficients[1],
            se=sqrt(diag(fit$var)[1])) %>%
                mutate(z=coef/se,
                       pval=2*pnorm(abs(z), lower.tail = F),
                       SNP=names(fit$coefficients[1])))
    }
    )        
})

# system time for coxph from survival package
system.time({
    all1 <- apply(geno[1:5,], 1, function(x) {
        fit.org <- coxph(Y ~ x + pheno$age + pheno$distatD)
        summary(fit.org)$coef[1]
    })     
})
