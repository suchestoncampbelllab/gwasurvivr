# setwd("~/Google Drive/Sucheston-Campbell Lab/survivR/")
# 
# pdata <-fread("~/Google Drive/Sucheston-Campbell Lab/DBMT_PhenoData/DBMT_PhenoData_EA_long_lessVar20171107.txt")
# 
# 
# ## write function
# vcf.file="./data/example.vcf"
# chunk.size=10000
# time="intxsurv_1Y"
# event="dead_1Y"
# covariates=c("distatD", "age")
# sample.file = fread("data/example.sample")[-1,]
# pheno.file = pdata %>% 
#     mutate(sample.ids=paste0("sample", 1:nrow(.))) %>%
#     dplyr::select(sample.ids, intxsurv_1Y, age, distatD, intxsurv_1Y, dead_1Y)
# sample.ids = sample.file %$% 
#                 ID_2 %>%
#                 sample(1000, replace = F)
# 
# output.name="test_survivR/chr21"
# 

devtools::install_github("suchestoncampbelllab/survivR")



pdata <-fread("~/Google Drive/Sucheston-Campbell Lab/DBMT_PhenoData/DBMTpheno_EA_long_20171023.txt")
vcf.file="../chr21_sub/chr21.25000000-26000000.dose.vcf.recode.vcf.gz"
chunk.size=10000
time="intxsurv_1Y"
event="dead_1Y"
covariates=c("distatD", "age")
pheno.file = pdata %>% 
    mutate(sample.ids=paste0("SAMP", 1:nrow(.))) %>%
    dplyr::select(sample.ids, intxsurv_1Y, age, distatD, intxsurv_1Y, dead_1Y)
sample.ids = paste0("SAMP", sample(1:1000, size=200))
output.name="chr21"
info.file = "../chr21.info"


geno <- fread("~/Google Drive/Sucheston-Campbell Lab/survivR/data/genotype.dosage")
pheno <- read.table("~/Google Drive/Sucheston-Campbell Lab/survivR/data/pheno_file.txt", header=T, sep="\t")



coxph.fit

# construct x -- 
X <- cbind(
    # as.numeric(geno[1, -1]), # only thing that is changing 
    pheno$age,
    as.integer(as.factor(pheno$distatD)) - 1L #,
    # pheno$age*(as.integer(pheno$distatD) - 1L)
)

dimnames(X) <- list(seq_len(nrow(X)), c( "age", "distatDearly"))#, "a:dist"))

# 
Y <- surv.fun

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
ROWNAMES <- as.character(seq_len(nrow(X)))


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

# get summary data
(coef <- data.frame(
    coef=fit$coefficients,
    se=sqrt(diag(fit$var))) %>%
    mutate(z=coef/se,
           pval=2*pnorm(abs(z), lower.tail = F)))


# compare with coxph function
surv.fun <- Surv(time = pheno$intxsurv_1Y, event=pheno$dead_1Y)
fit.org <- coxph(surv.fun ~ age + distatD, data=pheno)
summary(fit.org)


# system time for our coxph.fit
system.time({
    all <- apply(geno[1:5,-1], 1, function(x) {
        X <- cbind(X,x)
        fit <- coxph.fit(
            X, Y, STRATA, OFFSET, INIT, CONTROL, WEIGHTS, METHOD, ROWNAMES
        )
        (coef <- data.frame(
            coef=fit$coefficients,
            se=sqrt(diag(fit$var))) %>%
                mutate(z=coef/se,
                       pval=2*pnorm(abs(z), lower.tail = F)))
    })        
})

# system time for coxph from survival package
system.time({
    all1 <- apply(geno[1:5,-1], 1, function(x) {
        fit.org <- coxph(surv.fun ~ x + pheno$age + pheno$distatD)
        summary(fit.org)$coef
    })     
})
