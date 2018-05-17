coxExtract <- function(cox.out, snp, n.sample, n.event, print.covs="only"){
    if(print.covs == "only"){
        # calculate z-score
        z <- cox.out[,1]/cox.out[,2]
        # calculate p-value
        pval <- 2*pnorm(abs(z), lower.tail=FALSE)
        # calculate hazard ratio
        hr <- exp(cox.out[,1])
        # confidence interval HR
        lowerCI <- exp(cox.out[,1]-1.96*cox.out[,2])
        upperCI <- exp(cox.out[,1]+1.96*cox.out[,2])
        # putting everything back together
        sres <- cbind(pval, hr, lowerCI, upperCI, cox.out, z, n.sample, n.event)
        colnames(sres) <- c("PVALUE",
                            "HR",
                            "HR_lowerCI", 
                            "HR_upperCI",
                            "COEF",
                            "SE.COEF",
                            "Z",
                            "N", 
                            "NEVENT")
        rownames(sres) <- NULL 
    } else {
        cox.cols <- colnames(cox.out)
        coefs <- cox.cols[grepl("COEF", cox.cols)]
        cols <- sapply(strsplit(coefs, "_"), "[[", 2)
        colnames(cox.out) <- gsub("SERR","SE.COEF", cox.cols)
        serrs <- colnames(cox.out)[grepl("SE.COEF", colnames(cox.out))]
        # calculate z-score
        z <- cox.out[, coefs] / cox.out[, serrs]
        # calculate p-value
        pvalFun <- function(z) 2 * pnorm(abs(z), lower.tail=FALSE)
        if(is.matrix(z)){
            colnames(z) <- paste("Z", cols, sep = "_")
            pval <- apply(z, 2, pvalFun)
            colnames(pval) <- paste("PVALUE", cols, sep = "_")
            # calculate hazard ratio
            hr <- exp(cox.out[, coefs])
            colnames(hr) <- paste("HR", cols, sep = "_")
            # confidence interval HR
            lowerCI <- exp(cox.out[, coefs] - 1.96 * cox.out[, serrs])
            colnames(lowerCI) <- paste("HR_lowerCI", cols, sep = "_")
            upperCI <- exp(cox.out[, coefs] + 1.96 * cox.out[, serrs])
            colnames(upperCI) <- paste("HR_upperCI", cols, sep = "_")
        } else {
            z <- setNames(z, paste("Z", cols, sep = "_"))
            pval <- pvalFun(z)
            pval <- setNames(pval, paste("PVALUE", cols, sep = "_"))
            # calculate hazard ratio
            hr <- exp(cox.out[, coefs])
            hr <- setNames(hr,paste("HR", cols, sep = "_"))
            # confidence interval HR
            lowerCI <- exp(cox.out[, coefs] - 1.96 * cox.out[, serrs])
            lowerCI <- setNames(lowerCI, paste("HR_lowerCI", cols, sep = "_"))
            upperCI <- exp(cox.out[, coefs] + 1.96 * cox.out[, serrs])
            upperCI <- setNames(upperCI, paste("HR_upperCI", cols, sep = "_"))
        }
        sres <- cbind(pval,
                      hr,
                      lowerCI,
                      upperCI,
                      cox.out,
                      z,
                      N=n.sample,
                      NEVENT=n.event)
        rownames(sres) <- NULL
    }
    data.frame(cbind(snp,sres))
}
