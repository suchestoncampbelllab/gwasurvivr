coxExtractInt <- function(cox.out, snp, n.sample, n.event){
    
    cox.cols <- colnames(cox.out)
    serrs <- cox.cols[grepl("serr", cox.cols)]
    coefs <- cox.cols[grepl("coef", cox.cols)]
    cols <- sapply(strsplit(coefs,"_"), "[[", 2)
    
    # calculate z-score
    z <- cox.out[,coefs]/cox.out[,serrs]
    colnames(z) <- paste("z", cols, sep="_")
    # calculate p-value
    pval <- apply(z, 2, function(z) 2*pnorm(abs(z), lower.tail=FALSE) )
    colnames(pval) <- paste("pvalue", cols, sep="_")
    # calculate hazard ratio
    hr <- exp(cox.out[,coefs])
    colnames(hr) <- paste("hr", cols, sep="_")
    # confidence interval HR
    lowerCI <- exp(cox.out[,coefs]-1.96*cox.out[,serrs])
    colnames(lowerCI) <- paste("lowerCI", cols, sep="_")
    upperCI <- exp(cox.out[,coefs]+1.96*cox.out[,serrs])
    colnames(upperCI) <- paste("upperCI", cols, sep="_")
    
    # putting everything back together
    sres <- cbind(cox.out, hr, lowerCI, upperCI, z, pval, n.sample, n.event)
    rownames(sres) <- NULL # remove rownames so we don't have a duplicated rownames issue
    
    data.frame(cbind(snp,sres))
}