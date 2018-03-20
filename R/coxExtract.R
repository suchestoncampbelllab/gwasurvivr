#' @export 

coxExtract <- function(cox.out, snp, n.sample, n.event){
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
    sres <- cbind(cox.out, hr, lowerCI, upperCI, z, pval, n.sample, n.event)
    colnames(sres) <- c("coef", "se.coef", "exp.coef", "lb", "ub", "z", "p.value", "n", "nevents")
    rownames(sres) <- NULL # remove rownames so we don't have a duplicated rownames issue
    
    data.frame(cbind(snp,sres))
}



