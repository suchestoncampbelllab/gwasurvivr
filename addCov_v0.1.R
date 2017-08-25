addCov <- function(se, covFile){
        se <- se[,colnames(se) %in% covFile[[1]]]
        colData(se) <- merge(colData(se), covFile)
        colnames(se) <- colData(se)[[1]]
        return(se)
}


se.cov <- addCov(se.1, trm.cov)



