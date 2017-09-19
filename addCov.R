addCov <- function(se, covFile){
        library(SummarizedExperiment)
        se <- se[,colnames(se) %in% covFile[[1]]]
        colnames(covFile)[1] <- "ID_2"
        colData(se) <- DataFrame(plyr::join(data.frame(colData(se)), covFile))
        colnames(se) <- colData(se)$ID_2
        return(se)
}