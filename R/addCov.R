#'
#'
#'
#'
#'

addCov <- function(se, covfile){
        covfile <- read.table(covfile, header=T, sep="\t", stringsAsFactors=F)
        se <- se[,colnames(se) %in% covfile[[1]]]
        colnames(covfile)[1] <- "ID_2"
        colData(se) <- DataFrame(plyr::join(data.frame(colData(se)), covfile))
        colnames(se) <- colData(se)$ID_2
        return(se)
}