coxSurv <- function(cox_surv) {
  # browser()
  if(cox_surv$verbose) message("Analysis started on ",
                      format(Sys.time(), "%Y-%m-%d"),
                      " at ",
                      format(Sys.time(), "%H:%M:%S"))

  ############################################################################
  #### Phenotype data wrangling ##############################################
  cox.params <- coxPheno(pheno.file = cox_surv$covariate.file, 
                         covariates = cox_surv$covariates, 
                         id.column = cox_surv$id.column, 
                         inter.term = cox_surv$inter.term, 
                         time.to.event = cox_surv$time.to.event,
                         event = cox_surv$event, 
                         sample.ids = cox_surv$sample.ids, 
                         verbose = cox_surv$verbose)

  ############################################################################
  
  ############################################################################
  ##### Generate cluster obj #################################################
  cl <- create_cluster_obj(cox_surv$clusterObj)
  on.exit(stopCluster(cl), add=TRUE)
  ############################################################################
  
  results <- loadProcessWrite(x = cox_surv, cox.params = cox.params, cl = cl)
  
  if(cox_surv$verbose) closing_messages(snps_removed = results$snps_removed,
                               snps_analyzed = results$snps_analyzed,
                               out.file = cox_surv$out.file)
  
}



#create cluster object depending on user pref or OS type,
# also create option to input number of cores
create_cluster_obj <- function(clusterObj) {
  if(!is.null(clusterObj)){
    cl <- clusterObj
  }else if(.Platform$OS.type == "unix") {
    cl <- makePSOCKcluster(getOption("gwasurvivr.cores", 2L))
  } else {
    cl <- makeCluster(getOption("gwasurvivr.cores", 2L))
  }
  return(cl)
}



closing_messages <- function(snps_removed, snps_analyzed, out.file){
  
  message("Analysis completed on ",
          format(Sys.time(), "%Y-%m-%d"),
          " at ",
          format(Sys.time(), "%H:%M:%S"))
  
  message(snps_removed, " SNPs were removed from the analysis for ",
          "not meeting the threshold criteria.")
  
  message("List of removed SNPs can be found in ", 
          paste0(out.file, ".snps_removed"))
  
  message(snps_analyzed, " SNPs were analyzed in total")
  
  message("The survival output can be found at ",
          paste0(out.file, ".coxph"))
}



loadProcessWrite <- function (x, ...) {
  UseMethod("loadProcessWrite", x)
}

processSNPGenotypes <- function (x, ...) {
  UseMethod("processSNPGenotypes", x)
}

getGenoData <- function (x, ...) {
  UseMethod("getGenoData", x)
}


