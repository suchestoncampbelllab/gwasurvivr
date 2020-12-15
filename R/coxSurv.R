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

loadProcessWrite <- function (x, ...) {
  UseMethod("loadProcessWrite", x)
}


