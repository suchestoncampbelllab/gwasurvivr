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