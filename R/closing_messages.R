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