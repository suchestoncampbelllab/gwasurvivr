messageCompressionTime <- function(comp_time) {
  message("***** Compression time ******\n", 
          "User:", round(comp_time[1], 3),
          "\nSystem: ", round(comp_time[2], 3),
          "\nElapsed: ", round(comp_time[3], 3), 
          "\n*****************************"
  )
}